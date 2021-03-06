#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "args_parser.h"
#include <mpi.h>

const int MAXIMUM_FLOATS_PER_MESSAGE = 1000;

struct Submatrix{
    float *subA;
    float *subtemp;
    int n;
    int rank;
    int size;
    int ri;
    int rf;
    int subrows;
};

float stencil(float v1, float v2, float v3, float v4)
{
    return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error(float prev_error, float old, float new)
{
    float t= fabsf( new - old );
    return t > prev_error? t: prev_error;
}

void update_rows(struct Submatrix sm)
{
    if (sm.size > 1)
    {
        int previous    = sm.rank - 1;
        int next        = sm.rank + 1;
        int fst_row     = 0;
        int snd_row     = sm.n;
        int sndlast_row = sm.subrows * sm.n;
        int last_row    = (sm.subrows + 1) * sm.n;
        int floats_left = sm.n;
        int floats_to_send;

        // Split large messages into smaller ones to prevent MPI from hanging
        while (floats_left != 0)
        {
            if (floats_left > MAXIMUM_FLOATS_PER_MESSAGE)
                floats_to_send = MAXIMUM_FLOATS_PER_MESSAGE;
            else
                floats_to_send = floats_left;

            // First submatrix
            if (sm.rank == 0)
            {
                MPI_Send(&sm.subA[sndlast_row], floats_to_send, MPI_FLOAT,
                         next, 0, MPI_COMM_WORLD);
                MPI_Recv(&sm.subA[last_row], floats_to_send, MPI_FLOAT, next,
                         0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Last submatrix
            else if (sm.rank == sm.size - 1)
            {
                MPI_Send(&sm.subA[snd_row], floats_to_send, MPI_FLOAT,
                         previous, 0, MPI_COMM_WORLD);
                MPI_Recv(&sm.subA[fst_row], floats_to_send, MPI_FLOAT,
                         previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // Inner submatrix
            else
            {
                MPI_Send(&sm.subA[snd_row], floats_to_send, MPI_FLOAT,
                         previous, 0, MPI_COMM_WORLD);
                MPI_Send(&sm.subA[sndlast_row], floats_to_send, MPI_FLOAT,
                         next, 0, MPI_COMM_WORLD);
                MPI_Recv(&sm.subA[fst_row], floats_to_send, MPI_FLOAT,
                         previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&sm.subA[last_row], floats_to_send, MPI_FLOAT, next,
                         0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            floats_left -= floats_to_send;
        }
    }
}

float laplace_step(struct Submatrix sm)
{
    int i, j;
    float error=0.0f;

    int initial_row = 1;
    int final_row = sm.subrows + 1;
    if (sm.rank == 0)
        ++initial_row;
    if (sm.rank == sm.size - 1)
        --final_row;

    for (j = initial_row; j < final_row; j++)
        for (i = 1; i < sm.n - 1; i++)
        {
          sm.subtemp[j * sm.n + i]= stencil(sm.subA[j * sm.n + i + 1],
                                            sm.subA[j * sm.n + i - 1],
                                            sm.subA[(j - 1) * sm.n + i],
                                            sm.subA[(j + 1) * sm.n + i]);
          error = max_error(error, sm.subtemp[j * sm.n + i],
                            sm.subA[j * sm.n + i]);
        }
    return error;
}

void laplace_init(struct Submatrix sm)
{
    int i, n = sm.n;
    float V;
    const float pi  = 2.0f * asinf(1.0f);

    memset(sm.subA, 0, n * sizeof(float) * (sm.rf - sm.ri + 2));
    memset(sm.subtemp, 0, n * sizeof(float) * (sm.rf - sm.ri + 2));

    for (i = 1; i < sm.subrows + 1; i++)
    {
        V = sinf(pi * (i + sm.ri - 1) / (n - 1));
        sm.subA[i * n] = V;
        sm.subtemp[i * n] = V;

        sm.subA[i * n + n - 1] = V * expf(-pi);
        sm.subtemp[i * n + n - 1] = V * expf(-pi);
    }
}

float get_global_error(int rank, int size, float maximum_error)
{
    int i;
    if (rank == 0)
    {
        for (i = 1; i < size; i++)
        {
            float local_error;
            MPI_Recv(&local_error, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            if (local_error > maximum_error)
                maximum_error = local_error;
        }
        for (i = 1; i < size; i++)
            MPI_Send(&maximum_error, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(&maximum_error, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        MPI_Recv(&maximum_error, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
    return maximum_error;
}

void send_submatrix(struct Submatrix sm)
{
    MPI_Send(&sm.ri, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&sm.subrows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&sm.subA[sm.n], sm.subrows * sm.n, MPI_FLOAT, 0, 0,
             MPI_COMM_WORLD);
}

void recv_submatrix(int origin, int n, float *A)
{
    int ri;
    int subrows;
    MPI_Recv(&ri, 1, MPI_INT, origin, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&subrows, 1, MPI_INT, origin, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Recv(&A[ri * n], subrows * n, MPI_FLOAT, origin, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
}

void print_matrix(struct Submatrix sm, char file_dir[])
{
    int i, j;
    char *file_dir_temp = strdup(file_dir);
    char *filename      = strsep(&file_dir_temp, ".");
    char *extension     = strsep(&file_dir_temp, ".");
    char file[50];
    sprintf(file, "%s_%d.%s", filename, sm.rank, extension);
    FILE *f             = fopen(file, "w");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for(i = 0; i < sm.subrows + 2; i++)
    {
        for(j = 0; j < sm.n; j++)
            fprintf(f, "%f\t", sm.subA[i * sm.n + j]);
        fprintf(f, "\n");
    }

    fclose(f);
}

void plot_matrix(float *matrix, int n, char file_dir[], char sufix[])
{
    int i, j;
    char *file_dir_temp = strdup(file_dir);
    char *filename      = strsep(&file_dir_temp, ".");
    char *extension     = strsep(&file_dir_temp, ".");
    char file[50];
    sprintf(file, "%s_%splot.%s", filename, sufix, extension);
    FILE *f             = fopen(file, "w");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f, "library(plotly)\n");

    for(i = 0; i < n; i++)
    {
        fprintf(f, "z%d <- c(", i + 1);
        for(j = 0; j < n - 1; j++)
            fprintf(f, "%f,", matrix[i * n + j]);
        fprintf(f, "%f)\n", matrix[i * n + n - 1]);
    }

    fprintf(f, "z <- c(");
    for(i = 0; i < n-1; i++)
        fprintf(f, "z%d,", i + 1);
    fprintf(f, "z%d)\n", n);

    fprintf(f, "dim(z) <- c(%d,%d)\n", n, n);
    fprintf(f, "p <- plot_ly(showscale = FALSE) %%>%% add_surface(z)\np\n");

    fclose(f);
}

int main(int argc, char** argv)
{
    int i;
    struct timeval t0, t1;

    // get runtime arguments
    struct Args parsed_args = args_parser(argc, argv);
    int n                   = parsed_args.n;
    int iter_max            = parsed_args.iter_max;
    bool out                = parsed_args.out;
    char *file_dir          = malloc(50);
    strcpy(file_dir, parsed_args.file_dir);
    bool count              = parsed_args.count;

    if (count)
        gettimeofday(&t0, 0);

    // Initiate constants
    const float tol = 1.0e-5f;

    // Initiate variables
    float global_error = 1.0f;
    float local_error  = global_error;
    int iter           = 0;

    // Allocate memory for the results
    float *A = (float*) malloc(n * n * sizeof(float));

    // Initiate MPI region
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Determine initial and final rows indexes
    int ri      = n * rank / size;
    int rf      = n * (rank + 1) / size;
    int subrows = rf - ri;

    // Allocate memory to matrices and arrays
    float *subA    = (float*) malloc(n * sizeof(float) * (subrows + 2));
    float *subtemp = (float*) malloc(n * sizeof(float) * (subrows + 2));
    float *swap;

    // Initiate Submatrix struct
    struct Submatrix submatrix = {subA, subtemp, n, rank, size, ri, rf,
                                  subrows};

    // Set boundary conditions
    laplace_init(submatrix);

    // Set singular point
    for (i = 0; i < size; i++)
        if (n / 128 > ri && n / 128 < rf)
            submatrix.subA[(n / 128 + 1) * n + n / 128] = 1.0f;

    // Display configuration
    if (rank == 0)
        printf("Jacobi relaxation Calculation: %d x %d mesh,"
               " maximum of %d iterations, using %d processes\n",
               n, n, iter_max, size);

    // Print initial matrix
    if (out)
    {
        // Join submatrices
        if (rank == 0)
        {
            int i;
            for (i = 0; i < subrows * n; i++)
                A[i] = submatrix.subA[n + i];
            for (i = 1; i < size; i++)
                recv_submatrix(i, n, A);
            plot_matrix(A, n, parsed_args.file_dir, "initial");
        }
        else
            send_submatrix(submatrix);
    }

    // Initiate propagation
    while (global_error > tol * tol && iter < iter_max)
    {
        // Perform an iteration and calculate the associated error
        iter++;
        local_error = laplace_step(submatrix);
        global_error = get_global_error(rank, size, local_error);

        // Swap pointers A & temp
        swap = submatrix.subA;
        submatrix.subA = submatrix.subtemp;
        submatrix.subtemp = swap;

        // Update rows among processes
        update_rows(submatrix);
    }

    // Join submatrices
    if (rank == 0)
    {
        for (i = 0; i < subrows * n; i++)
            A[i] = submatrix.subA[n + i];
        for (i = 1; i < size; i++)
            recv_submatrix(i, n, A);
    }
    else
        send_submatrix(submatrix);

    // Save results
    if (out)
    {
            print_matrix(submatrix, parsed_args.file_dir);
            plot_matrix(A, n, parsed_args.file_dir, "final");
    }

    // End MPI region
    MPI_Finalize();

    // Print results
    if (rank == 0)
    {
        global_error = sqrtf(global_error);
        printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, global_error);
        printf("A[%d][%d]= %0.6f", n/128, n/128, A[(n/128)*n+n/128]);
        if (count)
        {
            gettimeofday(&t1, 0);
            long int diff = (t1.tv_sec - t0.tv_sec) *1000000L + t1.tv_usec -
                            t0.tv_usec;
            printf(", Running time: %ld", diff);
        }
        printf("\n");
    }

    free(file_dir);
    free(subA); free(subtemp);
    free(A);
}
