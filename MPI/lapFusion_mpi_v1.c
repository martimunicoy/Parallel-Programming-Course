#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "args_parser.h"
#include <mpi.h>

#define subA(i,j) subA[i * submatrix.n + j]
#define subtemp(i,j) subtemp[i * submatrix.n + j]

struct Submatrix {
    float *subA;
    float *subtemp;
    int n;
    int rank;
    int size;
    int ri;
    int rf;
    int sub_n;
};

float stencil ( float v1, float v2, float v3, float v4)
{
    return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
    float t= fabsf( new - old );
    return t>prev_error? t: prev_error;
}

//¿esta entendiendo bien el sen recv el buffer con esa structura? ¿¿¿¿¿¿¿????????
void update_rows(struct Submatrix submatrix)
{
    if (submatrix.rank == 0) /* First submatrix */
    {
        MPI_Send(&submatrix.subA(submatrix.sub_n - 1,0), submatrix.n,
                 MPI_FLOAT, submatrix.rank + 1, submatrix.rank + 1, MPI_COMM_WORLD);
        MPI_Recv(&submatrix.subA(submatrix.sub_n, 0), submatrix.n,
                 MPI_FLOAT, submatrix.rank + 1, submatrix.rank, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
    else if (submatrix.rank == submatrix.size - 1) /* Last submatrix */
    {
        MPI_Send(&submatrix.subA(1,0), submatrix.n, MPI_FLOAT,
                 submatrix.rank - 1, submatrix.rank - 1, MPI_COMM_WORLD);
        MPI_Recv(&submatrix.subA(0,0), submatrix.n, MPI_FLOAT,
                 submatrix.rank - 1, submatrix.rank, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
    else /* Inner submatrix */
    {
        MPI_Send(&submatrix.subA(submatrix.sub_n - 1,0), submatrix.n,
                 MPI_FLOAT, submatrix.rank + 1, submatrix.rank + 1, MPI_COMM_WORLD);
        MPI_Send(&submatrix.subA(1,0), submatrix.n, MPI_FLOAT,
                 submatrix.rank - 1, submatrix.rank - 1, MPI_COMM_WORLD);
        MPI_Recv(&submatrix.subA(0,0), submatrix.n, MPI_FLOAT,
                 submatrix.rank - 1, submatrix.rank, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Recv(&submatrix.subA(submatrix.sub_n, 0),
                 submatrix.n, MPI_FLOAT, submatrix.rank + 1, submatrix.rank, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }
}

float laplace_step(struct Submatrix  submatrix)
{
    int i, j;
    float error=0.0f;
    /*
    for (j = n*rank/size; j < n*(rank + 1)/size; j++)
        for (i=0; i < n; i++)
        {
          out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
          error = max_error( error, out[j*n+i], in[j*n+i] );
        }
    */
    return error;
}


void laplace_init(struct Submatrix  submatrix)
{
    int i, row;
    int n = submatrix.n;
    float V;
    const float pi  = 2.0f * asinf(1.0f);

    //Fill with zeros the allocated memory reserved for the matrix
    memset(submatrix.subA, 0, n * sizeof(float) * (submatrix.rf - submatrix.ri + 2));
    memset(submatrix.subtemp, 0, n * sizeof(float) * (submatrix.rf - submatrix.ri + 2));

    row = 1;
    for (i = submatrix.ri; i < submatrix.rf; ++i)
    {

        V = sinf(pi * i / (n - 1));
        submatrix.subA(row, 0) = V;
        submatrix.subtemp(row, 0) = V;

        submatrix.subA(row, n - 1) = V * expf(-pi);
        submatrix.subtemp(row, n - 1) = V * expf(-pi);
        row++;
    }
}

void print_matrix(struct Submatrix  submatrix, char file_dir[], int rank)
{
    int i, j;
    char *file_dir_temp = strdup(file_dir);
    char *filename = strsep(&file_dir_temp, ".");
    char *extension = strsep(&file_dir_temp, ".");
    sprintf(file_dir, "%s_%d.%s", filename, rank, extension);
    FILE *f = fopen(file_dir, "w");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for(i = 0; i < (submatrix.rf - submatrix.ri) + 2 ; i++)
    {
        for(j = 0; j < submatrix.n; j++)
            fprintf(f, "%f\t", submatrix.subA(i,j));
        fprintf(f, "\n");
    }

    fclose(f);
}

int main(int argc, char** argv)
{
    float *A, *temp;
    struct timeval t0, t1;

    // get runtime arguments
    struct Args parsed_args = args_parser(argc, argv);
    int n = parsed_args.n;
    int iter_max = parsed_args.iter_max;
    bool out = parsed_args.out;
    char *file_dir = malloc(50);
    strcpy(file_dir, parsed_args.file_dir);
    bool count = parsed_args.count;

    if (count)
        gettimeofday(&t0, 0);

    const float tol = 1.0e-5f;
    float error= 1.0f;

    int iter = 0;

    // Initiate MPI region
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Determine initial and final rows indexes
    int sub_n = n / size;
    int ri = rank * sub_n;
    int rf = (rank + 1) * sub_n;

    // Allocate memory to matrices and arrays
    float *subA    = (float *) malloc(n * sizeof(float) * (rf - ri + 2));
    float *subtemp = (float *) malloc(n * sizeof(float) * (rf - ri + 2));

    // Initiate Submatrix  struct
    struct Submatrix  submatrix = {subA, subtemp, n, rank, size, ri, rf, sub_n};

    printf("n: %d, ri: %d, rf: %d\n", n, ri, rf);

    // Set boundary conditions
    laplace_init(submatrix);

    printf("Initialized!\n");

    update_rows(submatrix);

    print_matrix(submatrix, parsed_args.file_dir, rank);

    /*
    // set singular point
    A[(n/128)*n+n/128] = 1.0f;

    if (rank == 0)
      printf("Jacobi relaxation Calculation: %d x %d mesh,"
           " maximum of %d iterations, using %d processors\n",
           n, n, iter_max, size);


    while ( error > tol*tol && iter < iter_max )
    {
      iter++;
      error= laplace_step (A, temp, n, rank, size);

      // swap pointers A & temp
      float *swap= A; A=temp; temp= swap;
    }
    */

    // End MPI region
    MPI_Finalize();

    /*
    error = sqrtf( error );

    if (out)
      print_matrix(A, n, file_dir);

    printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
    printf("A[%d][%d]= %0.6f", n/128, n/128, A[(n/128)*n+n/128]);

    if (count)
    {
      gettimeofday(&t1, 0);
      long int diff = (t1.tv_sec - t0.tv_sec) *1000000L + t1.tv_usec - t0.tv_usec;
      printf(", Running time: %ld", diff);
    }
    printf("\n");
    */
    free(file_dir);
    free(A); free(temp);
}
