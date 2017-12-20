#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "args_parser.h"
#include <mpi.h>

struct RowsSplitter{
    float *subA;
    float *subtemp;
    int n;
    int ri;
    int rf;
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

float laplace_step(float *in, float *out, int n, int rank, int size)
{
    int i, j;
    float error=0.0f;
    for (j = n*rank/size; j < n*(rank + 1)/size; j++)
        for (i=0; i < n; i++)
        {
          out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
          error = max_error( error, out[j*n+i], in[j*n+i] );
        }
    return error;
}

void laplace_init(struct RowsSplitter submatrix)
{
    int i, n = submatrix.n;
    float V;
    const float pi  = 2.0f * asinf(1.0f);

    memset(submatrix.subA, 0, n * sizeof(float) * (submatrix.rf - submatrix.ri) + 2 * n);
    memset(submatrix.subtemp, 0, n * sizeof(float) * (submatrix.rf - submatrix.ri) + 2 * n);

    for (i = 1; i < (submatrix.rf - submatrix.ri) + 1; i++)
    {
        V = sinf(pi * (i + submatrix.ri - 1) / (n - 1));
        submatrix.subA[i * n] = V;
        submatrix.subtemp[i * n] = V;

        submatrix.subA[i * n + n - 1] = V * expf(-pi);
        submatrix.subtemp[i * n + n - 1] = V * expf(-pi);
    }
}

void print_matrix(struct RowsSplitter submatrix, char file_dir[], int rank)
{
    int i, j;
    sprintf(file_dir, "%s_%d", file_dir, rank);
    FILE *f = fopen(file_dir, "w");

    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    for(i = 0; i < (submatrix.rf - submatrix.ri) + 2 ; i++)
    {
        for(j = 0; j < submatrix.n; j++)
            fprintf(f, "%f\t", submatrix.subA[i * submatrix.n + j]);
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
    int ri = n * rank / size;
    int rf = n * (rank + 1) / size;

    // Allocate memory to matrices and arrays
    float *subA    = (float *) malloc(n * sizeof(float) * (rf - ri) + 2 * n * sizeof(float));
    float *subtemp = (float *) malloc(n * sizeof(float) * (rf - ri) + 2 * n * sizeof(float));

    // Initiate RowsSplitter struct
    struct RowsSplitter submatrix = {subA, subtemp, n, ri, rf};

    printf("n: %d, ri: %d, rf: %d\n", n, ri, rf);

    // Set boundary conditions
    laplace_init(submatrix);

    printf("Initialized!\n");

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
