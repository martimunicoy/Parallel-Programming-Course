#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include "args_parser.h"

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n, int num_threads)
{
  int j;
  float error=0.0f;
  for ( j=1; j < n-1; j++ )
    #pragma omp parallel reduction(max:error)
    {
      int i, id, istart, iend;
      id = omp_get_thread_num();
      istart = id * n / num_threads;
      iend = (id+1) * n / num_threads;
      for ( i=istart; i < iend; i++ )
      {
        out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
        error = max_error( error, out[j*n+i], in[j*n+i] );
      }
    }
  return error;
}

void laplace_init(float *in, int n)
{
  int i;
  float V;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in, 0, n*n*sizeof(float));
  for (i=0; i<n; i++) {
    V = in[i*n] = sinf(pi*i / (n-1));
    in[i*n+n-1] = V * expf(-pi);
  }
}

void print_matrix(float *matrix, int n, char file_dir[50])
{
  int i, j;
  FILE *f = fopen(file_dir, "w");

  if (f == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

  for(i = 0; i < n; i++)
  {
    for(j = 0; j < n; j++) fprintf(f, "%f\t", matrix[i*n+j]);
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
  int num_threads = parsed_args.num_threads;
  bool out = parsed_args.out;
  char *file_dir = malloc(50);
  strcpy(file_dir, parsed_args.file_dir);

  gettimeofday(&t0, 0);

  const float tol = 1.0e-5f;
  float error= 1.0f;

  A    = (float*) malloc( n*n*sizeof(float) );
  temp = (float*) malloc( n*n*sizeof(float) );

  //  set boundary conditions
  laplace_init (A, n);
  laplace_init (temp, n);

  A[(n/128)*n+n/128] = 1.0f; // set singular point

  omp_set_num_threads(num_threads);

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations, using %d threads\n",
         n, n, iter_max, num_threads);

  int iter = 0;
  while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    error= laplace_step (A, temp, n, num_threads);
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
  }

  error = sqrtf( error );

  if (out)
    print_matrix(A, n, file_dir);

  gettimeofday(&t1, 0);

  long int diff = (t1.tv_sec - t0.tv_sec) *1000000L + t1.tv_usec - t0.tv_usec;

  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  printf("A[%d][%d]= %0.6f, Running time: %ld\n", n/128, n/128, A[(n/128)*n+n/128], diff);

  free(A); free(temp);
}
