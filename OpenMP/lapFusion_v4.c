#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <stdbool.h>
#include "args_parser.h"

int iter;
float error;
#pragma omp threadprivate(iter,error)

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n)
{
  int i, j;
  float local_error = 0.0f;
  #pragma omp for
  for ( j=1; j < n-1; j++ )
  {
    #pragma omp simd
    for ( i=0; i < n; i++ )
    {
      out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i],
                          in[(j+1)*n+i]);
      local_error = max_error( local_error, out[j*n+i], in[j*n+i] );
    }
  }
  return local_error;
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
    for(j = 0; j < n-1; j++) fprintf(f, "%f\t", matrix[i*n+j]);
    fprintf(f, "%f", matrix[i*n+j+1]);
    fprintf(f, "\n");
  }

  fclose(f);
}

int main(int argc, char* argv[])
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
  bool count = parsed_args.count;

  if (count) gettimeofday(&t0, 0);

  const float tol = 1.0e-5f;

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

  iter = 0;
  error = 1.0f;
  float errors_list[num_threads];

  #pragma omp parallel firstprivate(A,temp) copyin(iter,error)
  {
    while ( (error > tol*tol || iter == 0) && iter < iter_max )
    {
      iter++;
      error = laplace_step (A, temp, n);
      float *swap= A; A=temp; temp= swap; // swap pointers A & temp
    }
    errors_list[omp_get_thread_num()] = error;
  }

  int i;
  for(i=0; i<num_threads; i++)
    if(errors_list[i] > error) error=errors_list[i];
  error = sqrtf( error );

  if (out)
    print_matrix(A, n, file_dir);

  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  printf("A[%d][%d]= %0.6f", n/128, n/128, A[(n/128)*n+n/128]);

  if (count)
  {
    gettimeofday(&t1, 0);
    long int diff = (t1.tv_sec - t0.tv_sec) *1000000L + t1.tv_usec - t0.tv_usec;
    printf(", Running time: %ld\n", diff);
  }
  else printf("\n");

  free(file_dir);
  free(A); free(temp);
}
