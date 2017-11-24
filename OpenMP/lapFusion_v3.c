#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>
#include <stdbool.h>
#include <string.h>

const int DEFAULT_N = 4096;
const int DEFAULT_ITER_MAX = 1000;
const int ARGS_NUM = 3;
const char ARGS[ARGS_NUM][10] = {"-N", "-I", "-T"};

struct Args{
  int n;
  int iter_max;
  int num_threads;
};

bool startsWith(const char *string, const char *prefix){
  return strncmp(string, prefix, strlen(prefix)) == 0;
}

struct Args args_parser(int argc, char *argv[])
{
  int n = DEFAULT_N;
  int iter_max = DEFAULT_ITER_MAX;
  int num_threads = omp_get_num_procs();
  for (int i = 1; i < argc; i++){
    for (int j = 0; j < ARGS_NUM; j++){
      if (startsWith(*(argv + i), ARGS[j])){
        if (i + 2 <= argc){
          if (j == 0) n = atoi(*(argv + i + 1));
          if (j == 1) iter_max = atoi(*(argv + i + 1));
          if (j == 2) num_threads = atoi(*(argv + i + 1));
        }
        else{
          printf("Wrong command, using default values.\n");
        }
        break;
      }
    }
  }
  if (n<1 | n>DEFAULT_N){
    printf("\'N\' out of range, using default value (%d)\n", DEFAULT_N);
    n = DEFAULT_N;
  }
  if (iter_max<1 | iter_max>999999){
    printf("\'iter_max\' out of range, using default value (%d)\n", DEFAULT_ITER_MAX);
    iter_max = DEFAULT_ITER_MAX;
  }
  if (num_threads<1 | num_threads>8){
    printf("\'num_threads\' out of range, using default value (%d)\n", omp_get_num_procs());
    num_threads = omp_get_num_procs();
  }
  struct Args parsed_args = {n, iter_max, num_threads};
  return parsed_args;
}


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
  int i, j, id, jstart, jend;
  float error=0.0f;
  #pragma omp parallel reduction(max:error) private(i, j, id, jstart, jend)
  {
    id = omp_get_thread_num();
    jstart = id * (n-2) / num_threads + 1;
    jend = (id+1) * (n-2) / num_threads + 1;
    for ( j=jstart; j < jend; j++ )
    {
      for ( i=1; i < (n-1); i++ )
      {
        out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
        error = max_error( error, out[j*n+i], in[j*n+i] );
      }
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

int main(int argc, char* argv[])
{
  float *A, *temp;
  struct timeval t0, t1;

  // get runtime arguments
  struct Args parsed_args = args_parser(argc, argv);
  int n = parsed_args.n;
  int iter_max = parsed_args.iter_max;
  int num_threads = parsed_args.num_threads;

  gettimeofday(&t0, 0);

  const float tol = 1.0e-5f;
  float error= 1.0f;

  A    = (float*) malloc( n*n*sizeof(float) );
  temp = (float*) malloc( n*n*sizeof(float) );

  //  set boundary conditions
  laplace_init (A, n);
  laplace_init (temp, n);

  /*
  int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      printf("%f\t", A[i*n+j]);
    }
    printf("\n");
  }
  */

  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations, using %d threads\n",
         n, n, iter_max, num_threads);

  int iter = 0;

  omp_set_num_threads(num_threads);

  while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    error= laplace_step (A, temp, n, num_threads);
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
  }
  error = sqrtf( error );

  /*
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      printf("%f\t", A[i*n+j]);
    }
    printf("\n");
  }
  */

  gettimeofday(&t1, 0);

  long int diff = (t1.tv_sec - t0.tv_sec) *1000000L + t1.tv_usec - t0.tv_usec;

  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  printf("A[%d][%d]= %0.6f, Running time: %ld\n", n/128, n/128, A[(n/128)*n+n/128], diff);

  free(A); free(temp);
}
