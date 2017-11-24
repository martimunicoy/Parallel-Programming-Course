#include <stdio.h>
#include <omp.h>

int main()
{
    int t, p;
    t = omp_get_num_threads();
    p = omp_get_num_procs();
    printf("Number of threads = %d\n", t);
    printf("Number of processors = %d\n", p);
    #pragma omp parallel
    	printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}
