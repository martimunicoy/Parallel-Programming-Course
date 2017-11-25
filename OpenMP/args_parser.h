#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#define ARGS_NUM 5

const int DEFAULT_N = 4096;
const int DEFAULT_ITER_MAX = 1000;
const bool DEFAULT_OUT = false;
const char DEFAULT_FILE_DIR[50] = "lapFusion_results.tsv";
const char ARGS[ARGS_NUM][10] = {"-N", "-I", "-T", "-O", "-F"};

struct Args{
  int n;
  int iter_max;
  int num_threads;
  bool out;
  char *file_dir;
};

bool startsWith(const char *string, const char *prefix){
  return strncmp(string, prefix, strlen(prefix)) == 0;
}

bool endsWith(const char *string, const char *sufix){
  string = strrchr(string, '.');
  if (string != NULL)
    return strncmp(string, sufix, strlen(sufix)) == 0;
  return false;
}

struct Args args_parser(int argc, char *argv[])
{
  int n = DEFAULT_N;
  int iter_max = DEFAULT_ITER_MAX;
  int num_threads = omp_get_num_procs();
  bool out = DEFAULT_OUT;
  char *file_dir = malloc(50);
  strcpy(file_dir, DEFAULT_FILE_DIR);

  for (int i = 1; i < argc; i++){
    for (int j = 0; j < ARGS_NUM; j++){
      if (startsWith(*(argv + i), ARGS[j])){
        if (j == 3) out = true;
        if (i + 2 <= argc)
        {
          if (j == 0) n = atoi(*(argv + i + 1));
          if (j == 1) iter_max = atoi(*(argv + i + 1));
          if (j == 2) num_threads = atoi(*(argv + i + 1));
          if (j == 4) strcpy(file_dir, *(argv + i + 1));
        }
        else{
          printf("Wrong command, using default values.\n");
        }
        break;
      }
    }
  }
  if (n<1 | n>8192){
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
  if (!endsWith(file_dir, ".tsv")){
    printf("Wrong \'file_dir\', using default name (%s)\n", DEFAULT_FILE_DIR);
    strcpy(file_dir, DEFAULT_FILE_DIR);
  }
  struct Args parsed_args = {n, iter_max, num_threads, out, file_dir};
  return parsed_args;
}
