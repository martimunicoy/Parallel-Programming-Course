/*Program Examples: Finite Differences*/

// For printf() function:
#include <stdio.h>
// For malloc() function:
#include <stdlib.h>
// To work with bool functions:
#include <stdbool.h>
// To work with strlen() function:
#include <string.h>
// We need to use some math functions:
#include <math.h>
// For nanosleep() function:
#include <time.h>

const int DEFAULT_X = 30;
const int DEFAULT_T = 10;
const double L = 0.345678;
const int ARGS_NUM = 2;
const float SLEEP = 100000000L;
const char ARGS[ARGS_NUM][10] = {"-X", "-T"};


bool startsWith(const char *string, const char *prefix){
	return strncmp(string, prefix, strlen(prefix)) == 0;
}


struct Pair{
	int X;
	int T;
};


struct Pair arg_parser(int argc, char *argv[]){
	int X = DEFAULT_X;
	int T = DEFAULT_T;
	for (int i = 1; i < argc; i++){
		for (int j = 0; j < ARGS_NUM; j++){
			if (startsWith(*(argv + i), ARGS[j])){
				if (i + 2 <= argc){
					if (j == 0) X = atoi(*(argv + i + 1));
					if (j == 1) T = atoi(*(argv + i + 1));
				}
				else{
					printf("Wrong command, using default values.\n");
				}
				break;
			}
		}
	}
	if (X<1 | X>999999){
		printf("X out of range, using default value (%d)\n", DEFAULT_X);
		X = DEFAULT_X;
	}
	if (T<1 | T>999999){
		printf("T out of range, using default value (%d)\n", DEFAULT_T);
		T = DEFAULT_T;
	}
	struct Pair CO = {X, T};
	return CO;
}


struct Matrix{
	// Number of columns
	int x;
	// Number of rows
	int y;
	double **elements;
};


struct Matrix new_matrix(int x, int y){
	double **U;
	U = malloc(x * sizeof(double *));
	for (int i = 0; i < x; i++){
		U[i] = malloc(y * sizeof(double));
	}
	struct Matrix M = {x, y, U};
	return M;
}


void initialize(struct Matrix U, struct Pair CO){
	for (int j = 0; j < U.x; j++){
		U.elements[j][0] = 0.0f;
		U.elements[j][U.y - 1] = 0.0f;
	}
	for (int i = 1; i < U.y - 1; i++){
		U.elements[0][i] = sin(i * M_PI / CO.X);
	}
	for (int i = 1; i < U.y - 1; i++){
		 U.elements[1][i] = sin(i * M_PI / CO.X) * cos(M_PI / CO.T);
	}
}


void propagate(struct Matrix U, struct Pair CO){
	for (int j = 2; j < U.x; j++){
		for (int i = 1; i < U.y - 1; i++){
			//printf("%f %f %f %f = %f\n", U.elements[j-1][i], U.elements[j-1][i+1], U.elements[j-1][i-1], U.elements[j-2][i], 2 * (1 - L) * U.elements[j-1][i] + L * U.elements[j-1][i+1] + L * U.elements[j-1][i-1] - U.elements[j-2][i]);
			U.elements[j][i] = 2 * (1 - L) * U.elements[j-1][i] + L * U.elements[j-1][i+1] + L * U.elements[j-1][i-1] - U.elements[j-2][i];
			//U.elements[j][i] = 2 * U.elements[j-1][i] - U.elements[j-2][i] + L * L * (U.elements[j-1][i+1] + U.elements[j-1][i-1] - 2 * U.elements[j-1][i]);
		}
	}
}

void display(struct Matrix U, struct Pair CO){
	FILE *gnuplot = popen("gnuplot -persistent", "w");
	fprintf(gnuplot, "set xrange [0:%i]\n", CO.X);
	fprintf(gnuplot, "set yrange [-1.5:1.5]\n");
	for (int j = 0; j < U.x; j++){
		fprintf(gnuplot, "clear\n");
		fprintf(gnuplot, "plot '-'\n");
		int i;
		for (i = 0; i < U.y; i++){
		    fprintf(gnuplot, "%d %f\n", i, U.elements[j][i]);
		}
		fprintf(gnuplot, "e\n");
		fflush(gnuplot);
		nanosleep((const struct timespec[]){{0, SLEEP}}, NULL);
		}
}


void print_matrix(struct Matrix U){
	for (int i = 0; i < U.y; i++){
		for (int j = 0; j < U.x; j++) printf("\t%00.12f", U.elements[j][i]);
		printf("\n");
	}
}


int main(int argc, char *argv[]){
	// Initial variables from arg_parser
	struct Pair CO = arg_parser(argc, argv);

	// Dynamically store a matrix U in memory
	int x, y;
	x = CO.T + 2;
	y = CO.X + 1;
	struct Matrix U = new_matrix(x, y);

	// Initialize state
	initialize(U, CO);

	// Propagate simulation
	propagate(U, CO);

	// Display solution
	display(U, CO);

	// print_matrix(U);
	return 0;
}