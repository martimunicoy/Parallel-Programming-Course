/*Program Examples: Heat Transfer*/

// Include necessary libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>


const float DEFAULT_L = 0.345678;
const int DEFAULT_T = 10;
const int DEFAULT_N = 10;
const bool DEFAULT_V = 0;
const int ARGS_NUM = 4;
const char ARGS[ARGS_NUM][10] = {"-L", "-T", "-N", "-V"};


bool startsWith(const char *string, const char *prefix){
	return strncmp(string, prefix, strlen(prefix)) == 0;
}


struct Args{
	float L;
	int T;
	int N;
	bool V;
};


struct Args arg_parser(int argc, char *argv[]){
	float L = DEFAULT_L;
	int T = DEFAULT_T;
	int N = DEFAULT_N;
	bool V = DEFAULT_V;
	int e = 0;
	int com[3] = {-1, -1, -1};

	for (int i = 1; i < argc; i++){
		for (int j = 0; j < ARGS_NUM; j++){
			if (startsWith(*(argv + i), ARGS[j])){
				if (j == 3) V = 1;
				if (i + 2 <= argc){
					if (j == 0) L = atof(*(argv + i + 1));
					if (j == 1) T = atoi(*(argv + i + 1));
					if (j == 2) N = atoi(*(argv + i + 1));
				}
				else{
					com[e] = j;
					e++;
				}
				break;
			}
		}
	}

	for (int k = 0; k < e; k++){
		if (com[k] == 0) printf("Wrong command, using default value for L (%f)\n", DEFAULT_L);
		if (com[k] == 1) printf("Wrong command, using default value for T (%d)\n", DEFAULT_T);
		if (com[k] == 2) printf("Wrong command, using default value for N (%d)\n", DEFAULT_N);
	}

	if (L<=0 | L>0.5){
		printf("L out of range, using default value (%f)\n", DEFAULT_L);
		L = DEFAULT_L;
	}
	if (T<=0 | T>99999){
		printf("T out of range, using default value (%d)\n", DEFAULT_T);
		T = DEFAULT_T;
	}
	if (N<1 | N>99999){
		printf("N out of range, using default value (%d)\n", DEFAULT_N);
		N = DEFAULT_N;
	}
	struct Args CO = {L, T, N, V};
	return CO;
}


struct Array{
	int size;
	float * elements;
};


struct Array new_array(int size){
	float * U;
	U = malloc(size * sizeof(float));
	memset(U, 0.0, size * sizeof(float));
	struct Array A = {size, U};
	return A;
}


void print_array(struct Array U){
	for (int i = 0; i < U.size; i++){
		printf("%f ", U.elements[i]);
	}
	printf("\n");
}


void initiate(struct Array U){
	int i;
	for (i = 0; i < U.size - 1; i++){
		U.elements[i] = 20;
	}
	U.elements[i] = 40;
}


void simulate(struct Array * U, struct Array * newU, int counter, struct Args CO){
	if (counter > 0){
		for (int i = 1; i < (*U).size - 1; i++){
			(*newU).elements[i] = (1.0 - 2 * CO.L) * (*U).elements[i] + CO.L * ((*U).elements[i-1] + (*U).elements[i+1]);
		}
		counter--;

		if (CO.V){
			printf("Iteration %d: ", counter);
			print_array(*newU);
		}

		return simulate(newU, U, counter, CO);
	}
}


void simulate2(struct Array * U, struct Array * newU, int counter, struct Args CO){
	struct Array * ptr;
	struct Array * ptr1 = U;
	struct Array * ptr2 = newU;

	while (counter > 0){
		for (int i = 1; i < (*U).size - 1; i++){
			(*ptr2).elements[i] = (1.0 - 2 * CO.L) * (*ptr1).elements[i] + CO.L * ((*ptr1).elements[i-1] + (*ptr1).elements[i+1]);
		}
		counter--;

		if (CO.V){
			printf("Iteration %d: ", counter);
			print_array(*ptr2);
		}

		ptr = ptr1;
		ptr1 = ptr2;
		ptr2 = ptr;


	}
}


void copy_array(struct Array newU, struct Array U){
	for (int i = 0; i < U.size; i++){
		U.elements[i] = newU.elements[i];
	}
}


int main(int argc, char *argv[]){
	// Initial variables from arg_parser
	struct Args CO = arg_parser(argc, argv);

	// Dynamically store arrays U and newU in memory
	struct Array U = new_array(CO.N + 2);
	struct Array newU = new_array(CO.N + 2);

	// Put U and newU in initial state
	initiate(U);
	initiate(newU);

	// Simulation
	simulate2(&U, &newU, CO.T, CO);

	// Get back the best result
	if (CO.T%2 == 1){
		copy_array(newU, U);
	}

	// Print final U
	print_array(U);

	// Free allocated memory
	free(U.elements);
}