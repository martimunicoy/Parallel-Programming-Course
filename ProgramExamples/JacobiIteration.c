/*First Steps with C language: Solving Laplace Equation*/

#include <stdio.h>
#include <math.h>

/*Global variables*/

int m = 10;
int n = 20;
float tol = 0.001;

/*Functions*/

float mean(float u, float d, float r, float l){
	float result;
	result = (u + d + r + l)/4;

	return result;
}


void print_matrix(int m, int n, float A[][m]){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++) printf("\t%f", A[i][j]);
		printf("\n");
	}

	printf("\n");

}


void initialize_matrix(int m, int n, float A[][m]){
	for (int i = 0; i<n; i++){
		A[i][0] = sin(i * M_PI/(n - 1));
		A[i][m-1] = sin(i * M_PI/(n - 1)) * exp(-M_PI);
	}

	for (int i = 0; i < n; i++){
		for (int j = 1; j < m-1; j++) A[i][j] = 0.0;
	}
}


void copy_matrix(int m, int n, float A[][m], float B[][m]){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m ; j++) B[i][j] = A[i][j];
	}
}


float get_max_error(float max_error, float a, float anew){
	float error = sqrt(fabs(a - anew));

	if(error > max_error){
		return error;
	}else return max_error;
}


void jacobi_iteration(int m, int n, float A[][m]){
	int it = 0;
	float Anew[n][m];
	float max_error = tol + 1;

	printf(" Tolerance Goal:             \t%f\n", tol);
	printf(" ---------------------------------------\n");
	copy_matrix(m, n, A, Anew);

	while (max_error > tol){
		max_error = 0;
		for (int i = 1; i < n-1; i++){
			for (int j = 1; j < m-1; j++){
				float Anew_element = mean(A[i][j+1], A[i][j-1], A[i+1][j], A[i-1][j]);
				Anew[i][j] = Anew_element;
				max_error = get_max_error(max_error, A[i][j], Anew_element);
				//printf("max_error %f\n", max_error);
			}
		}
		copy_matrix(m, n, Anew, A);

		if (it%10 == 0) printf(" Error at iteration num. %i:\t%f\n", it, max_error);

		it = it + 1;

		//print_matrix(m, n, Anew);
	}

	printf(" ---------------------------------------\n");
	printf(" Maximum Error achieved:      \t%f\n", max_error);
	printf("\nRESULTS (after %i iterations):\n\n", it-1);
	print_matrix(m, n, A);
}

/*Main function*/ 

int main(){
	float A[n][m];
	initialize_matrix(m, n, A);
	jacobi_iteration(m, n, A);
	return 0;
}


/*Dumping functions*/

/*
float get_max_diff(float max_diff, float a, float anew){
	float diff = fabs(a - anew);
	if(diff > max_diff){
		return diff;
	}else return max_diff;
}


void zero_matrix(int m, int n, float A[][m]){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m ; j++){
			A[i][j] = 0;
		}
	}
}
*/