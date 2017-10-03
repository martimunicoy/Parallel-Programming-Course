/*First Steps with C language: Merge Sort (Divide & Conquer algorithmic strategy)*/

// For printf() function:
#include <stdio.h>
// For malloc() function:
#include <stdlib.h>
// For time() function:
#include <time.h>

/*Functions*/

void print_array(int ptr_len, float *ptr){
	for (int i = 0; i < ptr_len; i++){
		printf("%f\t", *(ptr + i));
	}
	printf("\n");
}


void zero_array(int ptr_len, float *ptr){
	for (int i = 0; i < ptr_len; i++){
		*(ptr + i) = 0.0;
	}
}


void rand_array(int ptr_len, float *ptr){
	for (int i = 0; i < ptr_len; i++){
		float r = rand() % 1000000;
		*(ptr + i) = r / 1000000;
	}
}


void copy_array(int len, float *original, float *clone){
	for (int i = 0; i < len; i++){
		*(clone + i) = *(original + i);
	}
}

void Merge(float *Y, float *X1, float *X2, int N1, int N2){
	// Both lists must be individually sorted. Check it:
	/*printf("Merging:\n");
	print_array(N1, X1);
	print_array(N2, X2);*/

	int i = 0;
	int j = 0;
	while(i < N2){
		while((*(X2 + i) > *(X1 + j)) & (j < N1)){
			*(Y + i + j) = *(X1 + j);
			j = j + 1;
		}
		*(Y + i + j) = *(X2 + i);
		i = i + 1;
	}
	while(j < N1){
		*(Y + i + j) = *(X1 + j);
		j = j + 1;
	}

	// See the results
	/*printf("\nResult\n");
	print_array(N1 + N2, Y);*/
}


void MergeSort(float *X, float *Y, int N){
	if (N==1){
		copy_array(N, X, Y);
		return;
	}

	float *T;
	T = (float*) malloc(N * sizeof(float));

	MergeSort(X, T, N/2);
	MergeSort(X + N/2, T + N/2, N - N/2);
	Merge(Y, T, T + N/2, N/2, N - N/2);	
}


int main(){
	// Initiate seed for random numbers
	srand((unsigned)time(NULL));

	// Initial variables
	int N = 20;

	// Dynamically store memory space
	float *X_ptr, *Y_ptr;
	X_ptr = (float*) malloc(N * sizeof(float));
	Y_ptr = (float*) malloc(N * sizeof(float));

	// Create an arrary of random numbers
	rand_array(N, X_ptr);

	// Call MergeSort function
	MergeSort(X_ptr, Y_ptr, N);

	// See the Results
	printf("X:\n");
	print_array(N, X_ptr);
	printf("Y:\n");
	print_array(N, Y_ptr);
}