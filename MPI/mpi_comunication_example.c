#include <mpi.h>
#include <stdio.h>
#include <string.h>

int main( int argc, char *argv[] ){
	/* message (other types possible) */
	int MAXSIZE = 60;
	char mess[MAXSIZE];
	/*Sheared variables by all procesors*/
	int myRank, numProc, count;

	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
	MPI_Comm_size( MPI_COMM_WORLD, &numProc );
	
	if (myRank != 0){
		/* all processes send to root
		*/
		/* create message */
		sprintf(mess, "Hello from %d", myRank);
		/* destination is root */
		int dest = 0;
		int tag = 32;
		count = strlen(mess) + 1; /* because of '\0' */
		MPI_Send(mess, count, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
	}
	else{
		/* root (0) process receives and prints messages */
		/* from each processor in rank order
		*/
		int source;
		int tag = 32;
		for(source = 1; source < numProc; source++){
			MPI_Status status;
			MPI_Recv(mess,MAXSIZE,MPI_CHAR,source,tag,MPI_COMM_WORLD, &status);
			MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
			printf("The proccesor %d is reciving the following message: %s\n", myRank, mess);
		}
	}


	MPI_Finalize();
	return 0;
}
