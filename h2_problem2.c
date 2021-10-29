#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#define ROOT 0

int comm_sz = 0;
int rank = 0;
MPI_Comm comm = MPI_COMM_WORLD;

int GenRandNumber(int max, int min);
void Quicksort(int *arr, int left, int right);
int Partition(int *arr, int left, int right);
void Swap(int *a, int *b);
int ComputePartner(int phase, int my_rank);


int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(comm, &comm_sz);
	MPI_Comm_rank(comm, &rank);

	double startTime = 0.0, totalTime = 0.0;

	//set random seed
	srand( time(NULL) + rank);
	fflush(stdout);
	//rank 0 read n and Broadcas to all processor
	int n = 0;
	if( rank == ROOT){
		printf("Input n:");
		scanf("%d", &n);
	}
	MPI_Bcast(&n, 1, MPI_INT, ROOT, comm);

	//global list has n keys
	int *global_list = malloc(n * sizeof(int));

	//the local list has n/comm_sz keys, and generater
	int i = 0;
	int num_keys = n / comm_sz;
	int *local_list = malloc(num_keys * sizeof(int));
	for( i = 0; i < num_keys; i++)
		local_list[i] = GenRandNumber(100, 1);
	
	//sort the local list
	Quicksort(local_list, 0, num_keys - 1);	
	
	//gather the local list and proc 0 print the global list
	MPI_Gather(local_list, num_keys, MPI_INT, global_list, num_keys, MPI_INT, ROOT,comm);
	if( rank == ROOT){
		printf("After local sort:");
		fflush(stdout);
		for( i = 0; i < n; i++){
			printf("%d%c", global_list[i], (i < n-1)? ' ' : '\n');
			fflush(stdout);	
		}
	}

	//odd even sort
	int phase = 0, partner = 0;
	for( phase = 0; phase < comm_sz; phase++){
		partner = ComputePartner(phase, rank);
		printf("phase %d) %d <-> %d\n", phase, rank, partner);
		fflush(stdout);
	}
	
	//gather the sorted local list and print the global list
	//MPI_Gather(local_list, num_keys, MPI_INT, global_list, num_keys, MPI_INT, ROOT, comm);
	if( rank == ROOT){
		printf("After odd-even sort:");
		fflush(stdout);
		for( i = 0; i < n; i++){
			printf("%d%c", global_list[i], (i < n-1)? ' ' : '\n');
			fflush(stdout);	
		}
	}

	free(global_list);
	free(local_list);
	MPI_Finalize();
	return 0;
}

int GenRandNumber(int max, int min)
{
	return (rand()%(max-min+1)) + min;
};

void Quicksort(int *arr, int left, int right)
{
	if(left < right){
		int q = Partition(arr, left, right);
		Quicksort(arr, left, q-1);
		Quicksort(arr, q+1, right);
	}
};

int Partition(int *arr, int left, int right)
{
	int i = left - 1;
	int j = 0;
	for(j = left; j < right; j++){
		if(arr[j] <= arr[right]){
			i++;
			Swap(&(arr[i]), &(arr[j]));
		}		
	}
	Swap(&(arr[i+1]), &(arr[right]));
	return i+1;
};

void Swap(int *a, int *b)
{
	int t = *a;
	*a = *b;
	*b = t;
	
};
int ComputePartner(int phase, int my_rank)
{
	int partner = 0;	
	if(phase % 2 == 0){				/* Even phase */
		if(my_rank % 2 != 0)			/* Odd rank */
			partner = my_rank - 1;
		else							/* Even rank */
			partner = my_rank + 1;
	}
	else{							/* Odd phase */
		if(my_rank % 2 != 0)			/* Odd rank */
			partner = my_rank + 1;
		else							/* Even rank */
			partner = my_rank - 1;
	}	
	if(partner == -1 || partner == comm_sz)
		partner = MPI_PROC_NULL;
	return partner;
};
