#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>


double * merge_array(int n, double * a, int m, double * b);
void     merge_sort(int n, double * a);
void     swap (double * a, double * b);

int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm);

int main (int argc, char *argv[]) {

	int rank, size;

	int n = 10000000, i;
	double m = 100.0;
	double * array;

	// Init + rank + size
	MPI_Init(&argc, &argv);
   	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   	MPI_Comm_size(MPI_COMM_WORLD, &size);
    array = (double *) calloc(n, sizeof(double));

	if(rank == 0) {
	   //initialise the array with random values, then scatter to all processors
	   srand(((unsigned)time(NULL)+rank));
	   for(i = 0; i < n; i++) {
	      array[i]=((double)rand()/RAND_MAX)*m;
	   }
	}

    // call and time evaluate MPI_Sort_direct
    double time = MPI_Wtime(), overallTime;
    MPI_Sort_bucket(n, array, m, 0, MPI_COMM_WORLD);
    time = MPI_Wtime()-time;
    MPI_Reduce(&time, &overallTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if(rank ==0) {
        //for(int i=0;i<n;i++)printf("%lf \n", array[i]);
        printf("Overall Exec time = %lf with %d procs\n", overallTime, size);
    }
	MPI_Finalize();
}

// MPI Functions

int MPI_Sort_bucket(int n, double * array, double max, int root, MPI_Comm comm){
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    double compT = 0, commT = 0, time;
    
    double * localBucket = (double *)calloc(n, sizeof(double));
    
    // bcast array
    time = MPI_Wtime();
    int rc = MPI_Bcast(array, n, MPI_DOUBLE, root, comm);
    if(rc != MPI_SUCCESS)return rc;
    time = MPI_Wtime() - time;
    commT += time;

    // filter array to get the elems of localBucket
    int count = 0;
    time = MPI_Wtime();
    for(int i=0;i<n;i++) {
        //check if array[i] is in the range of bucket
        if(array[i]>= rank*max/size && array[i]<(rank+1)*max/size){
            localBucket[count++] = array[i];
        }
    }
    time = MPI_Wtime() - time;
    compT += time;

    
    // sort local bucket
    time = MPI_Wtime();
    merge_sort(count, localBucket);
    time = MPI_Wtime() - time;
    compT += time;

    //gather counts
    int * counts = (int *)calloc(size, sizeof(int));
    int * displs = (int *)calloc(size, sizeof(int));
    time = MPI_Wtime();
    MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, root, comm);
    time = MPI_Wtime() - time;
    commT += time;

    if(rank == root){
        //calculate displs
        time = MPI_Wtime();
        displs[0] = 0;
        for(int i=1;i<size;i++){
            displs[i] = displs[i-1]+counts[i-1];
        }
        time = MPI_Wtime() - time;
        compT += time;
    }

    time = MPI_Wtime();
    rc = MPI_Gatherv(localBucket, count, MPI_DOUBLE, array, counts, displs, MPI_DOUBLE, root, comm);
    if(rc != MPI_SUCCESS)return rc;
    time = MPI_Wtime() - time;
    commT += time;

    printf("Processor %d: commT = %lf, compT = %lf \n", rank, commT, compT);

    return MPI_SUCCESS;

}

// function to merge the array a with n elements with the array b with m elements
// function returns the nerged array

double * merge_array(int n, double * a, int m, double * b){
    int i,j,k;
    double * c = (double *) calloc(n+m, sizeof(double));

    for(i=j=k=0;(i<n)&&(j<m);)
        if(a[i]<=b[j])c[k++]=a[i++];
        else c[k++]=b[j++];

    if(i==n)for(;j<m;)c[k++]=b[j++];
    else for(;i<n;)c[k++]=a[i++];

    return c;
}

// function to merge sort the array a with n elements

void merge_sort(int n, double * a){
    double * c;
    int i;

    if (n<=1) return;

    if(n==2) {

        if(a[0]>a[1])swap(&a[0],&a[1]);
        return;
    }
    merge_sort(n/2,a);merge_sort(n-n/2,a+n/2);

    c=merge_array(n/2,a,n-n/2,a+n/2);

    for(i=0;i<n;i++)a[i]=c[i];

    return;
}


// swap two doubles
void swap (double * a, double * b){
    double temp;
    temp=*a;*a=*b;*b=temp;
}
