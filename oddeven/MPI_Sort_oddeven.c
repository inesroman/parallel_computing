#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
      
//MPI methods
int MPI_Exchange(int n, double * array, int rank1, int rank2, MPI_Comm comm, double *compT, double *commT);
int MPI_Sort_oddeven(int n, double * array, int root, MPI_Comm comm);
int isArraySorted(int n, double * array, int root, MPI_Comm comm, double *compT, double *commT);
  
//all in previous labs
double * merge(int n, double * array, int m, double * b);
void merge_sort(int n, double * array);
void swap (double * array, double * b);
  
// function definitions
int main(int argc, char ** argv) {
  //setup
  int size, rank, result, i, *answer;
  int n = 10000000;
  double m = 100.0;
  double *array;
  
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  //allocate space for an array of doubles, size n
  array = (double *) calloc(n, sizeof(double));
    
  //fills array with random values on root proc
  if(rank == 0) {
    //get random values for array & output for testing
    srand(((unsigned) time(NULL) + rank));
    for(i = 0; i < n; i++) {
      array[i] = ((double) rand() / RAND_MAX) * m;
      //printf("Initial: %f\n", array[i]);
    }
  }

  //MPI_Sort does all the heavy work
  double time = MPI_Wtime(), overallTime;
  result = MPI_Sort_oddeven(n, array, 0, MPI_COMM_WORLD);
  if(result != MPI_SUCCESS) {
    return result;
  }
  time = MPI_Wtime()-time;
  MPI_Reduce(&time, &overallTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  
  //output ordered list for testing
  if(rank == 0) {
    //for(i = 0; i < n; i++) printf("Output : %f\n", array[i]);
    printf("Overall Exec time = %lf with %d procs\n", overallTime, size);
  }
  MPI_Finalize();
}
  
int MPI_Sort_oddeven(int n, double * array, int root, MPI_Comm comm) {
  
  // get rank and size of comm
  int size, rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  double compT = 0, commT = 0, time, time2;

  //allocate space for numElements/numProcessors amount of doubles
  double * localA = (double *)calloc(n/size, sizeof(double));
    
  //scatter a to local_a
  time = MPI_Wtime();
  int rc = MPI_Scatter(array, n/size, MPI_DOUBLE, localA, n/size, MPI_DOUBLE, root, comm);
  if(rc != MPI_SUCCESS) return rc;
  time = MPI_Wtime() - time;
  commT += time;
  
  //sort local_a using mergeSort
  time = MPI_Wtime();
  merge_sort(n/size, localA);
  time = MPI_Wtime() - time;
  compT += time;
  
  //odd-even iterations
  for(int step=0;step<size;step++) {
    if((step+rank)%2==0) {
      if(rank<size-1) {
        time = 0;
        time2 = 0;
        MPI_Exchange(n/size, localA, rank, rank+1, comm, &time, &time2);
        compT += time;
        commT += time2;
      }
    } else {
      if(rank>0) {
        time = 0;
        time2 = 0;
        MPI_Exchange(n/size, localA, rank-1, rank, comm, &time, &time2);
        compT += time;
        commT += time2;
      }
    }
    // test is localA sorted
    time = 0;
    time2 = 0;
    if(isArraySorted(n/size, localA, root, comm, &time, &time2)) {
      //printf("STOPPED AFTER %d REPETIONS\n", step);
      break;
    }
    compT += time;
    commT += time2;
  }

  //gather local_a
  time = MPI_Wtime();
  rc = MPI_Gather(localA, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, root, comm);
  if(rc != MPI_SUCCESS) return rc;
  time = MPI_Wtime() - time;
  commT += time;

  printf("Processor %d: commT = %lf, compT = %lf \n", rank, commT, compT);
  
  return MPI_SUCCESS;
}

int isArraySorted(int n, double * array, int root, MPI_Comm comm, double *compT, double *commT) {
  // get rank and size
  int size, rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  double time;
  // gather the first and last elements of array
  double * first = (double *)calloc(size, sizeof(double));
  double * last  = (double *)calloc(size, sizeof(double));
  time = MPI_Wtime();
  MPI_Gather(&array[0],   1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm);
  MPI_Gather(&array[n-1], 1, MPI_DOUBLE, last,  1, MPI_DOUBLE, root, comm);
  time = MPI_Wtime() - time;
  *commT += time;
  
  // if root then test array is sorted
  int answer = 1;
  if(rank==root) {
    time = MPI_Wtime();
    for(int i=0;i<size-1;i++) {
      if(last[i]>first[i+1]) {
        answer = 0;
      }
    }
    time = MPI_Wtime() - time;
    *compT += time;
  }
  // Bcast answer
  time = MPI_Wtime();
  MPI_Bcast(&answer, 1, MPI_INT, root, comm);
  time = MPI_Wtime() - time;
  *commT += time;
  return answer;
}
  
int MPI_Exchange(int n, double * array, int rank1, int rank2, MPI_Comm comm, double *compT, double *commT) {
  int rank, size, result, i, tag1 = 0, tag2 = 1;
  double * b = (double *) calloc(n, sizeof(double));
  double * c, time;
    
  MPI_Status status;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  //L8.6
  if(rank == rank1) {
    time = MPI_Wtime();
    result = MPI_Send(&array[0], n, MPI_DOUBLE, rank2, tag1, comm);
    result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank2, tag2, comm, &status);
    time = MPI_Wtime() - time;
    *commT += time;
    time = MPI_Wtime();
    c = merge(n, array, n, b);
    for(i = 0; i < n; i++) {
      array[i] = c[i];
    }
    time = MPI_Wtime() - time;
    *compT += time;
  }
  else if(rank == rank2) {
    time = MPI_Wtime();
    result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank1, tag1, comm, &status);
    result = MPI_Send(&array[0], n, MPI_DOUBLE, rank1, tag2, comm);
    time = MPI_Wtime() - time;
    *commT += time;
    time = MPI_Wtime();
    c = merge(n, array, n, b);
    for(i =0; i < n; i++) {
      array[i] = c[i + n];
    }
    time = MPI_Wtime() - time;
    *compT += time;
  }
  return MPI_SUCCESS;
}
  
//notes
double * merge(int n, double * a, int m, double * b) {
  int i, j, k;
  double * c = (double *) calloc(n + m, sizeof(double));

  for(i=j=k=0; (i < n) && (j < m);) {
    if(a[i] <= b[j]) {
      c[k++] = a[i++];
    } else {
      c[k++] = b[j++];
    }
  }
  if(i == n) {
    for(; j < m;) {
      c[k++] = b[j++];
    }
  } else {
    for(; i < n;) {
      c[k++] = a[i++];
    }
  }
  return c;
}
  
//notes
void merge_sort(int n, double * a) {
  double * c;
  int i;
  
  if (n <= 1) {
    return;
  }
  if(n == 2) {
    if(a[0] > a[1]) {
      swap(&a[0], &a[1]);
    }
    return;
  }
  
  merge_sort(n / 2, a);
  merge_sort(n - n / 2, a + n / 2);
  c = merge(n / 2, a, n - n / 2, a + n / 2);
  for(i = 0; i < n; i++) {
    a[i] = c[i];
  }
}
  
//notes
void swap (double * a, double * b) {
  double temp;
  temp = *a;
  *a = *b;
  *b = temp;
}
