#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

// MPI methods
int MPI_Exchange(int n, double * array, int rank1, int rank2, int direction, MPI_Comm comm, double *compT, double *commT);
int MPI_Sort_bitonic(int n, double * array, MPI_Comm comm);
  
// merge functions
double * merge(int n, double * array, int m, double * b);
void merge_sort(int n, double * array);
double * merge_inv(int n, double * array, int m, double * b);
int reorder(int n, int direction, double * array);
void reverse(int n, double * array);
void swap (double * array, double * b);

int main(int argc, char ** argv) {
  //setup
  int size, rank, result, i;
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
  result = MPI_Sort_bitonic(n, array, MPI_COMM_WORLD);
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

int MPI_Sort_bitonic(int n, double * array, MPI_Comm comm) {
  
  int rank, size, result, pair, direction;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  double compT, commT, time, time2;

  // allocate space for numElements/numProcessors amount of doubles
  double * localArray = (double *)calloc(n/size, sizeof(double));

  // scater array to localArray
  time = MPI_Wtime();
  result = MPI_Scatter(array, n/size, MPI_DOUBLE, localArray, n/size, MPI_DOUBLE, 0, comm);
  if(result != MPI_SUCCESS) return result;
  time = MPI_Wtime() - time;
  commT += time;

  // order localArray
  time = MPI_Wtime();
  merge_sort(n/size, localArray);
  time = MPI_Wtime() - time;
  compT += time;

  for(int i = 1; i <= log2(size); i++) {
    time = MPI_Wtime();
    if(rank&(1<<i)) {
      direction = -1;
    } else {
      direction = 1;
    }
    if(reorder(n/size, direction, localArray)) reverse(n/size, localArray);
    time = MPI_Wtime() - time;
    compT += time;
    //printf("level %d, rank %d: direction %d rank&(1<<i) %d\n", i, rank, direction, rank&(1<<i));
    for(int j = i-1; j >= 0; j--) {
      pair = rank^(1<<(j));
      //if(rank<pair) printf("round %d: rank %d, pair %d\n", j, rank, pair);
      time = 0;
      time2 = 0;
      if(pair<rank) MPI_Exchange(n/size, localArray, pair, rank, direction, comm, &time, &time2);
      else if(pair>rank) MPI_Exchange(n/size, localArray, rank, pair, direction, comm, &time, &time2);
      compT += time;
      commT += time2;
    }
  }

  // gather localArray
  time = MPI_Wtime();
  result = MPI_Gather(localArray, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, 0, comm);
  if(result != MPI_SUCCESS) return result;
  time = MPI_Wtime() - time;
  commT += time;

  printf("Processor %d: commT = %lf, compT = %lf \n", rank, commT, compT);

  return MPI_SUCCESS;
}

int MPI_Exchange(int n, double * array, int rank1, int rank2, int direction, MPI_Comm comm, double *compT, double *commT) {
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
    if(result != MPI_SUCCESS) return result;
    result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank2, tag2, comm, &status);
    if(result != MPI_SUCCESS) return result;
    time = MPI_Wtime() - time;
    *commT += time;
    time = MPI_Wtime();
    if(direction == 1) {
      c = merge(n, array, n, b);
    } else if(direction == -1) {
      c = merge_inv(n, array, n, b);
    }
    for(i = 0; i < n; i++) {
      array[i] = c[i];
    }
    time = MPI_Wtime() - time;
    *compT += time;
  }
  else if(rank == rank2) {
    time = MPI_Wtime();
    result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank1, tag1, comm, &status);
    if(result != MPI_SUCCESS) return result;
    result = MPI_Send(&array[0], n, MPI_DOUBLE, rank1, tag2, comm);
    if(result != MPI_SUCCESS) return result;
    time = MPI_Wtime() - time;
    *commT += time;
    time = MPI_Wtime();
    if(direction == 1) {
      c = merge(n, array, n, b);
    } else if(direction == -1) {
      c = merge_inv(n, array, n, b);
    }
    for(i =0; i < n; i++) {
      array[i] = c[i + n];
    }
    time = MPI_Wtime() - time;
    *compT += time;
  }
  return MPI_SUCCESS;
}

double * merge(int n, double * a, int m, double * b) {
    int i, j, k;
    double * c = (double *) calloc(n + m, sizeof(double));
  
    for(i=j=k=0; (i < n) && (j < m);) {
      if(a[i] <= b[j]) {
        c[k++] = a[i++];
      }
      else {
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
  
double * merge_inv(int n, double * a, int m, double * b) {
  int i, j, k;
  double * c = (double *)calloc(n + m, sizeof(double));

  for (i = 0, j = 0, k = 0; (i < n) && (j < m);) {
      if (a[i] >= b[j]) {
          c[k++] = a[i++];
      }
      else {
          c[k++] = b[j++];
      }
  }
  if (i == n) {
      for (; j < m;) {
          c[k++] = b[j++];
      }
  }
  else {
      for (; i < n;) {
          c[k++] = a[i++];
      }
  }
  return c;
}

int reorder(int n, int direction, double * array) {
  if(direction == 1) {
    for(int i = 1; i < n; i++) {
      if(array[i] > array[i-1]) return 0;
      else if(array[i] < array[i-1]) return 1;
    }
  } else if(direction == -1) {
    for(int i = 1; i < n; i++) {
      if(array[i] > array[i-1]) return 1;
      else if(array[i] < array[i-1]) return 0;
    }
  }
  return 0;
}

void reverse(int n, double *array) {
  double temp;
  int i, j;

  for (i = 0, j = n - 1; i < j; i++, j--) {
    temp = array[i];
    array[i] = array[j];
    array[j] = temp;
  }
}

void swap (double * a, double * b) {
  double temp;
  temp = *a;
  *a = *b;
  *b = temp;
}
