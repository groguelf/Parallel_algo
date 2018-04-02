#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int **initialize_matrix(int number){
  int **mat = calloc(number, sizeof(int *));
  for (int i = 0; i < number; i++){
    mat[i] = calloc(number, sizeof(int));
  }
  return mat;
}

void print_matrix(int **matrix, int size){
  for (int i = 0; i < size; i++){
    for(int j = 0; j < size; j++){
      printf("%d     ", matrix[i][j]);
    }
    printf("\n");
  }
}

int main(int argc, char **argv){
  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  float float_root = sqrt(size);
  int root = (int) float_root;
  if (float_root * float_root != (float) size){
    if (rank == 0){
      printf("You have to use a perfect square number of processes.\n\n");
    }
    MPI_Finalize();
    return 0;
  }

  int nnodes = size;
  int ndims = 2;
  int dims[ndims];
  int periods[ndims];
  int reorder = 0;
  MPI_Comm comm, new_comm;
  int remain_dims[ndims];

  for (int i = 0; i < ndims; i++){
    dims[i] = root;
    periods[i] = 0;
    remain_dims[i] = 0;
  }

  remain_dims[1] = 1;

  MPI_Dims_create(nnodes, ndims, dims);
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,periods, reorder, &comm);
  MPI_Cart_sub(comm, remain_dims, &new_comm);

  int coords[2];
  MPI_Cart_coords(comm, rank, 2, coords);
  printf("\nrank = %d, coords = (%d, %d)\n", rank, coords[0], coords[1]);

  int sub_matrix_size = atoi(argv[1]);
  int **a = initialize_matrix(sub_matrix_size);
  int **b = initialize_matrix(sub_matrix_size);
  int **c = initialize_matrix(sub_matrix_size);
  int **diag = initialize_matrix(sub_matrix_size);

  for (int i = 0; i < sub_matrix_size; i++){
    a[i][i] = rank;
    b[0][i] = rank;
  }

  for (int k = 0; k < root; k++){
    if ((root + coords[1] - coords[0])%root == k){
      diag = a;
    }
    // MPI_Bcast(&diag, root, MPI_INT, ?, new_comm);
    if (k == 0){
      printf("k = %d\n", k);
      printf("rank = %d\n", rank);
      printf("diag\n");
      print_matrix(diag, sub_matrix_size);
    }
  }

  free(a);
  free(b);
  free(c);
  free(diag);
  MPI_Finalize();
  return 0;
}
