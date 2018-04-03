#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int **initialize_matrix(int matrix_size){
  int **mat = calloc(matrix_size, sizeof(int *));
  for (int i = 0; i < matrix_size; i++){
    mat[i] = calloc(matrix_size, sizeof(int));
  }
  return mat;
}

void free_matrix(int **mat, int matrix_size){
  for (int i = 0; i < matrix_size; i++){
    free(mat[i]);
  }
  free(mat);
}

void mul_matrices(int **A, int **B, int **C, int matrix_size){
  for (int i = 0; i < matrix_size; i++){
    for (int j = 0; j < matrix_size; j++){
      for (int k = 0; k < matrix_size; k++){
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void copy(int **from, int **to, int matrix_size){
	for (int i = 0; i < matrix_size; i++){
    for (int j = 0; j < matrix_size; j++){
      to[i][j] = from[i][j];
    }
	}
}


void print_matrix(int **matrix, int matrix_size){
  for (int i = 0; i < matrix_size; i++){
    for(int j = 0; j < matrix_size; j++){
      fprintf(stderr, "%d     ", matrix[i][j]);
    }
    fprintf(stderr, "\n");
  }
}

int main(int argc, char **argv){
  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int root = (int) sqrt(size);
  if (root * root != size){
    if (rank == 0){
      fprintf(stderr, "You have to use a perfect square number of processes.\n\n");
    }
    MPI_Finalize();
    return 0;
  }

  int nnodes = size;
  int ndims = 2;
  int dims[ndims];
  int periods[ndims];
  int reorder = 0;
  MPI_Comm comm, hor_comm, ver_comm;
  int hor_dims[ndims], ver_dims[ndims];

  for (int i = 0; i < ndims; i++){
    dims[i] = root;
    periods[i] = 0;
  }
  hor_dims[0] = 0; hor_dims[1] = 1;
  ver_dims[0] = 1; ver_dims[1] = 0;
  int source, dest;
  MPI_Status status;

  MPI_Dims_create(nnodes, ndims, dims);
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,periods, reorder, &comm);
  MPI_Cart_sub(comm, hor_dims, &hor_comm);
  MPI_Cart_sub(comm, ver_dims, &ver_comm);

  int coords[2];
  MPI_Cart_coords(comm, rank, 2, coords);
  printf("\nrank = %d, coords = (%d, %d)\n", rank, coords[0], coords[1]);

  int sub_matrix_size = atoi(argv[1]);
  int **A = initialize_matrix(sub_matrix_size);
  int **B = initialize_matrix(sub_matrix_size);
  int **C = initialize_matrix(sub_matrix_size);
  int **old_A = initialize_matrix(sub_matrix_size);

  for (int i = 0; i < sub_matrix_size; i++){
    A[i][i] = rank;
    B[0][i] = rank;
  }

  for (int k = 0; k < root; k++){
    copy(A, old_A, sub_matrix_size);
    fprintf(stderr, "\nk = %d, rank = %d and copy done\n", k, rank);
    MPI_Bcast(A, size, MPI_INT, (coords[0]+k)%root, hor_comm);
    fprintf(stderr, "\nk = %d, rank = %d and cast done\n", k, rank);
		mul_matrices(A, B, C, sub_matrix_size);
    fprintf(stderr, "\nk = %d, rank = %d and mult done\n", k, rank);
		MPI_Cart_shift(ver_comm, coords[1], 1, &source, &dest);
		MPI_Sendrecv_replace(B, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, ver_comm, &status);
		copy(old_A, A, sub_matrix_size);
  }

  fprintf(stderr, "rank = %d\n", rank);
  print_matrix(C, sub_matrix_size);

  free_matrix(A, sub_matrix_size);
  free_matrix(B, sub_matrix_size);
  free_matrix(C, sub_matrix_size);
  free_matrix(old_A, sub_matrix_size);
  MPI_Finalize();
  return 0;
}
