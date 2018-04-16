#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void mul_matrices(int *A, int *B, int *C, int matrix_size){
    for (int i = 0; i < matrix_size; i++){
        for (int j = 0; j < matrix_size; j++){
            for (int k = 0; k < matrix_size; k++){
                C[i*matrix_size + j] += A[i*matrix_size + k] * B[k*matrix_size + j];
            }
        }
    }
}

void print_matrix(int *matrix, int matrix_size){
    for (int i = 0; i < matrix_size; i++){
        for(int j = 0; j < matrix_size; j++){
            printf("%d     ", matrix[i*matrix_size + j]);
        }
        printf("\n");
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
    MPI_Comm comm, hor_comm, ver_comm;
    int hor_dims[ndims], ver_dims[ndims];

    for (int i = 0; i < ndims; i++){
        dims[i] = root;
        periods[i] = 1;
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

    int sub_matrix_size = atoi(argv[1]);
    int *A = calloc(sub_matrix_size*sub_matrix_size, sizeof(int));
    int *B = calloc(sub_matrix_size*sub_matrix_size, sizeof(int));
    int *C = calloc(sub_matrix_size*sub_matrix_size, sizeof(int));
    int *old_A = calloc(sub_matrix_size*sub_matrix_size, sizeof(int));

    // Fill A and B with some values
    // You can put what you want
    for (int i = 0; i < sub_matrix_size; i++){
        for (int j = 0; j < sub_matrix_size; j++){
            A[sub_matrix_size*i + j] = rank;
        }
        B[i] = rank;
    }

    double start, end;
    start = MPI_Wtime();

	// Preskewing of A    
    MPI_Cart_shift(hor_comm, coords[0], -coords[0], &source, &dest);
    MPI_Sendrecv_replace(A, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, hor_comm, &status);

	// Preskewing of B    
    MPI_Cart_shift(ver_comm, coords[1], -coords[1], &source, &dest);
    MPI_Sendrecv_replace(B, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, ver_comm, &status);
    
    for (int k = 0; k < root; k++){
	    mul_matrices(A, B, C, sub_matrix_size);

        // Horizontal shift of A
        MPI_Cart_shift(hor_comm, coords[0], -1, &source, &dest);
        MPI_Sendrecv_replace(A, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, hor_comm, &status);

        // Vertical shift of B
	    MPI_Cart_shift(ver_comm, coords[1], -1, &source, &dest);
	    MPI_Sendrecv_replace(B, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, ver_comm, &status);
    }

	// Postskewing of A    
    MPI_Cart_shift(hor_comm, coords[0], coords[0], &source, &dest);
    MPI_Sendrecv_replace(A, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, hor_comm, &status);

	// Postskewing of B    
    MPI_Cart_shift(ver_comm, coords[1], coords[1], &source, &dest);
    MPI_Sendrecv_replace(B, sub_matrix_size*sub_matrix_size, MPI_INT, dest, 0, source, 0, ver_comm, &status);

    end = MPI_Wtime();

    //printf("rank = %d\n", rank);
    //print_matrix(C, sub_matrix_size);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        printf("\nRealized in %f s.\n", end - start);
    }

    free(A);
    free(B);
    free(C);
    free(old_A);
    MPI_Finalize();
    return 0;
}
