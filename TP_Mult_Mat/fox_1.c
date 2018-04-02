#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void fox(int *, int);
void fill_matrix(int *, int, int );
void print_matrix(int *);
void broadcast_row(int);
void mul_matrices(int *, int *, int *);
void copy_A(int *, int *);

MPI_Status status;
MPI_Comm comm_cart;
MPI_Comm hor_comm;
MPI_Comm ver_comm;

int nnodes, ndims, dims[2], period[2], reorder, coords[2], hor_dims[2], ver_dims[2];


int *A;
int *old_A;
int *B;
int *C;
int my_rank, p, q;

int main(int argc, char *argv[])
{
	
	int i, k;
	int source, dest;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	A = (int *) malloc(sizeof(int) * p);
	old_A = (int *) malloc(sizeof(int) * p);
	B = (int *) malloc(sizeof(int) * p);
	C = (int *) malloc(sizeof(int) * p);

	nnodes = p;
	q = sqrt(p);
	ndims = 2;
	period[0] = 1; period[1] = 1;
	reorder = 0;
	hor_dims[0] = 0; hor_dims[1] = 1;
	ver_dims[0] = 1; ver_dims[1] = 0;

	MPI_Dims_create(nnodes, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &comm_cart);
	MPI_Cart_coords(comm_cart, my_rank, 2, coords);
	MPI_Cart_sub(comm_cart, hor_dims, &hor_comm);
	MPI_Cart_sub(comm_cart, ver_dims, &ver_comm);


	fill_matrix(A, p, my_rank);
	fill_matrix(B, p, my_rank);

	for (i=0; i< p; i++) {
		C[i] = 0;
	}

	for (k = 0; k < q; k++) {
		copy_A(A, old_A);
		broadcast_row(k);
		mul_matrices(A,B,C);
		MPI_Cart_shift(ver_comm, coords[1], 1, &source, &dest); 
		MPI_Sendrecv_replace(B, p, MPI_INT, dest, 0,source, 0, ver_comm, &status);
		copy_A(old_A, A);
	}
	
	

	if (my_rank == 2){
		print_matrix(C);
	}
		
	// if (my_rank == 0) 
		// printf("Finished\n" );
	MPI_Finalize();

	return 0;
}

void mul_matrices(int *A, int *B, int *C)
{
    int i, j, k;
    int matrix_dimension = sqrt(p);
    for (i = 0; i < matrix_dimension; i++)
    {
        for (j = 0; j < matrix_dimension; j++)
        {
            // C[(i*matrix_dimension) + j] = 0;
            for (k = 0; k < matrix_dimension; k++)
                C[i*matrix_dimension + j] += A[(i*matrix_dimension) + k]* B[(k*matrix_dimension) + j];
        }
    }
}

void broadcast_row(int k)
{
	int root = (coords[0]+k) % q;
	if (coords[0] == ((coords[1]+k) % q)) {
		MPI_Bcast(
		    A,
		    p,
		    MPI_INT,
		    root,
		    hor_comm
    	);
	} else {
		MPI_Bcast(
		    A,
		    p,
		    MPI_INT,
		    root,
		    hor_comm
    	);
	}

}

void fill_matrix(int *matrix, int size, int my_rank) {
	int i;
	for (i = 0; i < size; i++) {
		matrix[i] = my_rank*size + i;
	}
}

void print_matrix(int *matrix)
{
	int i, j;
	for (i = 0; i < q; i++) {
		for (j = 0; j <q; j++)
			printf("%d ", matrix[i*q + j]);

		printf("\n");
	}
}

void copy_A(int *from, int *to)
{
	int i;
	for (i = 0; i < p; i++) {
		to[i] = from[i];
	}
}