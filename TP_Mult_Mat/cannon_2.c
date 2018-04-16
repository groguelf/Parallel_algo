#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define MATRIX_SIZE 64

void fill_matrix(int *, int);
void print_matrix(int *, int);
void broadcast_row(int);
void mul_matrices(int *, int *, int *);
void copy_A(int *, int *);

MPI_Status status;
MPI_Comm comm_cart;
MPI_Comm hor_comm;
MPI_Comm ver_comm;

int ndims, dims[2], period[2], reorder, coords[2], hor_dims[2], ver_dims[2];

int *A;
int *B;
int *C;
int *a;
int *old_a;
int *b;
int *c;
int my_rank, p, q;

int main(int argc, char *argv[])
{
	int i, j, k;
	int source, dest;	
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	q = sqrt(p);
	if (q*q != p) {
		printf("Number of Processes must be a perfect square\n");
		exit(-1);
	}

	int sub_matrix_size = MATRIX_SIZE/p;
	int sub_matrix_dim = sqrt(sub_matrix_size);

	int matrix_dimension = sqrt(MATRIX_SIZE);
	A = (int *) malloc(sizeof(int) * MATRIX_SIZE);
	B = (int *) malloc(sizeof(int) * MATRIX_SIZE);
	C = (int *) malloc(sizeof(int) * MATRIX_SIZE);

	a = (int *) malloc(sizeof(int) * MATRIX_SIZE/p);
	old_a = (int *) malloc(sizeof(int) * MATRIX_SIZE/p);
	b = (int *) malloc(sizeof(int) * MATRIX_SIZE/p);
	c = (int *) malloc(sizeof(int) * MATRIX_SIZE/p);

	ndims = 2;
	period[0] = 1; period[1] = 1;
	reorder = 0;
	hor_dims[0] = 0; hor_dims[1] = 1;
	ver_dims[0] = 1; ver_dims[1] = 0;

	MPI_Dims_create(p, ndims, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &comm_cart);
	MPI_Cart_coords(comm_cart, my_rank, 2, coords);
	MPI_Cart_sub(comm_cart, hor_dims, &hor_comm);
	MPI_Cart_sub(comm_cart, ver_dims, &ver_comm);

	if (my_rank == 0){
		fill_matrix(A, 0);
		fill_matrix(B, 0);
	}


	MPI_Datatype type, temp_type;
	MPI_Type_vector(sub_matrix_dim, sub_matrix_dim, matrix_dimension, MPI_INT, &temp_type);
    MPI_Type_create_resized(temp_type, 0, sizeof(int), &type);
    MPI_Type_commit(&type);
    
    int sendcounts[p];
    for(i=0; i<p; i++)
    	sendcounts[i] = 1;

    int displs[p];
    for (i=0; i<q; i++) {
        for (j=0; j < q; j++) {
        	displs[i*q + j] = i * sub_matrix_size * q + j*sub_matrix_dim;
        }
    }

    MPI_Scatterv(A, sendcounts, displs, type, a, sub_matrix_size, MPI_INT, 0, comm_cart);
    MPI_Scatterv(B, sendcounts, displs, type, b, sub_matrix_size, MPI_INT, 0, comm_cart);

	// MPI_Scatter(A, MATRIX_SIZE/p, MPI_INT, a, MATRIX_SIZE/p, MPI_INT,0,comm_cart);
	// MPI_Scatter(A, MATRIX_SIZE/p, MPI_INT, a, MATRIX_SIZE/p, MPI_INT,0,comm_cart);

	// if (my_rank == 1)
	// 	print_matrix(a, 4);

    for (i =0; i < sub_matrix_size; i++){
    	c[i] = 0;
    }

	for (k = 0; k < q; k++) {
		copy_A(a, old_a);
		broadcast_row(k);
		mul_matrices(a,b,c);
		MPI_Cart_shift(ver_comm, coords[1], -1, &source, &dest); 
		MPI_Sendrecv_replace(b, sub_matrix_size, MPI_INT, dest, 0,source, 0, ver_comm, &status);
		copy_A(old_a, a);
		MPI_Barrier(comm_cart);
	} 
	MPI_Barrier(comm_cart);

	MPI_Gatherv(
	    c,
	    sub_matrix_size,
	    MPI_INT,
	    C,
	    sendcounts,
	    displs,
	    type,
	    0,
	    comm_cart
	);

	if (my_rank == 0)
		print_matrix(C,8);
	// free(A);
	// free(a);
	MPI_Finalize();

	return 0;
}

void mul_matrices(int *a, int *b, int *c)
{
    int i, j, k;
    int matrix_dimension = sqrt(MATRIX_SIZE/p);
    for (i = 0; i < matrix_dimension; i++)
    {
        for (j = 0; j < matrix_dimension; j++)
        {
            for (k = 0; k < matrix_dimension; k++)
                c[i*matrix_dimension + j] += a[(i*matrix_dimension) + k]* b[(k*matrix_dimension) + j];
        }
    }
}

void broadcast_row(int k)
{
	int root = (coords[0]+k) % q;
	MPI_Bcast(
	    a,
	    MATRIX_SIZE/p,
	    MPI_INT,
	    root,
	    hor_comm
	);
}

void fill_matrix(int *matrix, int my_rank) {
	int i,j;
	int matrix_dimension = sqrt(MATRIX_SIZE); 
	for (i = 0; i < MATRIX_SIZE; i++){
		matrix[i] = i;
	} 
}

void print_matrix(int *matrix, int matrix_dimension)
{
	int i, j;
	for (i = 0; i < matrix_dimension; i++) {
		for (j = 0; j <matrix_dimension; j++)
			printf("%d ", matrix[i*matrix_dimension + j]);

		printf("\n");
	}
}

void copy_A(int *from, int *to)
{
	int i;
	for (i = 0; i < MATRIX_SIZE/p; i++) {
		to[i] = from[i];
	}
}