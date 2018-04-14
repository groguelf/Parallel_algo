#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define SUB_SIZE 4

void fill_matrix(int *, int );
void print_matrix(int *);
void broadcast_row(int);
void mul_matrices(int *, int *, int *);
void copy_A(int *, int *);

MPI_Status status;
MPI_Comm comm_cart;
MPI_Comm hor_comm;
MPI_Comm ver_comm;

int ndims, dims[2], period[2], reorder, coords[2], hor_dims[2], ver_dims[2];


int *a;
int *old_a;
int *b;
int *c;
int my_rank, p, q;

int main(int argc, char *argv[])
{
	int i, k;
	int source, dest;	
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	q = sqrt(p);
	if (q*q != p) {
		printf("Number of Processes must be a perfect square\n");
		exit(-1);
	}

	a = (int *) malloc(sizeof(int) * SUB_SIZE);
	old_a = (int *) malloc(sizeof(int) * SUB_SIZE);
	b = (int *) malloc(sizeof(int) * SUB_SIZE);
	c = (int *) malloc(sizeof(int) * SUB_SIZE);

	
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


	fill_matrix(a, my_rank);
	fill_matrix(b, my_rank);

	for (i=0; i< SUB_SIZE; i++) {
		c[i] = 0;
	}

	for (k = 0; k < q; k++) {
		copy_A(a, old_a);
		broadcast_row(k);
		mul_matrices(a,b,c);
		MPI_Cart_shift(ver_comm, coords[1], -1, &source, &dest); 
		MPI_Sendrecv_replace(b, SUB_SIZE, MPI_INT, dest, 0,source, 0, ver_comm, &status);
		copy_A(old_a, a);
	}

	if (my_rank == 0)
		print_matrix(c); 
		
	MPI_Finalize();

	return 0;
}

void mul_matrices(int *a, int *b, int *c)
{
    int i, j, k;
    int matrix_dimension = sqrt(SUB_SIZE);
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
	    SUB_SIZE,
	    MPI_INT,
	    root,
	    hor_comm
	);
}

void fill_matrix(int *matrix, int my_rank) {
	int i;
	for (i = 0; i < SUB_SIZE; i++) {
		matrix[i] = my_rank*SUB_SIZE + i;
	}
}

void print_matrix(int *matrix)
{
	int i, j;
	int matrix_dimension = sqrt(SUB_SIZE);
	for (i = 0; i < matrix_dimension; i++) {
		for (j = 0; j <matrix_dimension; j++)
			printf("%d ", matrix[i*matrix_dimension + j]);

		printf("\n");
	}
}

void copy_A(int *from, int *to)
{
	int i;
	for (i = 0; i < SUB_SIZE; i++) {
		to[i] = from[i];
	}
}