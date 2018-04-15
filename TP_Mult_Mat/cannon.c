#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define SUB_SIZE 16

void fill_matrix(int *, int );
void print_matrix(int *);
void mul_matrices(int *, int *, int *);
void copy_A(int *, int *);
void skewing(int *, int , int , MPI_Comm ); 

MPI_Status status;
MPI_Comm comm_cart;
MPI_Comm hor_comm;
MPI_Comm ver_comm;

int nnodes, ndims, dims[2], period[2], reorder, coords[2], hor_dims[2], ver_dims[2];


int *a;
int *b;
int *c;
int my_rank, p, q;
int source, dest;

int main(int argc, char *argv[])
{
	
	int i, k;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	q = sqrt(p);
	if (q*q != p) {
		printf("Number of Processes must be a perfect square\n");
		exit(-1);
	}

	a = (int *) malloc(sizeof(int) * SUB_SIZE);
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
	for (i=0; i< p; i++) {
		c[i] = 0;
	}

	/* Pre-skewing */
	skewing(a, 0, -1 * coords[0], hor_comm);
	skewing(b, 1, -1 * coords[1], ver_comm);


	for (k = 0; k < q; k++) {
		mul_matrices(a, b, c);
		/*Horizontal shift of A*/
		MPI_Cart_shift(hor_comm, 0, -1, &source, &dest); 
		MPI_Sendrecv_replace(a, SUB_SIZE, MPI_INT, dest, 0,source, 0, hor_comm, &status);
		/*Vertical shift of B*/
		MPI_Cart_shift(ver_comm, 0, -1, &source, &dest); 
		MPI_Sendrecv_replace(b, SUB_SIZE, MPI_INT, dest, 0,source, 0, ver_comm, &status);		
				
	}

	/* post-skewing */
	skewing(a, 0, coords[0], hor_comm);
	skewing(b, 1, coords[1], ver_comm);

	
	if (my_rank == 0){
		print_matrix(c);
	}
		
	MPI_Finalize();

	return 0;
}

void skewing(int *matrix, int direction, int steps, MPI_Comm comm) 
{
	MPI_Cart_shift(comm, direction, steps, &source, &dest); 
	MPI_Sendrecv_replace(matrix, SUB_SIZE, MPI_INT, dest, 0,source, 0, comm, &status);
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