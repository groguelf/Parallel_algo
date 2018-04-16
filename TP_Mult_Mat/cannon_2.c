#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void fill_matrix(int *, int);
void print_matrix(int *, int);
void broadcast_row(int, int);
void mul_matrices(int *, int *, int *, int);
void copy_A(int *, int *, int);
void skewing(int *, int , int , MPI_Comm, int); 

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
int source, dest;

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

	int matrix_dimension = atoi(argv[1]);
	if (matrix_dimension % 2 != 0) {
		printf("Matrix matrix dimension should be even\n");
		exit(-2);	
	}

	if ((matrix_dimension % q) != 0) {
		printf("the square root of number of Processes should devide matrix dimension\n");
		exit(-2);
	}

	int matix_size = matrix_dimension*matrix_dimension;
	int sub_matrix_size = (matrix_dimension*matrix_dimension)/p;
	int sub_matrix_dim = sqrt(sub_matrix_size);

	A = (int *) malloc(sizeof(int) * matix_size);
	B = (int *) malloc(sizeof(int) * matix_size);
	C = (int *) malloc(sizeof(int) * matix_size);

	a = (int *) malloc(sizeof(int) * sub_matrix_size);
	b = (int *) malloc(sizeof(int) * sub_matrix_size);
	c = (int *) malloc(sizeof(int) * sub_matrix_size);

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
		fill_matrix(A, matix_size);
		fill_matrix(B, matix_size);
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

    for (i =0; i < sub_matrix_size; i++){
    	c[i] = 0;
    }

    double start, end;
    start = MPI_Wtime();

    MPI_Scatterv(A, sendcounts, displs, type, a, sub_matrix_size, MPI_INT, 0, comm_cart);
    MPI_Scatterv(B, sendcounts, displs, type, b, sub_matrix_size, MPI_INT, 0, comm_cart);

	// MPI_Scatter(A, MATRIX_SIZE/p, MPI_INT, a, MATRIX_SIZE/p, MPI_INT,0,comm_cart);
	// MPI_Scatter(A, MATRIX_SIZE/p, MPI_INT, a, MATRIX_SIZE/p, MPI_INT,0,comm_cart);

	// if (my_rank == 1)
	// 	print_matrix(a, 4);


		/* Pre-skewing */
	skewing(a, 0, -1 * coords[0], hor_comm, sub_matrix_size);
	skewing(b, 1, -1 * coords[1], ver_comm, sub_matrix_size);


	for (k = 0; k < q; k++) {
		mul_matrices(a, b, c, sub_matrix_dim);
		/*Horizontal shift of A*/
		MPI_Cart_shift(hor_comm, 0, -1, &source, &dest); 
		MPI_Sendrecv_replace(a, sub_matrix_size, MPI_INT, dest, 0,source, 0, hor_comm, &status);
		/*Vertical shift of B*/
		MPI_Cart_shift(ver_comm, 0, -1, &source, &dest); 
		MPI_Sendrecv_replace(b, sub_matrix_size, MPI_INT, dest, 0,source, 0, ver_comm, &status);		
				
	}

	/* post-skewing */
	skewing(a, 0, coords[0], hor_comm, sub_matrix_size);
	skewing(b, 1, coords[1], ver_comm, sub_matrix_size);

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

	
	if (my_rank == 0){
		//print_matrix(C, matrix_dimension);
        end = MPI_Wtime();
        printf("\nRealized in %f s.\n", end - start);
	}
	free(A);
	free(B);
	free(C);
	free(a);
	free(b);
	free(c);
	MPI_Finalize();

	return 0;
}

void skewing(int *matrix, int direction, int steps, MPI_Comm comm, int sub_matrix_size) 
{
	MPI_Cart_shift(comm, direction, steps, &source, &dest); 
	MPI_Sendrecv_replace(matrix, sub_matrix_size, MPI_INT, dest, 0,source, 0, comm, &status);
}

void mul_matrices(int *a, int *b, int *c, int matrix_dimension)
{
    int i, j, k;
    for (i = 0; i < matrix_dimension; i++)
    {
        for (j = 0; j < matrix_dimension; j++)
        {
            for (k = 0; k < matrix_dimension; k++)
                c[i*matrix_dimension + j] += a[(i*matrix_dimension) + k]* b[(k*matrix_dimension) + j];
        }
    }
}

void broadcast_row(int k, int sub_matrix_size)
{
	int root = (coords[0]+k) % q;
	MPI_Bcast(
	    a,
	    sub_matrix_size,
	    MPI_INT,
	    root,
	    hor_comm
	);
}

void fill_matrix(int *matrix, int matrix_size) {
	int i;
	for (i = 0; i < matrix_size; i++){
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

void copy_A(int *from, int *to, int matrix_size)
{
	int i;
	for (i = 0; i < matrix_size; i++) {
		to[i] = from[i];
	}
}
