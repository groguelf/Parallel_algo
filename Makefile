all: triangular_matrix

triangular_matrix: triangular_matrix.o
	gcc -fopenmp -o triangular_matrix triangular_matrix.o

triangular_matrix.o: triangular_matrix.c
	gcc -c -O0 -fopenmp triangular_matrix.c

clean:
	rm -f triangular_matrix triangular_matrix.o *~
