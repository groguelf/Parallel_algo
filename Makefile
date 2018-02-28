all: triangular_matrix bubble qsort

triangular_matrix: triangular_matrix.o
	gcc -fopenmp -o triangular_matrix triangular_matrix.o

bubble: bubble.o
	gcc -fopenmp -o bubble bubble.o

qsort: qsort.o
	gcc -fopenmp -o qsort qsort.o

%.o: %.c
	gcc -c -O0 -fopenmp $^

clean:
	rm -f triangular_matrix bubble *.o *~
