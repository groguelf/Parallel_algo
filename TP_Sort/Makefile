all: triangular_matrix bubble qsort mergesort

triangular_matrix: triangular_matrix.o
	gcc -fopenmp -o triangular_matrix triangular_matrix.o

bubble: bubble.o
	gcc -fopenmp -o bubble bubble.o

qsort: qsort.o
	gcc -fopenmp -o qsort qsort.o -lm

mergesort: mergesort.o
	gcc -fopenmp -o mergesort mergesort.o

%.o: %.c
	gcc -c -O2 -fopenmp $^

clean:
	rm -f triangular_matrix bubble qsort mergesort *.o *~
