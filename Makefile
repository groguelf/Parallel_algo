all: triangular_matrix bubble

triangular_matrix: triangular_matrix.o
	gcc -fopenmp -o triangular_matrix triangular_matrix.o

bubble: bubble.o
	gcc -fopenmp -o bubble bubble.o

%.o: %.c
	gcc -c -O0 -fopenmp $^

clean:
	rm -f triangle_matrix bubble *.o *~
