CC=mpicc
CFLAGS=-g -Wall -Wextra -std=c99
LFLAGS=-lm

PROGRAM=exercise_1

all: $(PROGRAM)

exercise_1:
	$(CC) $(CFLAGS) -o exercise_1 exercise_1.c $(LFLAGS)

clean:
	rm -f $(PROGRAM) *~
