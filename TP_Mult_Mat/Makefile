CC=mpicc
CFLAGS=-g -Wall -Wextra -std=c99
LFLAGS=-lm

PROGRAM=exercise_1 fox_1 cannon

all: $(PROGRAM)

exercise_1:
	$(CC) $(CFLAGS) -o exercise_1 exercise_1.c $(LFLAGS)

fox_1:
	$(CC) $(CFLAGS) -o fox_1 fox_1.c $(LFLAGS)

cannon:
	$(CC) $(CFLAGS) -o cannon cannon.c $(LFLAGS)

clean:
	rm -f $(PROGRAM) *~
