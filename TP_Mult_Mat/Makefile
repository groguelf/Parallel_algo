CC=mpicc
CFLAGS=-g -Wall -Wextra -std=c99
LFLAGS=-lm

PROGRAM=exercise_1 exercise_2 fox_1 cannon fox_2 cannon_2

all: $(PROGRAM)

exercise_1:
	$(CC) $(CFLAGS) -o exercise_1 exercise_1.c $(LFLAGS)

exercise_2:
	$(CC) $(CFLAGS) -o exercise_2 exercise_2.c $(LFLAGS)

fox_1:
	$(CC) $(CFLAGS) -o fox_1 fox_1.c $(LFLAGS)

fox_2:
	$(CC) $(CFLAGS) -o fox_2 fox_2.c $(LFLAGS)

cannon:
	$(CC) $(CFLAGS) -o cannon cannon.c $(LFLAGS)

cannon_2:
	$(CC) $(CFLAGS) -o cannon_2 cannon_2.c $(LFLAGS)

clean:
	rm -f $(PROGRAM) *~
