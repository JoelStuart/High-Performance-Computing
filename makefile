CC=mpiicc
CFLAGS=-I. -lm -O2

%.o: %.c 
	$(CC) -c -fopenmp -o $@ $< $(CFLAGS)

proj: proj.o 
	mpiicc -fopenmp -o proj proj.o -I. -lm -O2