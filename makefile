OPT  = -O3 -w 

all: clean aristas
aristas: aristas.c aux.o 
	mpicc $(OPT) -o aristas aristas_mpi.c aux.o
aux.o: aux.c
	gcc $(OPT) -c aux.c
clean:
	rm -f *.o aristas
