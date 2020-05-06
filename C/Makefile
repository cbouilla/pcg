CFLAGS = -O3 -Wall -Wextra -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp -g
LDLIBS = -lgmp -lm

all: main

main.o: fonctions.h pcg_oneseq.h
main: fonctions.o pcg_oneseq.o main.o

.PHONY: clean

clean:
	rm -rf *.o main
