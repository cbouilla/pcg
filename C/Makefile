# CC = icc
# CFLAGS = -O3 -xHost -qopt-zmm-usage=high -qopenmp
# LDFLAGS = -fopenmp -qopenmp

# gcc
CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp



LDLIBS = -lm

all: main

main.o: fonctions.h pcg_oneseq.h
main: fonctions.o pcg_oneseq.o main.o

.PHONY: clean

clean:
	rm -rf *.o main
