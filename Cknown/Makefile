# icc (avec "source /usr/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64" avant)
#CC = icc
#CFLAGS = -O3 -qopenmp
#LDFLAGS = -qopenmp

CC = gcc
CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp



LDLIBS = -lm

all: main

main.o: fonctions.h pcg_oneseq.h
main: fonctions.o pcg_oneseq.o main.o

.PHONY: clean

clean:
	rm -rf *.o main
