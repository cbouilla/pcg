# icc (avec "source /usr/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64" avant)
#CC = icc
#CFLAGS = -O3 -xHost -qopt-zmm-usage=high -qopenmp
#LDFLAGS = -fopenmp -qopenmp

# gcc
CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp
LDLIBS = -lm

all: main benchmark test test_falsenegative benchmark_omp

main.o: CC=mpicc
main: CC=mpicc	

fonctions.o: fonctions.h
fonctions_bonus.o: fonctions.h
main.o: fonctions.h pcg_setseq.h
test.o: fonctions.h pcg_setseq.h
test_falsenegative.o: fonctions.h pcg_setseq.h
benchmark.o: fonctions.h pcg_setseq.h
benchmark_omp.o: fonctions.h pcg_setseq.h

main: fonctions.o main.o
benchmark: fonctions.o benchmark.o
benchmark_omp: fonctions.o benchmark_omp.o
test: fonctions.o test.o fonctions_bonus.o
test_falsenegative: fonctions.o test_falsenegative.o

.PHONY: clean check

clean:
	rm -rf *.o main test test_falsenegative benchmark

check: test test_falsenegative
	prove -v ./test ./test_falsenegative

bench: benchmark benchmark_omp
	./benchmark
	./benchmark_omp