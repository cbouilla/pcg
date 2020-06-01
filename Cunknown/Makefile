# icc (avec "source /usr/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64" avant)
#CC = icc
#CFLAGS = -O3 -xHost -qopt-zmm-usage=high -qopenmp
#LDFLAGS = -fopenmp -qopenmp

# gcc
CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp
LDLIBS = -lm

all: main pre_validator success_proba

main.o: CC=mpicc
main: CC=mpicc	

fonctions.o: fonctions.h
fonctions_bonus.o: fonctions.h
main.o: fonctions.h pcg_setseq.h
success_proba.o: fonctions.h pcg_setseq.h

main: fonctions.o main.o
success_proba: fonctions.o success_proba.o
pre_validator: fonctions.o pre_validator.o

.PHONY: clean check

clean:
	rm -rf *.o main test test_falsenegative benchmark benchmark_omp pre_validator

