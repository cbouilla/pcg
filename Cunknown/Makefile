CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp -g
LDFLAGS = -fopenmp
LDLIBS = -lm

all: main benchmark test

fonctions.o: fonctions.h
main.o: fonctions.h pcg_setseq.h
test.o: fonctions.h pcg_setseq.h
test_falsenegative.o: fonctions.h pcg_setseq.h
benchmark.o: fonctions.h pcg_setseq.h

main: fonctions.o main.o
benchmark: fonctions.o benchmark.o
test: fonctions.o test.o fonctions_bonus.o
test_falsenegative: fonctions.o test_falsenegative.o

.PHONY: clean check

clean:
	rm -rf *.o main test test_falsenegative benchmark

check: test test_falsenegative
	prove -v ./test ./test_falsenegative