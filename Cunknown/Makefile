CFLAGS = -O3 -Wall -Wextra -Werror -march=native -mtune=native -fopenmp
LDFLAGS = -fopenmp
LDLIBS = -lgmp -lm

all: main test

fonctions.o: fonctions.h
main.o: fonctions.h pcg_oneseq.h
test.o: fonctions.h pcg_oneseq.h
main: fonctions.o pcg_oneseq.o main.o
test: fonctions.o pcg_oneseq.o test.o
test_falsenegative.o: fonctions.h pcg_oneseq.h
test_falsenegative: fonctions.o pcg_oneseq.o test_falsenegative.o

.PHONY: clean check

clean:
	rm -rf *.o main test test_falsenegative

check: test test_falsenegative
	prove -v ./test ./test_falsenegative