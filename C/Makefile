CFLAGS = -O3 -Wall -Wextra -march=native -mtune=native
LDLIBS = -lgmp -lm

all: main

main: fonctions.o pcg_oneseq.o main.o

.PHONY: clean

clean:
	rm -rf *.o main
