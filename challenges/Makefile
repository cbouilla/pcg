PCG = $(HOME)/prison/pcg-c-0.94
CPPFLAGS += -I$(PCG)/include
CFLAGS   += -std=c99 -Wall -Wextra
LDFLAGS  += -L$(PCG)/src
LDLIBS   += -lpcg_random

all: pcg64-challenge pcg64s-challenge

clean:
	rm -f *.o pcg64-challenge pcg64s-challenge
