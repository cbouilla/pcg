PCG = $(HOME)/prison/pcg-c-0.94
CPPFLAGS += -I$(PCG)/include
CFLAGS   += -std=c99 -Wall -Wextra
LDFLAGS  += -L$(PCG)/src
LDLIBS   += -lpcg_random

all: pcg64-challenge pcg64s-challenge pcg64s-check pcg64-check

clean:
	rm -f *.o pcg64-challenge pcg64s-challenge pcg64-check pcg64s-check
