#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "pcg_variants.h"

int main()
{ 
	pcg128_t seed = PCG_128BIT_CONSTANT(0xffbb5e70239d0d0dULL, 0x5a21d22cdb7dfb4cULL);
	pcg128_t increment = PCG_128BIT_CONSTANT(0x3d2204fa591b526aULL, 0x07cc1987323f4721ULL);
 
	printf("increment : %016" PRIx64 " %016" PRIx64 "\n", (uint64_t) (increment >> 64), (uint64_t) increment);
	printf("Seed      : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seed >> 64), (uint64_t) seed);

	pcg64_random_t rng;
	pcg64_srandom_r(&rng, seed, increment);

	printf("Predictor input:\n");
	for (int i = 0; i < 48; i++) {
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64_random_r(&rng));
	}
	printf("\n");
	printf("Remaining of the sequence:\n");
	for (int i = 48; i < 64; i++)
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64_random_r(&rng));

	return EXIT_SUCCESS;
}
