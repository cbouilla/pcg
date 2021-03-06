#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "pcg_variants.h"

int main()
{
	pcg128_t seed = (((pcg128_t) 0x5b3d418b8b51364b) << 64) + ((pcg128_t) 0x74eed710223e763f);
	printf("Seed : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seed >> 64), (uint64_t) seed);

	pcg64s_random_t rng;
	pcg64s_srandom_r(&rng, seed);

	printf("Predictor input:\n");
	for (int i = 0; i < 3; i++)
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64s_random_r(&rng));
	printf("\n");

	printf("Remaining of the sequence :\n");
	for (int i = 3; i < 16; i++)
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64s_random_r(&rng));

	return EXIT_SUCCESS;
}
