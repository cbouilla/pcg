#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "pcg_variants.h"

int main()
{
    // Seed the "single-stream" generator with external entropy
    pcg128_t seed;    
    FILE *f = fopen("/dev/urandom", "r");
    if (f == NULL) {
        perror("Something went wrong when opening /dev/urandom");
        exit(EXIT_FAILURE);
    }
    if (fread(&seed, sizeof(seed), 1, f) != 1)  {
        perror("Something went wrong when reading /dev/urandom");
        exit(EXIT_FAILURE);
    }
    printf("Seed : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seed >> 64), (uint64_t) seed);

    pcg64s_random_t rng;
    pcg64s_srandom_r(&rng, seed);
    
    printf("Predictor input:\n");
    for (int i = 0; i < 3; i++)
        printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64s_random_r(&rng));
    printf("\n");
    
    printf("Remaining of the sequence (predictor output, in principle):\n");
    for (int i = 3; i < 10; i++)
        printf("X[%2d] = 0x%016" PRIx64 ";\n", i, pcg64s_random_r(&rng));

    return EXIT_SUCCESS;
}