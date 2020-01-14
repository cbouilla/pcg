#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "pcg_variants.h"

#define DEBUG 1

int main()
{
    // Seed the full-blown generator with external entropy
    pcg128_t seeds[2];
    
    FILE *f = fopen("/dev/urandom", "r");
    if (f == NULL) {
        perror("Something went wrong when opening /dev/urandom");
        exit(EXIT_FAILURE);
    }
    if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
        perror("Something went wrong when reading /dev/urandom");
        exit(EXIT_FAILURE);
    }

    printf("Seed[0] : %016" PRIx64 " %016" PRIx64 "\n", (uint64_t) (seeds[0] >> 64), (uint64_t) seeds[0]);
    printf("Seed[1] : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seeds[1] >> 64), (uint64_t) seeds[1]);

    pcg64_random_t rng;
    pcg64_srandom_r(&rng, seeds[0], seeds[1]);
    
    printf("Predictor input:\n");
    for (int i = 0; i < 32; i++) {
        printf("X[%2d] = 0x%016" PRIx64 "\n", i, pcg64_random_r(&rng));
        if (i == 0 && DEBUG) {
            int known_low = 13;
            printf("known_low = %d\n", known_low);
            printf("W_0 = %04lld\n", (unsigned long long) (rng.state % (1 << known_low)));
            printf("W_c = %04lld\n", (unsigned long long) (rng.inc % (1 << known_low)));
        }
    }
    printf("\n");
    printf("Remaining of the sequence (predictor output, in principle):\n");
    for (int i = 32; i < 48; i++)
        printf("X[%2d] = 0x%016" PRIx64 "\n", i, pcg64_random_r(&rng));

    return EXIT_SUCCESS;
}