#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

__uint128_t state;

uint64_t lehmer64() {
  state *= 0xda942042e4dd58b5;
  return state >> 64;
}

int main()
{
	FILE *f = fopen("/dev/urandom", "r");
	if (f == NULL) {
		perror("Something went wrong when opening /dev/urandom");
		exit(EXIT_FAILURE);
	}
	if (fread(&state, sizeof(state), 1, f) != 1) {
		perror("Something went wrong when reading /dev/urandom");
		exit(EXIT_FAILURE);
	}
	state |= 1;
	printf("Seed : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (state >> 64), (uint64_t) state);

	printf("Predictor input:\n");
	for (int i = 0; i < 3; i++)
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, lehmer64());
	printf("\n");

	printf("Remaining of the sequence (predictor output, in principle):\n");
	for (int i = 3; i < 16; i++)
		printf("X[%2d] = 0x%016" PRIx64 ";\n", i, lehmer64());

	return EXIT_SUCCESS;
}