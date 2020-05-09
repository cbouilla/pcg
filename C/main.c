#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// avec ICC, ce truc ne marche pas sur ppti-gpu-1 si on met -xHost (zmm?). Il faut investiguer...

void result_found(pcg128_t S, double start)
{
	printf("\nInternal state found (%.1fs)\n", wtime() - start);
	
	struct pcg_state_128 rng;
	rng.state = (S - c) * a_inv;
	pcg128_t seed = (rng.state - c - a*c) * a_inv;

	printf("Seed : %016" PRIx64 " %016" PRIx64 "\n", (uint64_t) (seed >> 64), (uint64_t) seed);
	for(int i = 0 ; i < 16 ; i++)
		printf("X[%2d] = %016" PRIx64 "\n", i, pcg_oneseq_128_xsl_rr_64_random_r(&rng));
}


int main()
{
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
	
    /********** INPUT ***********/
    unsigned long long X[nbiter];

    // challenge output given by M. O'Neil
    X[ 0] = 0x6ec191a37a421087;
    X[ 1] = 0xec140ace169176fc;
    X[ 2] = 0x85994d489913af70;

    double start = wtime();
    unsigned long long done = 0;
    
    #pragma omp parallel for
    for (unsigned long long W0 = 0; W0 < (1<<known_low) ; W0++){
	
	/*Variables privées*/
	pcg128_t S[nbiter];
	// u64 urX[nbiter];
	int rot[nbiter] = {0, 0, 0};
	struct task_t task;
	setup_task(W0, X, &task);

	if (omp_get_thread_num() == 0 && (W0 % 64) == 0) {
	    printf("\rW0 = %llx / %llx --- %.1f / s", done, 1ull<<known_low, done / (wtime() - start));
	    fflush(stdout);
	}


	for (int r = 0 ; r < 1<<(3*known_up) ; r++) {
	    /***** Modification de rot et unrotX *****/
	    rot[0] = (rot[0] + 1) % k;
	    int i = 0;
	    while (rot[i] == 0 && i < nbiter) {
		i++;
		rot[i] = (rot[i] + 1) %k;
	    }
	    // unrotate(urX, X, rot);
	    
	    /***** Résolution *****/
	    if (solve(S, rot, &task))
		result_found(S[0], start);
	}
	
	#pragma omp atomic
	done++;;
    }
    printf("Total time = %.1f\n", wtime() - start);
    return EXIT_SUCCESS;
}