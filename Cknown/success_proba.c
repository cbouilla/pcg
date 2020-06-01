#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "fonctions.h"

/* auxiliary script to determine the success probability of the reconstruction procedure
   with small ell */

FILE *f; 

bool run_test()
{
	pcg128_t S[nbiter];

	if (fread(&S[0], sizeof(pcg128_t), 1, f) != 1) {
		perror("Something went wrong when reading /dev/urandom");
		exit(EXIT_FAILURE);
	} 	
        
	// prepare input
	u64 X[nbiter];
    	for (int i = 1; i < nbiter; i++)
    		S[i] = a * S[i-1] + c;
    	for (int i = 0; i < nbiter; i++)
    		X[i] = pcg_output_xsl_rr_128_64(S[i]);

    	// cheat
      	u64 W0 = (u64) (S[0] % (1 <<known_low));
	int rot[nbiter];
	for (int i = 0; i < nbiter; i++)
		rot[i] = (int) (S[i] >> 122);
	
	struct task_t task;
	setup_task(W0, X, &task);
	refresh_task(rot, &task);

	pcg128_t S_target[nbiter];
	if (solve(S_target, rot, &task)) {
	    	assert(S_target[0] == S[0]);
		return true;
	}
	return false;
}


int main()
{
    init_var_globales();
    f = fopen("/dev/urandom", "r");
    if (f == NULL) {
		perror("Something went wrong when opening /dev/urandom");
		exit(EXIT_FAILURE);
    }

    int good = 0;
    for (int foobar = 0; foobar < 1000000; foobar++) {
    	if (run_test())
    		good++;
    }
    printf("Total good = %d\n", good);
    return EXIT_SUCCESS;
}