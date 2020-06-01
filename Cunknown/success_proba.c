#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "fonctions.h"

FILE *f; 
struct task_t task;


bool run_test()
{
 	
	pcg128_t S[nboutput];
	pcg128_t c;

	if (fread(&S[0], sizeof(pcg128_t), 1, f) != 1) {
		perror("Something went wrong when reading /dev/urandom");
		exit(EXIT_FAILURE);
	} 	
        if (fread(&c, sizeof(pcg128_t), 1, f) != 1) {
		perror("Something went wrong when reading /dev/urandom");
		exit(EXIT_FAILURE);
	}
	c |= 1;

	// prepare input
	u64 X[nboutput];
    	for (int i = 1; i < nboutput; i++)
    		S[i] = a * S[i-1] + c;
    	for (int i = 0; i < nboutput; i++)
    		X[i] = pcg_output_xsl_rr_128_64(S[i]);

    	// cheat
      	u64 W0 = (u64) (S[0] % (1 <<known_low));
      	u64 Wc = (u64) (c    % (1 <<known_low));
	   assert((Wc % 2) == 1);

	prepare_task(&X, W0, Wc, &task);
	 
    for (int i = 0; i < nbiter; i++)
        task.rot[i] = (int) (S[i] >> 122);

	 /* check */
	bool ok = false;
        if (solve_isgood(&task)) {
            u64 DS640;
            u64 Y0;
            solve(&task, &DS640, &Y0);
            ok = true;
        }
        
        /* reset goodY */
    finish_task(&X, &task);
	return ok;
}


int main()
{
    
    /*  INITIALISATION DES PARAMETRES  */    
    init_var_globales();

    f = fopen("/dev/urandom", "r");
    if (f == NULL) {
		perror("Something went wrong when opening /dev/urandom");
		exit(EXIT_FAILURE);
    }
    init_task(&task);

    int good = 0;
    for (int foobar = 0; foobar < 1000000; foobar++) {
    	if (run_test())
    		good++;
    }
    printf("Total good = %d\n", good);
    return EXIT_SUCCESS;
}