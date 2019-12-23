#include <assert.h>
#include <math.h>

#include "fonctions.h"

int testValid (FILE* f, int n) 
{
    int cpt = 0;
    struct task_t task;
    init_task(&task);

    for (int i = 0 ; i < n ; i++) {
        /* read random seed */
        pcg128_t seeds[2];
        if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
            perror("Something went wrong when reading /dev/urandom");
            exit(EXIT_FAILURE);
        }
        
        /* setup PRNG */
        pcg128_t vraiS[nboutput];
        u64 X[nboutput];
        pcg(vraiS, X, seeds[0], seeds+1, nboutput);
        
        /* extract "right" guessed values */
        u64 W0 = (u64) (vraiS[0] % (1 << known_low));
        u64 WC = (u64) (seeds[1] % (1 << known_low));
        prepare_task(X, W0, WC, &task);
        for(int i = 0; i < nbiter; i++)
            task.rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
        
        /**** Polynômes en WC et W0 utilisés dans la résolution ****/
        
        int a = solve_isgood(task.goodY, task.rot, task.tabTmp, task.sumPolY, task.sumPolTest); 
        cpt += a;
        if (a) {
            u64 DS640;
            u64 Y0;
            solve(&DS640, &Y0, task.goodY, task.rot, task.tabTmp, task.sumPolY, task.sumPolTest);
        }
        
        /* reset goodY */
        finish_task(X, &task);
    }
    return cpt;
}



int main()
{    
    /*  INITIALISATION DES PARAMETRES  */   
    init_var_globales();
    
    
    FILE *f = fopen("/dev/urandom", "r");
    if (f == NULL) {
        perror("Something went wrong when opening /dev/urandom");
        exit(EXIT_FAILURE);
    }
    
    printf("1..1\n");
    static const int nbtests = 1000000;
    int successes = testValid(f, nbtests);
    if (successes != nbtests)
        printf("not ok 1 - #success = %d / %d\n", successes, nbtests);
    else
        printf("ok 1 - all %d tests OK\n", nbtests);

    exit(0);
}
