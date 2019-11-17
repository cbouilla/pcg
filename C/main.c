#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main()
{
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** INPUT ***********/
    unsigned long long X[nbiter];
    X[ 0] = 0x3c8673c55ead80bc;
    X[ 1] = 0xfd27cfa6f84ef5f1;
    X[ 2] = 0x68a098da45558d92;

    double t1 = wtime();
    unsigned long long done = 0;
    
    #pragma omp parallel for
    for (unsigned long long W0 = 0; W0 < (1<<known_low) ; W0++){//(1<<known_low)//W0=209818  
        
        /*Variables privées*/
        pcg128_t S[nbiter];
        pcg128_t polW[nbiter];
        unsigned long long urX[nbiter];
        int rot[nbiter];
        unsigned long long sumPol[nbiter];
        unsigned long long sumPolY[nbiter];
        
        getPolW(polW, W0);
        getSumPol(sumPol, sumPolY, polW);
        
        if (omp_get_thread_num() == 0 && (W0 % 64) == 0) {
            printf("\rW0 = %llx / %llx --- %.1f / s", done, 1ull<<known_low, done / (wtime() - t1));
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
            unrotate(urX, X, rot);
            
            /***** Résolution *****/
            if (solve(S, urX, rot, sumPol, sumPolY)) {
                printf("\nInternal state found (%.1fs)\n", wtime() - t1);
                struct pcg_state_128 rng;
                rng.state = S[0];
                for(int i = 1 ; i < 10 ; i++)
                    printf("X[%2d] = %016" PRIx64 "\n", i, pcg_oneseq_128_xsl_rr_64_random_r(&rng));
            }
        }
        
        #pragma omp atomic
        done++;;
    }
    printf("Total time = %.1f\n", wtime() - t1);
    return EXIT_SUCCESS;
}