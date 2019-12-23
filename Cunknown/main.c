#include <stdlib.h>
#include <stdio.h>

#include "fonctions.h"

int main()
{
    /*  INITIALISATION DES PARAMETRES  */
    init_var_globales();
        
    /********** prepare test input ***********/
    u64 X[nboutput];
    pcg128_t S0 = (((pcg128_t) 5995207026785010249u) << k) + ((pcg128_t) 179350442155841024u);
    pcg128_t c = ((((pcg128_t) 6364136223846793005u) << k) + 1442695040888963407u) >> 1;
    pcg128_t vraiS[nboutput];
    pcg(vraiS, X, S0, &c, nboutput);

    double t1 = wtime();
    u64 W0 = 5018, WC = 335;

    /* either one task per thread, or we paralellize the "r" loop below */
    struct task_t task;
    init_task(&task);

    for (W0 = 5018; W0 < /*(1<<known_low)*/ 5019 ; W0++){//W0=5018 
        for (WC = 335 ; WC < /*(1<<known_low)*/ 336 ; WC++){//WC = 335
            
            prepare_task(X, W0, WC, &task);

            for(int r = 0 ; r < 1<<(nbiter * known_up) ; r++) {
                
                /***** Modification de rot et unrotX *****/
                task.rot[0] = (task.rot[0] + 1) % k;
                int i = 0;
                while (task.rot[i] == 0 && i < nbiter){
                    i++;
                    task.rot[i] = (task.rot[i] + 1) %k;
                }
                
                if (solve_isgood(task.goodY, task.rot, task.tabTmp, task.sumPolY, task.sumPolTest)) {
                    u64 DS640;
                    u64 Y0;
                	solve(&DS640, &Y0, task.goodY, task.rot, task.tabTmp, task.sumPolY, task.sumPolTest);
                    printf("candidat DS64 trouvé !!\n");
                    printf("%llu\n", DS640);
                    printf("temps pour trouver la solution = %f\n", wtime() - t1 );
                }
                /*if(DS640 == 7304601715607344736u){
                    printf("On a le bon !\n");
                    printf("DS640 = %llu\n", DS640);
                }*/
            }
            getGoodY(task.goodY, X, task.lowSumPol, 0);
        }
    }
    printf("temps total = %f\n", wtime() - t1);
    return(0);
}
