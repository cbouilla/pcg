#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "fonctions.h"

#define WORK_FACTOR (1ull << 27)

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
    
    u64 W0 = 5018, WC = 335;

    /**** Polynômes en WC et W0 utilisés dans la résolution ****/
    struct task_t task;
    init_task(&task);
    prepare_task(X, W0, WC, &task);

    printf("known_low = %d\n", known_low);
    printf("Taille de GoodY = %d Ko\n", (1 << (known_low + known_up)) / 1024 / 8);
    printf("Début du benchmark (%llu iterations)\n", WORK_FACTOR);
    double t1 = wtime();

    for (u64 r = 0 ;r < WORK_FACTOR; r++) {
        
        /***** Modification de rot et unrotX *****/
        task.rot[0] = (task.rot[0] + 1) % k;
        int i = 0;
        while (task.rot[i] == 0 && i < nbiter) {
            i++;
            task.rot[i] = (task.rot[i] + 1) % k;
        }
        
        if (solve_isgood(&task)) {
            printf("candidat DS64 trouvé !!\n");
            printf("temps pour trouver la solution = %f\n", wtime() - t1);
        }
    }

    double t = wtime() - t1;
    printf("Durée benchmark = %.2fs\n", t);
    printf("Itérations/s = %.1fM/s\n", WORK_FACTOR / t / 1e6);
    printf("Concrètement = %d tasks of size %.1f h-CPU\n", 1 << known_low, 
        t / WORK_FACTOR * (1ull << (nbiter * known_up + known_low - 1)) / 3600);
    printf("Attaque complète = %.0fK h-CPU\n", t / WORK_FACTOR * (1ull << (nbiter * known_up + 2*known_low - 1)) / 3600 / 1e3);
    
    exit(0);
}
