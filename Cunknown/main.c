#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "fonctions.h"


void do_task(u64 current, struct task_t *task, const u64 *X)
{
    u64 W0 = current >> known_low;
    u64 WC = current % (1 << known_low);
    printf("Doing task %lld / %lld\n", W0, WC);
    
    prepare_task(X, W0, WC, task);
    for (u64 r = 0; r < 1 << (nbiter * known_up); r++) {
        task->rot[0] = (task->rot[0] + 1) % k;
        int i = 0;
        while (task->rot[i] == 0 && i < nbiter) {
            i++;
            task->rot[i] = (task->rot[i] + 1) % k;
        }
            
        if (solve_isgood(task)) {
            u64 DS640;
            u64 Y0;
            solve(task, &DS640, &Y0);

            printf("thread %d, candidat DS64 trouvé !! %lld\n", omp_get_thread_num(), DS640);
        }
    }
    finish_task(X, task);
}


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
    
    // restart();

    u64 range_start = (5018 << known_low) + 335;
    u64 range_end = range_start + 10;

    double t1 = wtime();

    /* init all tasks */
    int T = omp_get_max_threads();
    struct task_t task[T];
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        init_task(&task[tid]);
    }

    u64 current = range_start;
    while (current < range_end) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            if (current + tid < range_end) {
                do_task(current + tid, &task[tid], X);
            }
        }
        // checkpoint();
        current += T;
    }


    /*if(DS640 == 7304601715607344736u){
        printf("On a le bon !\n");
        printf("DS640 = %llu\n", DS640);
    }*/

    printf("temps total = %f\n", wtime() - t1);
    return(0);
}
