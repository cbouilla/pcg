#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    pcg128_t S0 = (((pcg128_t) 5995207026785010249u) << k) + ((pcg128_t) 179350442155841024u);
    pcg128_t c = ((((pcg128_t) 6364136223846793005u) << k) + 1442695040888963407u) >> 1;
    
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    
   // printVal(S0, c);
    
    pcg(vraiS, X, S0, &c, nboutput);
    double t1 = wtime();
    unsigned long long done = 0;
    unsigned long long W0 = 5018, WC = 335;
    
    //#pragma omp parallel for
    for (W0 = 5018; W0 < /*(1<<known_low)*/ 5019 ; W0++){//(1<<known_low)//W0=5018 
        for(WC = 335 ; WC < /*(1<<known_low)*/ 336 ; W0++){//WC = 335
            /*Variables privées*/
            unsigned long long uX[nbiter];
            int rot[nbiter];
            for(int i = 0 ; i < nbiter ; i++)
                    rot[i] = 0;
            unsigned long long DS64[nboutput];
            unsigned long long Y[nbiter];
            
            /*if (omp_get_thread_num() == 0) {
                printf("\rW0 = %llx / %llx --- %.1f / s", done, 1ull<<known_low, done / (wtime() - t1));
                fflush(stdout);
            }*/
            for(int r = 0 ; r < 1<<(nbiter * known_up) ; r++){//cette boucle n'est pas parralélisable direct !
                /*unsigned long long uX[nbiter];
                int rot[nbiter];
                for(int i = 0 ; i < nbiter ; i++)
                    rot[i] = 0;
                unsigned long long DS64[nboutput];
                unsigned long long Y[nbiter];*/
                
                /***** Modification de rot et unrotX *****/
                rot[0]=(rot[0] + 1) % k;
                int i = 0;
                while(rot[i] == 0 && i < nbiter){
                    i++;
                    rot[i]=(rot[i] + 1) %k;
                }
                
                unrotateX(uX, X, rot);
                
                /***** Résolution *****/
                getY(Y, W0, WC, rot, uX);
                
                FindDS64(DS64, uX,  rot, W0, WC);
                if(DS64[0] == 7304601715607344736u)
                    printf("On a le bon !\n");
                if(testDS640(DS64[0], X, Y[0], W0, WC, 3)){
                    printf("candidat DS64 trouvé !!\n");
                    printf("%llu\n", DS64[0]);
                    printf("temps pour trouver la solution = %f\n", wtime() - t1 );
                }
            }
            
            //#pragma omp atomic
           // done++;;
        }
    }
    printf("temps total = %f\n", wtime() - t1);
    return(0);
}