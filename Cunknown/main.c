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
    
    //printVal(S0, c);
    
    pcg(vraiS, X, S0, &c, nboutput);
    double t1 = wtime();
    //unsigned long long done = 0;
    unsigned long long W0 = 5018, WC = 335;
    
    #pragma omp parallel for
    for (W0 = 5018; W0 < /*(1<<known_low)*/ 5019 ; W0++){//W0=5018 
        for(WC = 335 ; WC < /*(1<<known_low)*/ 336 ; WC++){//WC = 335

            /**** Polynômes en WC et W0 utilisés dans la résolution ****/
            unsigned long long lowSumPol[nbiter + nbtest];
            unsigned long long sumPolY[nbiter];
            unsigned long long sumPolTest[nbtest];
            for(int i = 0 ; i < nbiter ; i++){
                lowSumPol[i]  = (W0 * ((unsigned long long) powA[i]) + WC * ((unsigned long long) polA[i]));
                sumPolY[i] = (polA[i] * WC + powA[i] * W0) >> (k - known_up);
            }
            for(int i = 0 ; i < nbtest ; i++){
                lowSumPol[nbiter + i] = (W0 * ((unsigned long long) powA[i + nbiter]) + WC * ((unsigned long long) polA[i + nbiter]));
                sumPolTest[i] = W0 * ((unsigned long long) (powA[i + nbiter] >> known_low) - 1) + WC * ((unsigned long long) (polA[i + nbiter] >> known_low) - 1);
            }

            /*Variables privées*/
            int rot[nbiter];
            for(int i = 0 ; i < nbiter ; i++)
                    rot[i] = 0;
            unsigned long long DS640;
            unsigned long long Y0;


            for(int r = 0 ; r < 1<<(nbiter * known_up) ; r++){//cette boucle n'est pas parralélisable direct !
                
                /***** Modification de rot et unrotX *****/
                rot[0]=(rot[0] + 1) % k;
                int i = 0;
                while(rot[i] == 0 && i < nbiter){
                    i++;
                    rot[i]=(rot[i] + 1) %k;
                }
                
                if(solve(&DS640, &Y0, X, rot, lowSumPol, sumPolY, sumPolTest)){
                    printf("candidat DS64 trouvé !!\n");
                    printf("%llu\n", DS640);
                    printf("temps pour trouver la solution = %f\n", wtime() - t1 );
                }
                /*if(DS640 == 7304601715607344736u){
                    printf("On a le bon !\n");
                    printf("DS640 = %llu\n", DS640);
                }*/
            }
            
            //#pragma omp atomic
           // done++;;
        }
    }
    printf("temps total = %f\n", wtime() - t1);
    return(0);
}