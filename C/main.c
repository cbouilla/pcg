#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    pcg128_t vraiS[nbiter];
    unsigned long long X[nbiter]; //= {15511475948252377726llu, 6234043309952295463llu, 9845631557281214488llu};
    pcg128_t S0 = (((pcg128_t) 5995207026785010249) << k) + ((pcg128_t) 179350442155841024); //valeur aléatoire
    pcg(vraiS, X, S0, nbiter);
    printf("W0 : %llu, rot = %d, %d, %d \n", (unsigned long long) (vraiS[0]% (1<<known_low)), (int) (vraiS[0] >> (2*k - known_up)),(int) (vraiS[1] >> (2*k - known_up)),(int) (vraiS[2] >> (2*k - known_up)));
    unsigned long long W0;
    
    FILE *f;
    f = fopen("result.txt","w");
    
    float temps;
    clock_t t1, t2;
    t1 = clock();
    #pragma omp parallel for
    for(W0 = 0 ; W0 < (1<<known_low) ; W0++){//(1<<known_low)//W0=209818  
        
        /*Variables privées*/
        pcg128_t S[nbiter];
        pcg128_t polW[nbiter];
        unsigned long long urX[nbiter];
        int rot[nbiter];
        unsigned long long sumPol[nbiter];
        unsigned long long sumPolY[nbiter];
        
        if((W0 % (1<<10)) == 0){
            printf("num de th :%d\n",omp_get_thread_num());
            printf("W0 : %llu\n", W0);
            fprintf(f, "W0 : %llu\n", W0);
            fflush(f);
        }
        getPolW(polW, W0);
        getSumPol(sumPol,sumPolY, polW);
        for(int r = 0 ; r < 1<<(3*known_up) ; r++){
            //modif de rot :
            rot[0]=(rot[0] + 1) % k;
            int i = 0;
            while(rot[i] == 0 && i < nbiter){
                i++;
                rot[i]=(rot[i] + 1) %k;
            }
            unrotate(urX, X, rot);
            
            if(solve(S, urX, rot,sumPol,sumPolY)){
                fprintf(f,"S :\n");
                printf("S trouvé !!\n");
                for(int i = 0 ; i < nbiter ; i++)
                    fprintf(f,"%016llx %016llx\n", (unsigned long long) (S[i]>>64), (unsigned long long) S[i]);
                t2 = clock();
                temps = (float)(t2-t1)/CLOCKS_PER_SEC;
                fprintf(f,"temps pour trouver la solution = %f\n", temps);
                fflush(f);
            }
        }
    }
    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    fprintf(f,"temps total = %f\n", temps);
    return(0);
}