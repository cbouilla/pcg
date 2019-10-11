#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    
    unsigned long long X[nbiter] = {9067441263659890769llu, 12674224149338242009llu, 2612107274013664962llu};
    
    unsigned long long W0;
    int rot[nbiter];
    unsigned long long sumPol[nbiter];
    unsigned long long sumPolY[nbiter];
    mpz_t S[nbiter];
    mpz_t polW[nbiter];
    for(int i = 0 ; i < nbiter ; i++){
        mpz_init(S[i]);
        mpz_init(polW[i]);
    }
    
    float temps;
    clock_t t1, t2;
    t1 = clock();
    for(W0 = 0 ; W0 < (1<<known_low) ; W0++){//(1<<known_low)
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
            //decaler X à rajouter avec vraie sortie de pcg
            if(solve(S, X, rot,sumPol,sumPolY)){
                printf("S :\n");
                for(int i = 0 ; i < nbiter ; i++)
                    gmp_printf("Si :%Zd, rot :%d %d %d, W0 : %llu\n",S[i],rot[0], rot[1], rot[2],W0);
            }
        }
    }
    t2 = clock();
    temps = (float)(t2-t1)/CLOCKS_PER_SEC;
    printf("temps = %f\n", temps);
    return(0);
}