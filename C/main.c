#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    
    unsigned long long W0 = 98735;
    unsigned long long rot[nbiter] = {7, 54, 50};
    unsigned long long X[nbiter] = {9067441263659890769llu, 12674224149338242009llu, 2612107274013664962llu};
    mpz_t polW[nbiter];
    getPolW(polW, W0);
    /*for(int i = 0 ; i < nbiter ; i++)
        gmp_printf("%Zd\n",polW[i]);*/
    
    unsigned long long sumPol[nbiter];
    unsigned long long sumPolY[nbiter];
    getSumPol(sumPol,sumPolY, polW);
    
    unsigned long long Y[nbiter];
    getY(Y, sumPol, rot, X);
    printf("Y :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Y[i]);
    unsigned long long Yprim[nbiter];
    getYprim(Yprim, Y, sumPolY);
    printf("Yprim :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Yprim[i]);
    
    unsigned long long Sprim[nbiter];
    findSprim(Sprim, Yprim);
    printf("Sprim :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Sprim[i]);
    
    mpz_t S[nbiter];
    findS(S, Sprim, X, sumPol);
    printf("S :\n");
    for(int i = 0 ; i < nbiter ; i++)
        gmp_printf("%Zd\n",S[i]);
    
    printf("resultat du test %d\n",test(S,X));
    return(0);
}