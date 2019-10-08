#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    int k = 64;
    int known_up = 6;
    int known_low = 18;
    int nbiter = 3;
    
    //m le modulo 2**128
    mpz_t m;
    mpz_init_set_ui(m,1);
    mpz_mul_2exp(m,m,k*2);
    
    //n le modulo 2**64
    
    //multiplier a OK !
    mpz_t a;
    mpz_init_set_ui(a,2549297995355413924); 
    mpz_mul_2exp(a,a,k);
    mpz_add_ui(a,a,4865540595714422341);
    
    //increment c OK !
    mpz_t c;
    mpz_init_set_ui(c,6364136223846793005); 
    mpz_mul_2exp(c,c,k);
    mpz_add_ui(c,c,1442695040888963407);
    
    //increment polynome polC OK !
    mpz_t* polC = malloc(nbiter*sizeof(mpz_t));
    mpz_init_set_ui(polC[0], 0);
    mpz_t powA;
    mpz_init_set_ui(powA,1);
    for(int i = 1; i < nbiter ; i++){
        mpz_init_set(polC[i],polC[i-1]);
        mpz_addmul(polC[i],powA,c);
        mpz_mod(polC[i],polC[i],m);
        mpz_mul(powA, powA, a);
        mpz_mod (powA, powA, m);
    }
        
    //Matrice génératrice réduite grâce à Sage à rentrer (et son inverse)
    //j'ai pris -2**63 <= long long int <2**63, faut voir si ça fait des bugs, peut etre ce sera moins efficace que si on utilise unsigned long long avec les erreurs de retenue...
    
    unsigned long long Greduite[9] = {-1241281756092, -5001120657083,  8655886039732,
     3827459685972, -2117155768935,  3303731088004,
     -728312298332,  5479732607037,  6319848582548};

    float invG[9] = {-9.25221813226351e-14,  2.32272588749499e-13,  5.30001389997814e-15,
    -7.81560146462246e-14, -4.52719459143047e-15,  1.09411828555506e-13,
     5.71040610630735e-14,  3.06929503514353e-14,  6.39750297008745e-14};    
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    unsigned long long W0 = 234270;
    unsigned long long rot[3] = {25, 53, 44};
    unsigned long long X[3] = {9681998970132172112, 3255734256114478613, 12077861056297022149};
    mpz_t* polW = malloc(nbiter*sizeof(mpz_t));
    getPolW(polW, W0, a, m, nbiter);
    /*for(int i = 0 ; i < nbiter ; i++)
        gmp_printf("%Zd\n",polW[i]);*/
    unsigned long long* Y = malloc(nbiter*sizeof(unsigned long long));
    getY(Y, polW, polC, rot, X, nbiter, known_low, known_up, k);
    printf("Y :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Y[i]);
    unsigned long long* Yprim = malloc(nbiter*sizeof(unsigned long long));
    getYprim(Yprim, Y, polW, polC, nbiter, known_low, known_up, k);
    printf("Yprim :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Yprim[i]);
    
    unsigned long long* Sprim = malloc(nbiter*sizeof(unsigned long long));
    findSprim(Sprim, Yprim, Greduite, invG, known_low, known_up, k, nbiter);
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu\n",Sprim[i]);
    
    mpz_t* S = malloc(nbiter*sizeof(mpz_t));
    findS(S, Sprim, X, polC, polW, known_low, k, nbiter);
    for(int i = 0 ; i < nbiter ; i++)
        printf("%Zd\n",S[i]);
    
    //memory liberation
    free(polC);
    free(polW);
    free(Y);
    free(Yprim);
    free(S);
    
    return(0);
}