#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <math.h>
#include "fonctions.h"

////////////////Fonction conversion mpz - int64 //////////////
unsigned long long mpz_get_ull(mpz_t n)
{
    unsigned int lo, hi;
    mpz_t tmp;

    mpz_init( tmp );
    mpz_mod_2exp( tmp, n, 64 );   /* tmp = (lower 64 bits of n) */

    lo = mpz_get_ui( tmp );       /* lo = tmp & 0xffffffff */ 
    mpz_div_2exp( tmp, tmp, 32 ); /* tmp >>= 32 */
    hi = mpz_get_ui( tmp );       /* hi = tmp & 0xffffffff */

    mpz_clear( tmp );

    return (((unsigned long long)hi) << 32) + lo;
}

long long mpz_get_sll(mpz_t n) //a verif
{
    return (long long)mpz_get_ull(n); /* just use unsigned version */
}



////////////////Fonctions calcul matriciel//////////////

//Matrice avce coeffs négatifs mais résult mod 2**64 (je sais pas si ça marche)
void prodMatVecU(unsigned long long* res, unsigned long long * M, unsigned long long* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

void prodMatMatU(unsigned long long* res, unsigned long long* M1, unsigned long long* M2, int n){//juste pour verif
    int i, j, l;
    for(i=0 ; i<n ; i++){
        for(j=0 ; j<n ; j++){
            res[i * + j] = 0;
            for(l=0; l<n ; l++)
                res[i * n + j]+= M1[i * n + l] * M2[l * n + j];
        }
    }
}


////////////////Fonctions pour la récupération de S//////////////

void getPolW(mpz_t *polW, unsigned long long W0, mpz_t a, mpz_t m){ //OK !
    mpz_init_set_si(polW[0], W0);
    for(int i = 1 ; i < nbiter ; i++){
        mpz_init(polW[i]);
        mpz_addmul(polW[i], polW[i-1], a);
        mpz_mod(polW[i], polW[i], m);
    }
}

void getSumPol(unsigned long long* sumPol,unsigned long long* sumPolY, mpz_t* polC, mpz_t* polW){
    for(int i = 0 ; i < nbiter ; i++){
        mpz_t tmp;
        mpz_init(tmp);
        mpz_add(tmp, polC[i], polW[i]);
        sumPol[i] = mpz_get_ull(tmp);
        mpz_tdiv_q_2exp(tmp, tmp, k - known_up);
        sumPolY[i] = mpz_get_ull(tmp) % (1 << (known_low + known_up));
    }
}

void getY(unsigned long long* Y, unsigned long long* sumPol, unsigned long long* rot, unsigned long long* X){ //OK !
    for(int i = 0 ; i < nbiter ; i++){
        Y[i] = (((sumPol[i] ^ X[i]) % (1 << known_low)) << known_up ) + (rot[i] ^ (X[i] >> (k - known_up)));
    }
}

void getYprim(unsigned long long* Yprim, unsigned long long* Y, unsigned long long* sumPolY){//OK ! avec erreurs arrondi
     //Utiliser sumPolY
    mpz_t tmp;
    mpz_init(tmp);
    for(int i = 0 ; i < nbiter ; i++){
        Yprim[i] = Y[i] - sumPolY[i];
        Yprim[i] = Yprim[i] % (1 << (known_low + known_up));
    }
}

void findSprim(unsigned long long* Sprim, unsigned long long* Yprim, unsigned long long* Greduite, float* invG){ //OK !
    unsigned long long* tmp1 = malloc(nbiter * sizeof(unsigned long long));
    for(int i = 0 ; i < nbiter ; i++)
        tmp1[i] = Yprim[i] << (k - known_low - known_up);
    float* tmp2 = malloc(nbiter * sizeof(float));
    prodMatVecFFU(tmp2, invG, tmp1, nbiter);
    for(int i = 0 ; i < nbiter ; i++)
        tmp1[i] = (unsigned long long) roundf(tmp2[i]);
    prodMatVecU(Sprim, Greduite, tmp1, nbiter);
    free(tmp1);
    free(tmp2);
}
 
void findS(mpz_t* S, unsigned long long* Sprim, unsigned long long* X, unsigned long long* sumPol){
    unsigned long long tmp ;
    for(int i = 0 ; i < nbiter ; i++){
        tmp = (Sprim[i] << known_low) + sumPol[i];
        mpz_init_set_si(S[i], tmp ^ X[i]);
        mpz_mul_2exp(S[i], S[i], k);
        mpz_add_ui(S[i], S[i], tmp);
        mpz_mod_2exp(S[i], S[i], 2 * k);
    }
}

int test(mpz_t* S, unsigned long long* X){
    unsigned long long tmp;
    for(int i = 0 ; i < nbiter ; i++){
        tmp = mpz_get_ull(S[i]);
        mpz_tdiv_q_2exp(S[i], S[i], k);
        tmp = tmp ^ mpz_get_ull(S[i]);
        if(tmp != X[i])
            return 0;
    }
    return 1;
}
/* Sorties du générateur sans rotation 
void sortiesGenerateur(mpz_t* X, mpz_t* S, mpz_t m,  mpz_t a, mpz_t c, int nbiter){
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_init(S[0])
    mpz_urandomb (S[0], state, 2*k);
    for(int i = 1 ; i < nbiter ; i++){
        mpz_init_set(S[i],c)
        mpz_addmul(S[i], a, S[i-1])
        mpz_mod(
    
}*/
