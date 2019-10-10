#include <stdint.h>
#include <gmp.h>

/***** Variables globales *****/
#define k 64
#define known_up 6
#define known_low 18
#define nbiter 3

void prodMatVecUIU(unsigned long long* res, unsigned long long* M, unsigned long long* v, int n);

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n);

void prodMatMatI(unsigned long long* res, unsigned long long* M1, unsigned long long* M2, int n);


////////////////Fonctions pour la récupération de S//////////////

void getPolW(mpz_t *polW, unsigned long long W0, mpz_t a, mpz_t m);

void getSumPol(unsigned long long* sumPol,unsigned long long* sumPolY, mpz_t* polC, mpz_t* polW);

void getY(unsigned long long* Y, unsigned long long* sumPol, unsigned long long* rot, unsigned long long* X);

void getYprim(unsigned long long* Yprim, unsigned long long* Y, unsigned long long* sumPolY);

void findSprim(unsigned long long* Sprim, unsigned long long* Yprim, unsigned long long * Greduite, float* invG);

void findS(mpz_t* S, unsigned long long* Sprim, unsigned long long* X, unsigned long long*sumPol);

int test(mpz_t* S, unsigned long long* X);

