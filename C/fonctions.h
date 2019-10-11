#include <stdint.h>
#include <gmp.h>

/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 18
#define nbiter 3

mpz_t m;
mpz_t a;
mpz_t c;
mpz_t polC[nbiter];

extern unsigned long long Greduite[9];

extern float invG[9];

/***** Fonctions *****/
void init_var_globales();

void prodMatVecUIU(unsigned long long* res, unsigned long long* M, unsigned long long* v, int n);

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n);

void prodMatMatI(unsigned long long* res, unsigned long long* M1, unsigned long long* M2, int n);


////////////////Fonctions pour la récupération de S//////////////
void rotate(unsigned long long* X,int* rot);
void unrotate(unsigned long long* X,int* rot);

void getPolW(mpz_t *polW, unsigned long long W0);

void getSumPol(unsigned long long* sumPol,unsigned long long* sumPolY, mpz_t* polW);

void getY(unsigned long long* Y, unsigned long long* sumPol, int* rot, unsigned long long* X);

void getYprim(unsigned long long* Yprim, unsigned long long* Y, unsigned long long* sumPolY);

void findSprim(unsigned long long* Sprim, unsigned long long* Yprim);

void findS(mpz_t* S, unsigned long long* Sprim, unsigned long long* X, unsigned long long*sumPol);

int test(mpz_t* S, unsigned long long* X);

int solve(mpz_t* S, unsigned long long* X, int* rot,unsigned long long* sumPol,unsigned long long* sumPolY);

