#include <stdint.h>
#include <gmp.h>

void prodMatVecUIU(unsigned long long* res, unsigned long long* M, unsigned long long* v, int n);

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n);

void prodMatMatI(unsigned long long* res, unsigned long long* M1, unsigned long long* M2, int n);


////////////////Fonctions pour la récupération de S//////////////

void getPolW(mpz_t *polW, unsigned long long W0, mpz_t a, int k, int nbiter);

void getY(unsigned long long* Y, mpz_t* polW, mpz_t* polC, unsigned long long* rot, unsigned long long* X, int nbiter, int known_low, int known_up, int k);

void getYprim(unsigned long long* Yprim, unsigned long long* Y, mpz_t* polW, mpz_t* polC,
 int nbiter, int known_low, int known_up, int k);

void findSprim(unsigned long long* Sprim, unsigned long long* Yprim, unsigned long long * Greduite, float* invG, int known_low, int known_up, int k, int nbiter);

void findS(mpz_t* S, unsigned long long* Sprim, unsigned long long* X, mpz_t* polC, mpz_t* polW, int known_low, int k, int nbiter);

void Solve(mpz_t* S, unsigned long long* X, unsigned long long W0, unsigned long long* rot, unsigned long long* Greduite, float* invG, mpz_t* polC, mpz_t a, int known_low, int known_up, int k, int nbiter);