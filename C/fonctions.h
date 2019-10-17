#include <stdint.h>
#include "pcg_oneseq.h"
#include <omp.h>


/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 18
#define nbiter 3

int nb_thread;
pcg128_t a;
pcg128_t c;
pcg128_t polC[nbiter];

extern unsigned long long Greduite[9];

extern float invG[9];

/***** Fonctions *****/
void init_var_globales();
double wtime();


////////////////Fonctions pour la récupération de S//////////////
void rotate(unsigned long long* rX, unsigned long long* X,int* rot);
void unrotate(unsigned long long* urX, unsigned long long* X, int* rot);
void getPolW(pcg128_t *polW, unsigned long long W0);

void getSumPol(unsigned long long* sumPol,unsigned long long* sumPolY, pcg128_t* polW);

int solve(pcg128_t* S, unsigned long long* X, const int* rot, const unsigned long long* sumPol, const unsigned long long* sumPolY);

void pcg(pcg128_t *S, unsigned long long* X, pcg128_t S0, int n);
