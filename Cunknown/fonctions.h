#include <stdint.h>
//#include "pcg_oneseq.h"
#include "pcg_setseq.h" //inclus dans pcg_oneseq
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 11
#define nbiter 5
#define nboutput 31
#define nbtest 3

int nb_thread;
pcg128_t a;
pcg128_t powA[nboutput];
pcg128_t polA[nboutput];

extern unsigned long long Greduite[16];
extern float invG[16];

/***** Fonctions *****/
void init_var_globales();
double wtime();

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n);

////////////////Fonctions pour la récupération de S//////////////
void rotateX(unsigned long long* rX, unsigned long long* X,int* rot);
void unrotateX(unsigned long long* urX, unsigned long long* X, int* rot);
unsigned long long unrotate1(unsigned long long Xi);


void getY(unsigned long long *Y, unsigned long long W0, unsigned long long WC, int* rot, unsigned long long* uX);
void getYprim(unsigned long long *Yprim, unsigned long long *Y, unsigned long long W0, unsigned long long WC);
void getDY(unsigned long long *DY, unsigned long long* Yprim);
void FindDS64(unsigned long long* DS64, unsigned long long *Y, unsigned long long* uX,int* rot, unsigned long long* lowSumPol, unsigned long long* sumPolY);

unsigned long long FindDS640(unsigned long long* Y, unsigned long long* uX, int* rot,unsigned long long *lowSumPol,unsigned long long* sumPolY);

int testDS640(unsigned long long DS640,  unsigned long long* X, unsigned long long Y0, unsigned long long* sumPolTest, unsigned long long* lowSumPol);

int testValid(FILE* f, int n);
//void pcgone(pcg128_t *S, unsigned long long* X, pcg128_t S0, int n);
void pcg(pcg128_t *S, unsigned long long* X, pcg128_t S0, pcg128_t* c, int n);

/***** Tests *****/
int testFonctions();
void printVal(pcg128_t S0, pcg128_t c);
