#include <stdint.h>
#include <stdbool.h>
//#include "pcg_oneseq.h"
#include "pcg_setseq.h" //inclus dans pcg_oneseq
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef long long          i64;
typedef unsigned long long u64;

/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 11
#define nbiter 5
#define nboutput 31
#define nbtest 3

extern pcg128_t a;
extern pcg128_t powA[nboutput];
extern pcg128_t polA[nboutput];

extern unsigned long long Greduite[16] __attribute__((aligned(32)));
extern double invG[16] __attribute__((aligned(32)));

/***** Fonctions *****/
void init_var_globales();
double wtime();

static inline void prodMatVecFFU(double* res, double* M, unsigned long long* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

////////////////Fonctions pour la récupération de S//////////////
static inline void rotateX(unsigned long long* rX, const unsigned long long* X, const int* rot){ //pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        rX[i]= (X[i] >> rot[i]) | (X[i] << ((- rot[i]) & (k-1)));
}

static inline void unrotateX(unsigned long long* urX, const unsigned long long* X, const int* rot){//pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        urX[i]= (X[i] >> ((- rot[i]) & (k-1))) | (X[i] << rot[i]);
}

static inline unsigned long long unrotate1(unsigned long long Xi)
{
    return (Xi >> (k-1)) | (Xi << 1);
}

static inline unsigned long long unrotate(unsigned long long Xi, int i)
{
    return (Xi >> (k-i)) | (Xi << i);
}

char* setupGoodY();
void getGoodY(char* goodY, unsigned long long* tabX, unsigned long long* lowSumPol, int v);
void getTabTmp(unsigned long long* tabTmp, unsigned long long* X, unsigned long long* lowSumPol, unsigned long long* sumPolY);

void getY(unsigned long long *Y, unsigned long long W0, unsigned long long WC, int* rot, unsigned long long* uX);
void getYprim(unsigned long long *Yprim, unsigned long long *Y, unsigned long long W0, unsigned long long WC);
void getDY(unsigned long long *DY, unsigned long long* Yprim);
void FindDS64(unsigned long long* DS64, unsigned long long *Y0, unsigned long long* uX,int* rot, unsigned long long* lowSumPol, unsigned long long* sumPolY);

unsigned long long FindDS640(unsigned long long* Y, unsigned long long* uX, int* rot,unsigned long long *lowSumPol,unsigned long long* sumPolY);

int testDS640(unsigned long long DS640,  unsigned long long* X, unsigned long long Y0, unsigned long long* sumPolTest, unsigned long long* lowSumPol);

void solve(unsigned long long* DS640, unsigned long long* Y0, char* goodY, int* rot, unsigned long long* tabTmp, unsigned long long* sumPolY, unsigned long long* sumPolTest);

bool solve_isgood(const char* goodY, const int* rot, const unsigned long long* tabTmp, const unsigned long long* sumPolY, const unsigned long long* sumPolTest);

int testValid(FILE* f, int n);
//void pcgone(pcg128_t *S, unsigned long long* X, pcg128_t S0, int n);
void pcg(pcg128_t *S, unsigned long long* X, pcg128_t S0, pcg128_t* c, int n);

/***** Tests *****/
int testFonctions();
void printVal(pcg128_t S0, pcg128_t c);

