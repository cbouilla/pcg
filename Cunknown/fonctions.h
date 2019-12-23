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

extern u64 Greduite[16];
extern double invG[16];

struct task_t {
    u64 lowSumPol[nbiter + nbtest];
    u64 sumPolY[nbiter];
    u64 sumPolTest[nbtest];
    u64 tabTmp[k * nbiter];
    int rot[nbiter];
    char *goodY;
};

/***** Fonctions *****/
void init_var_globales();
void init_task(struct task_t *t);
void prepare_task(const u64 *X, u64 W0, u64 WC, struct task_t *t);
double wtime();

static inline void prodMatVecFFU(double* res, double* M, u64* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

////////////////Fonctions pour la récupération de S//////////////
static inline void rotateX(u64* rX, const u64* X, const int* rot){ //pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        rX[i]= (X[i] >> rot[i]) | (X[i] << ((- rot[i]) & (k-1)));
}

static inline void unrotateX(u64* urX, const u64* X, const int* rot){//pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        urX[i]= (X[i] >> ((- rot[i]) & (k-1))) | (X[i] << rot[i]);
}

static inline u64 unrotate1(u64 Xi)
{
    return (Xi >> (k-1)) | (Xi << 1);
}

static inline u64 unrotate(u64 Xi, int i)
{
    return (Xi >> (k-i)) | (Xi << i);
}

char* setupGoodY();
void getGoodY(char* goodY, const u64* X, const u64* lowSumPol, int v);
void getTabTmp(u64* tabTmp, const u64* X, const u64* lowSumPol, const u64* sumPolY);

void getY(u64 *Y, u64 W0, u64 WC, int* rot, u64* uX);
void getYprim(u64 *Yprim, u64 *Y, u64 W0, u64 WC);
void getDY(u64 *DY, u64* Yprim);
void FindDS64(u64* DS64, u64 *Y0, u64* uX,int* rot, u64* lowSumPol, u64* sumPolY);

u64 FindDS640(u64* Y, u64* uX, int* rot,u64 *lowSumPol,u64* sumPolY);

int testDS640(u64 DS640,  u64* X, u64 Y0, u64* sumPolTest, u64* lowSumPol);

void solve(u64* DS640, u64* Y0, char* goodY, int* rot, u64* tabTmp, u64* sumPolY, u64* sumPolTest);

bool solve_isgood(const char* goodY, const int* rot, const u64* tabTmp, const u64* sumPolY, const u64* sumPolTest);

int testValid(FILE* f, int n);
//void pcgone(pcg128_t *S, u64* X, pcg128_t S0, int n);
void pcg(pcg128_t *S, u64* X, pcg128_t S0, pcg128_t* c, int n);

/***** Tests *****/
int testFonctions();
void printVal(pcg128_t S0, pcg128_t c);

