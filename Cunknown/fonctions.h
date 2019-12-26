#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "pcg_setseq.h"

typedef long long          i64;
typedef unsigned long long u64;

/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 13
#define nbiter 5
#define nbtest 4
#define nboutput (nbiter + nbtest)

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
void finish_task(const u64 *X, struct task_t *t);
double wtime();
bool solve_isgood(const struct task_t *task);
void solve(const struct task_t *task, u64* DS640, u64* Y0);
void pcg(pcg128_t *S, u64* X, pcg128_t S0, pcg128_t* c, int n);

///// tout ceci ne sert que dans test.c

void getY(u64 *Y, u64 W0, u64 WC, int* rot, u64* uX);
void getYprim(u64 *Yprim, u64 *Y, u64 W0, u64 WC);
void getDY(u64 *DY, u64* Yprim);
void FindDS64(u64* DS64, u64 *Y0, u64* uX,int* rot, u64* lowSumPol, u64* sumPolY);
u64 FindDS640(u64* Y, u64* uX, int* rot,u64 *lowSumPol,u64* sumPolY);
int testDS640(u64 DS640,  u64* X, u64 Y0, u64* sumPolTest, u64* lowSumPol);