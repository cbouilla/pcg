#include <stdint.h>
#include <stdbool.h>

#include "pcg_oneseq.h"

#include <omp.h>


/***** Macro et Variables globales *****/
#define k 64
#define known_up 6
#define known_low 18
#define nbiter 3

int nb_thread;
pcg128_t a;
pcg128_t a_inv;
pcg128_t c;
pcg128_t polC[nbiter];


typedef unsigned long long u64;

struct task_t {
        u64 sumPol[nbiter];
        u64 sumPolY[nbiter];
        u64 X[nbiter][64];
        double Yprime[nbiter][64];
};


/***** Fonctions *****/
void init_var_globales();
double wtime();

////////////////Fonctions pour la récupération de S//////////////
void rotate(unsigned long long* rX, const unsigned long long* X, const int* rot);
void unrotate(u64* urX, const u64* X, const int* rot);
void getPolW(pcg128_t *polW, u64 W0);
void getSumPol(u64* sumPol, u64* sumPolY, const pcg128_t* polW);

void setup_task(u64 W0, const u64 *X, struct task_t *task);


bool solve(pcg128_t* S, const int* rot, const struct task_t *task);
