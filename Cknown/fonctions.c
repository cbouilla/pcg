#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "fonctions.h"


unsigned long long Greduite[9] = {
    -1241281756092, -5001120657083,  8655886039732,
     3827459685972, -2117155768935,  3303731088004,
     -728312298332,  5479732607037,  6319848582548
 };

double invG[9] = {
    -9.25221813226351e-14,  2.32272588749499e-13,  5.30001389997814e-15,
    -7.81560146462246e-14, -4.52719459143047e-15,  1.09411828555506e-13,
     5.71040610630735e-14,  3.06929503514353e-14,  6.39750297008745e-14
};

void init_var_globales(){
    //multiplier a OK !
    a = (((pcg128_t) 2549297995355413924) << 64) + ((pcg128_t) 4865540595714422341);
    
    a_inv = (((pcg128_t) 566787436162029664) << 64) + ((pcg128_t) 11001107174925446285ull);

    //increment c OK !
    c = (((pcg128_t)6364136223846793005) << 64) + 1442695040888963407;
    
    //nombre de threads 
    nb_thread = omp_get_max_threads();
    
    //increment polynome polC OK !
    polC[0] = 0;
    pcg128_t powA =1;
    for (int i = 1; i < 3 ; i++){
        polC[i] = polC[i-1] + powA * c;
        powA *= a;
    }

    for (int i = 0; i < 9; i++)
        invG[i] *= 1ll << (64 - known_low - 6);
}


//////////////////// chrono //////////////////
double wtime()
{
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}


////////////////Fonctions pour la récupération de S//////////////

static u64 rot(u64 x, int i)
{
        return (x >> (64 - i)) | (x << i);
}


void getPolW(pcg128_t *polW, u64 W0)
{ //OK !
    polW[0] = W0;
    for(int i = 1 ; i < 3 ; i++){
        polW[i] = polW[i-1] * a;
    }
}

void getSumPol(u64* sumPol, u64* sumPolY, const pcg128_t* polW)
{
    pcg128_t sum;
    for(int i = 0 ; i < 3 ; i++){
        sum = polC[i] + polW[i];
        sumPol[i] = sum;
        sumPolY[i] = (sum >> 58) % (1 << (known_low + 6));
    }
}

/* cf. https://stackoverflow.com/questions/17035464/a-fast-method-to-round-a-double-to-a-32-bit-int-explained#comment61972557_17035583 */
static inline long long light_crazy_round(double x)
{
    union { double d; long long l; } magic; 
    magic.d = x;
    magic.l <<= 13; 
    magic.l >>= 13;
    return magic.l;
}

static inline long long crazy_round(double x)
{
    return light_crazy_round(x + 6755399441055744.0);
}


void setup_task(u64 W0, const u64 *urX, struct task_t *task)
{
        pcg128_t polW[3]; // temporary byproduct
        u64 sumPol[3]; // temporary byproduct
        u64 sumPolY[3]; // temporary byproduct
        
        // initialise polW à partir de W0
        getPolW(polW, W0);

        // initialise sumPol et sumPolY à partir de polW
        getSumPol(sumPol, sumPolY, polW);

        u64 Yprime[3][64];
        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 64; j++) {
                        u64 X = rot(urX[i], j);
                        u64 Y = (((sumPol[i] ^ X) % (1 << known_low)) << 6) + (j ^ (X >> 58));
                        Yprime[i][j] = (Y - sumPolY[i]) % (1 << (known_low + 6));
                        task->X[i][j] = X;
                }

        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 64; k++)
                                task->Ginv_Yprime[i][j][k] = invG[i * 3 + j] * Yprime[j][k];

        task->sumPol0 = sumPol[0];
}

void refresh_task(const int *rot, struct task_t *task)
{
    for (int i = 0; i < 3; i++) {
        task->tmp2[i] = 0;
        for(int j = 1; j < 3; j++)
            task->tmp2[i] += task->Ginv_Yprime[i][j][rot[j]];
    }
}

bool solve(pcg128_t* S, const int* rot, const struct task_t *task)
{  
    u64 Sprim0 = 0;
    for (int i = 0; i < 3; i++) {
        double tmp2 = task->tmp2[i] + task->Ginv_Yprime[i][0][rot[0]];
        Sprim0 += Greduite[i] * crazy_round(tmp2);
    }

    // SHIFT, ADD
    u64 Smod = (Sprim0 << known_low) + task->sumPol0;
    S[0] = (((pcg128_t) (Smod ^ task->X[0][rot[0]])) << 64) + ((pcg128_t) Smod);

    // clock PCG two times, check that we reconstruct the correct (de-rotated) output
    for (int i = 1; i < 3; i++) {
        S[i] = S[i - 1] * a + c;
        u64 XX = S[i] ^ (S[i] >> 64);
        if (XX != task->X[i][rot[i]])
            return false;
    }
    return true;
}
 