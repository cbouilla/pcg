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
    a = (((pcg128_t) 2549297995355413924) << k) + ((pcg128_t) 4865540595714422341);
    
    a_inv = (((pcg128_t) 566787436162029664) << k) + ((pcg128_t) 11001107174925446285ull);

    //increment c OK !
    c = (((pcg128_t)6364136223846793005) << k) + 1442695040888963407;
    
    //nombre de threads 
    nb_thread = omp_get_max_threads();
    
    //increment polynome polC OK !
    polC[0] = 0;
    pcg128_t powA =1;
    for (int i = 1; i < nbiter ; i++){
        polC[i] = polC[i-1] + powA * c;
        powA *= a;
    }

    for (int i = 0; i < 9; i++)
        invG[i] *= 1ll << (k - known_low - known_up);
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

void rotate(unsigned long long* rX, const unsigned long long* X, const int* rot)
{ //pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        rX[i]= (X[i] >> rot[i]) | (X[i] << ((- rot[i]) & 63));
}

void unrotate(u64* urX, const u64* X, const int* rot)
{//pas verifié, repris de pcg_random
    int rot2[nbiter];
    for( int i = 0 ; i < nbiter ; i++)
        rot2[i] = (k - rot[i]) % k;
    rotate(urX, X, rot2);
}

void getPolW(pcg128_t *polW, u64 W0)
{ //OK !
    polW[0] = W0;
    for(int i = 1 ; i < nbiter ; i++){
        polW[i] = polW[i-1] * a;
    }
}

void getSumPol(u64* sumPol, u64* sumPolY, const pcg128_t* polW)
{
    pcg128_t sum;
    for(int i = 0 ; i < nbiter ; i++){
        sum = polC[i] + polW[i];
        sumPol[i] = sum;
        sumPolY[i] = (sum >> (k - known_up)) % (1 << (known_low + known_up));
    }
}

/* cf. https://stackoverflow.com/questions/17035464/a-fast-method-to-round-a-double-to-a-32-bit-int-explained#comment61972557_17035583 */
static inline long long crazy_round(double x)
{
    union { double d; long long l; } magic; 
    magic.d = x + 6755399441055744.0; 
    magic.l <<= 13; 
    magic.l >>= 13;
    return magic.l;
}

static inline long long light_crazy_round(double x)
{
    union { double d; long long l; } magic; 
    magic.d = x;
    magic.l <<= 13; 
    magic.l >>= 13;
    return magic.l;
}


void setup_task(u64 W0, const u64 *urX, struct task_t *task)
{
        pcg128_t polW[nbiter]; // temporary byproduct
        u64 sumPol[nbiter]; // temporary byproduct
        u64 sumPolY[nbiter]; // temporary byproduct
        
        // initialise polW à partir de W0
        getPolW(polW, W0);

        // initialise sumPol et sumPolY à partir de polW
        getSumPol(sumPol, sumPolY, polW);

        for (int i = 0; i < nbiter; i++)
                for (int j = 0; j < 64; j++) {
                        u64 X = rot(urX[i], j);
                        u64 Y = (((sumPol[i] ^ X) % (1 << known_low)) << 6) + (j ^ (X >> 58));
                        u64 Yprime = (Y - sumPolY[i]) % (1 << (known_low + 6));
                        
                        task->X[i][j] = X;
                        task->Yprime[i][j] = (double) Yprime; // conversion en double
                }
        task->sumPol0 = sumPol[0];
}


bool solve(pcg128_t* S, const int* rot, const struct task_t *task)
{  
    double tmp2[nbiter];
    // total : 9 MUL (double) + 10 ADD (double)
    for (int i = 0; i < nbiter; i++) {
        tmp2[i] = 0;
        for(int j = 0; j < nbiter; j++)
            tmp2[i] += invG[i * nbiter + j] * task->Yprime[j][rot[j]];
        tmp2[i] += 6755399441055744.0;
    }
 
    // total : 6 SHIFT
    u64 tmp3[nbiter];
    for(int i = 0 ; i < nbiter ; i++)
        tmp3[i] = light_crazy_round(tmp2[i]);

    u64 Sprim0 = 0;
    // total : 3 ADD, 3 MUL
    for(int j=0 ; j<nbiter ; j++)
        Sprim0 += Greduite[j] * tmp3[j];

    // SHIFT, ADD
    u64 Smod = (Sprim0 << known_low) + task->sumPol0;
    S[0] = (((pcg128_t)(Smod ^ task->X[0][rot[0]])) << 64) + ((pcg128_t) Smod);

    // clock PCG two times, check that we reconstruct the correct (de-rotated) output
    for (int i = 1; i < nbiter; i++) {
        S[i] = S[i - 1] * a + c;
        u64 XX = S[i] ^ (S[i] >> 64);
        if (XX != task->X[i][rot[i]])
            return false;
    }
    return true;
}
 