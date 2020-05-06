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


void rotate(unsigned long long* rX, unsigned long long* X,int* rot){ //pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        rX[i]= (X[i] >> rot[i]) | (X[i] << ((- rot[i]) & 63));
}

void unrotate(unsigned long long* urX, unsigned long long* X, int* rot){//pas verifié, repris de pcg_random
    int rot2[nbiter];
    for( int i = 0 ; i < nbiter ; i++)
        rot2[i] = (k - rot[i]) % k;
    rotate(urX, X, rot2);
}

void getPolW(pcg128_t *polW, unsigned long long W0){ //OK !
    polW[0] = W0;
    for(int i = 1 ; i < nbiter ; i++){
        polW[i] = polW[i-1] * a;
    }
}

void getSumPol(unsigned long long* sumPol,unsigned long long* sumPolY, pcg128_t* polW){
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


/* sumPol/sumPolY are constant over many iterations. X (unrotated) and rot vary each time. */
int solve(pcg128_t* S, const unsigned long long* X, const int* rot, const unsigned long long* sumPol, const unsigned long long* sumPolY)
{
    unsigned long long Y[nbiter];
    double Yprim[nbiter];
    unsigned long long tmp3[nbiter];
    
    // total : 6 XOR, 6 AND, 6 SHIFT, 3 ADD, 3 SUB, 3 u64->double
    for (int i = 0 ; i < nbiter ; i++) {
        Y[i] = (((sumPol[i] ^ X[i]) % (1 << known_low)) << known_up ) + (rot[i] ^ (X[i] >> (k - known_up)));
        Yprim[i] = (Y[i] - sumPolY[i]) % (1 << (known_low + known_up));    
    }
    
    double tmp2[nbiter];
    // total : 9 MUL (double) + 10 ADD (double)
    for (int i=0 ; i<nbiter ; i++) {
        tmp2[i] = 0;
        for(int j=0 ; j<nbiter ; j++)
            tmp2[i] += invG[i * nbiter + j] * Yprim[j];
        tmp2[i] += 6755399441055744.0;
    }
 
    // total : 6 SHIFT
    for(int i = 0 ; i < nbiter ; i++)
        tmp3[i] = light_crazy_round(tmp2[i]);

    unsigned long long Sprim0 = 0;
    // total : 3 ADD, 3 MUL
    for(int j=0 ; j<nbiter ; j++)
        Sprim0 += Greduite[j] * tmp3[j];

    // SHIFT, ADD
    unsigned long long Smod = (Sprim0 << known_low) + sumPol[0];
    S[0] = (((pcg128_t)(Smod ^ X[0])) << k) + ((pcg128_t) Smod);

    // clocking PCG two times
    for (int i = 1 ; i < nbiter ; i++) {
        S[i] = S[i-1] * a + c;
        unsigned long long XX = S[i] ^ (S[i] >> 64);
        if (XX != X[i])
            return 0;
    }
    return 1;
}
