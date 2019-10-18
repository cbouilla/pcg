#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "fonctions.h"


unsigned long long Greduite[9] = {-1241281756092, -5001120657083,  8655886039732,
    3827459685972, -2117155768935,  3303731088004,
    -728312298332,  5479732607037,  6319848582548};

float invG[9] = {-9.25221813226351e-14,  2.32272588749499e-13,  5.30001389997814e-15,
    -7.81560146462246e-14, -4.52719459143047e-15,  1.09411828555506e-13,
    5.71040610630735e-14,  3.06929503514353e-14,  6.39750297008745e-14};

void init_var_globales(){
    //multiplier a OK !
    a = (((pcg128_t) 2549297995355413924) << k) + ((pcg128_t) 4865540595714422341);
    
    //increment c OK !
    c = (((pcg128_t)6364136223846793005) << k) + 1442695040888963407;
    
    //nombre de threads 
    nb_thread = omp_get_max_threads();
    
    //increment polynome polC OK !
    polC[0] = 0;
    pcg128_t powA =1;
    for(int i = 1; i < nbiter ; i++){
        polC[i] = polC[i-1] + powA * c;
        powA *= a;
    }
}


//////////////////// chrono //////////////////
double wtime()
{
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}



////////////////Fonctions calcul matriciel//////////////

//Matrice avce coeffs négatifs mais résult mod 2**64 (je sais pas si ça marche)
void prodMatVecU(unsigned long long* res, unsigned long long * M, unsigned long long* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

void prodMatVecFFU(float* res, float* M, unsigned long long* v, int n){
    int i, j;
    for(i=0 ; i<n ; i++){
        res[i] = 0;
        for(j=0 ; j<n ; j++)
            res[i]+= M[i * n + j] * v[j];
    }
}

void prodMatMatU(unsigned long long* res, unsigned long long* M1, unsigned long long* M2, int n){//juste pour verif
    int i, j, l;
    for(i=0 ; i<n ; i++){
        for(j=0 ; j<n ; j++){
            res[i * + j] = 0;
            for(l=0; l<n ; l++)
                res[i * n + j]+= M1[i * n + l] * M2[l * n + j];
        }
    }
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

void getY(unsigned long long* Y, unsigned long long* sumPol, int* rot, unsigned long long* X){ //OK !
    for(int i = 0 ; i < nbiter ; i++){
        Y[i] = (((sumPol[i] ^ X[i]) % (1 << known_low)) << known_up ) + (rot[i] ^ (X[i] >> (k - known_up)));
    }
}

void getYprim(unsigned long long* Yprim, unsigned long long* Y, unsigned long long* sumPolY){//OK ! avec erreurs arrondi
     //Utiliser sumPolY
    for(int i = 0 ; i < nbiter ; i++){
        Yprim[i] = Y[i] - sumPolY[i];
        Yprim[i] = Yprim[i] % (1 << (known_low + known_up));
    }
}

void findSprim(unsigned long long* Sprim, unsigned long long* Yprim){ //OK !
    unsigned long long tmp1[nbiter];
    for(int i = 0 ; i < nbiter ; i++)
        tmp1[i] = Yprim[i] << (k - known_low - known_up);
    float tmp2[nbiter];
    prodMatVecFFU(tmp2, invG, tmp1, nbiter);
    for(int i = 0 ; i < nbiter ; i++)
        tmp1[i] = (unsigned long long) llroundf(tmp2[i]);
    prodMatVecU(Sprim, Greduite, tmp1, nbiter);
}
 
void findS(pcg128_t* S, unsigned long long* Sprim, unsigned long long* X, unsigned long long* sumPol){
    unsigned long long Smod ;
    for(int i = 0 ; i < nbiter ; i++){
        Smod = (Sprim[i] << known_low) + sumPol[i];
        S[i] = (((pcg128_t)(Smod ^ X[i])) << k) + ((pcg128_t) Smod);
    }
}

int test(pcg128_t* S, unsigned long long* X){
    pcg128_t Si;
    Si = S[0];
    for(int i = 1 ; i < nbiter ; i++){
        Si = Si * a + c ; //mod 2^128 auto
        if(Si != S[i])
            return 0;
    }
    return 1;
}

int solve(pcg128_t* S, unsigned long long* X, int* rot,unsigned long long* sumPol,unsigned long long* sumPolY){
    unsigned long long Y[nbiter];
    getY(Y, sumPol, rot, X);
    unsigned long long Yprim[nbiter];
    getYprim(Yprim, Y, sumPolY);
    unsigned long long Sprim[nbiter];
    findSprim(Sprim, Yprim);
    findS(S, Sprim, X, sumPol);
    return test(S,X);
}

void pcg(pcg128_t* S, unsigned long long* X, pcg128_t S0, int n){
    struct pcg_state_128* rng = malloc(sizeof(struct pcg_state_128));
    pcg_oneseq_128_srandom_r(rng, S0);
    for(int i = 0 ; i < n ; i++){
        S[i] = rng->state;
        X[i] = pcg_output_xsl_rr_128_64(rng->state);
        pcg_oneseq_128_step_r(rng);
    }
}
        

