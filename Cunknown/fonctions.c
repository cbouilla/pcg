#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "fonctions.h"

unsigned long long Greduite[16] =
{-186304953996472, -216211368070119,  110964501361298,  131252974561432,
-126056243766680,   99587582169277,   -5646098666150, -233919070109448,
   7937589136904, -214303762177807, -268280113597118,  -98716819647784,
  93078431381544,   -1707551230219,  149382085707466, -134620659538888}; 

float invG[16] =
{-2.04279952328856e-15, -2.93791683689260e-15,  6.74861642263602e-16,  2.61840245666112e-15,
-1.98886046520741e-15,  1.10621147658696e-15, -2.12731308263381e-15, -2.30132753131773e-15,
 1.44762674070265e-15, -1.54769860122306e-16, -1.55490854567009e-15,  2.82055195172104e-15,
 2.19171469167945e-16, -2.21708502251561e-15, -1.23181629826941e-15, -2.45886210697988e-15};


void init_var_globales(){
    //multiplier a OK !
    a = (((pcg128_t) 2549297995355413924) << k) + ((pcg128_t) 4865540595714422341);

    //nombre de threads 
    nb_thread = omp_get_max_threads();
    
    //increment polynome polC OK !
    polA[0] = 0;
    powA[0] = 1;
    for(int i = 1; i < nboutput ; i++){
        polA[i] = polA[i-1] + powA [i-1];
        powA[i] = powA[i-1] * a;
    }
    for(int i = 0; i < nboutput ; i++){
        lowPolA[i] = polA[i];
        lowPowA[i] = powA[i];
        YZPowA[i] = powA[i] >> known_low;
        YZPolA[i] = polA[i] >> known_low;
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


void rotateX(unsigned long long* rX, unsigned long long* X,int* rot){ //pas verifié, repris de pcg_random
    for(int i = 0 ; i < nbiter ; i++)
        rX[i]= (X[i] >> rot[i]) | (X[i] << ((- rot[i]) & (k-1)));
}

void unrotateX(unsigned long long* urX, unsigned long long* X, int* rot){//pas verifié, repris de pcg_random
    int rot2[nbiter];
    for( int i = 0 ; i < nbiter ; i++)
        rot2[i] = (k - rot[i]) % k;
    rotateX(urX, X, rot2);
}

unsigned long long unrotate1(unsigned long long Xi){
    return (Xi >> (k-1)) | (Xi << 1);
}

void getY(unsigned long long *Y, unsigned long long W0, unsigned long long WC, int* rot, unsigned long long* uX){
    for(int i = 0 ; i < nbiter ; i++){
        Y[i] = ((((unsigned long long) ((polA[i] * WC + powA[i] * W0) % (1 << known_low))) ^ (uX[i] % (1 << known_low))) << known_up) + (rot[i] ^ (uX[i] >> (k - known_up)));
    }
}

void getYprim(unsigned long long *Yprim, unsigned long long *Y, unsigned long long W0, unsigned long long WC){
    for(int i = 0 ; i < nbiter ; i++)
        Yprim[i] = (unsigned long long) (Y[i] - ((polA[i] * WC + powA[i] * W0) >> (k - known_up))) % (1<<(known_up + known_low));
}

void getDY(unsigned long long *DY, unsigned long long* Yprim){
    for(int i = 0 ; i < nbiter - 1 ; i++)
        DY[i] = (Yprim[i+1] - Yprim[i]) % (1<<(known_low + known_up));
}

void FindDS64(unsigned long long* DS64,unsigned long long* uX,int* rot,unsigned long long W0,unsigned long long WC){
    unsigned long long Y[nbiter];
    getY(Y, W0, WC, rot, uX);
    
    unsigned long long Yprim[nbiter];
    getYprim(Yprim, Y, W0,WC);
    
    
    unsigned long long tmp[nbiter-1];
    getDY(tmp, Yprim);
    
    for(int i = 0 ; i < nbiter - 1 ; i++)
        tmp[i] = tmp[i] << (k - known_up - known_low);
    float u[nbiter-1];
    prodMatVecFFU(u, invG, tmp, nbiter-1);
    for(int i = 0 ; i < nbiter-1 ; i++)
        tmp[i] = (unsigned long long) llroundf(u[i]);
    prodMatVecU(DS64, Greduite, tmp, nbiter-1);
}
    

    
int testDS640(unsigned long long DS640,  unsigned long long* X, unsigned long long Y0, unsigned long long W0, unsigned long long WC, int n){
    for(int i = nbiter ; i < n + nbiter ; i++){
        unsigned long long Xi = X[i];
        unsigned long long tmp = polA[i] * DS640; //ATTENTION cast pcg128_t
        tmp += W0 * ((unsigned long long) (powA[i] >> known_low) - 1) + WC * ((unsigned long long) (polA[i] >> known_low) - 1);//ATTENTION cast pcg128_t
        unsigned long long Yi1 = (Y0 + (tmp >> (k - known_low - known_up))) % (1<< (known_low + known_up)); //avec ou sans retenue OK!
        unsigned long long Yi2 = Yi1 + 1;
        unsigned long long Wi = (W0 * ((unsigned long long) powA[i]) + WC * ((unsigned long long) polA[i])) % (1<< known_low);//ATTENTION cast pcg128_t
        int test1, test2,test = 0;
        for(int i = 0 ; i < k ; i++){
            test1 = (((Xi ^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) && ((i ^ (Xi >> (k - known_up))) == Yi1 % (1 << known_up));
            test2 = (((Xi ^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) && ((i ^ (Xi >> (k - known_up))) == Yi2 % (1 << known_up));
            if (test1 || test2){
                test = 1;
            }
            Xi = unrotate1(Xi);
        }
        if(!test) {
            //printf("erreur a i = %d\n", i);
            return 0;
        }
    }
    return 1;
    
}

void pcg(pcg128_t *S, unsigned long long* X, pcg128_t S0, pcg128_t* c, int n)
{
    struct pcg_state_setseq_128 rng;
    pcg_setseq_128_srandom_r(&rng, S0, *c);
    for (int i = 0 ; i < n ; i++) {
        X[i] = pcg_output_xsl_rr_128_64(rng.state);
        S[i] = rng.state;
        pcg_setseq_128_step_r(&rng);
    }
    *c = rng.inc;
}
