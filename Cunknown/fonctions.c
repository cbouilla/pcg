#include <assert.h>
#include <time.h>
#include <sys/time.h>

#include "fonctions.h"

unsigned long long Greduite[16] =
{-186304953996472, -216211368070119,  110964501361298,  131252974561432,
-126056243766680,   99587582169277,   -5646098666150, -233919070109448,
   7937589136904, -214303762177807, -268280113597118,  -98716819647784,
  93078431381544,   -1707551230219,  149382085707466, -134620659538888}; 

double invG[16] =
{-2.04279952328856e-15, -2.93791683689260e-15,  6.74861642263602e-16,  2.61840245666112e-15,
-1.98886046520741e-15,  1.10621147658696e-15, -2.12731308263381e-15, -2.30132753131773e-15,
 1.44762674070265e-15, -1.54769860122306e-16, -1.55490854567009e-15,  2.82055195172104e-15,
 2.19171469167945e-16, -2.21708502251561e-15, -1.23181629826941e-15, -2.45886210697988e-15};


void init_var_globales()
{
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
    for (int i = 0; i < 16; i++)
    	invG[i] *= 1ull << (k - known_up - known_low);
}


//////////////////// chrono //////////////////
double wtime()
{
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

char * setupGoodY()
{
	char* goodY = malloc((1<<(known_low + known_up)) * sizeof(char) * nbtest / 8);
	for (unsigned long long y = 0 ; y < nbtest * (1<< (known_low + known_up)) / 8 ; y++)
	    goodY[y] = 0;
	return goodY;
}

static inline void setbit(char *goodY, int idx, int v)
{
	// idx = idx / 4;
	int j = idx / 8;
	int l = idx % 8;
	if (v == 1)
		goodY[j] |= (1 << l);
	else
		goodY[j] &= ~(1 << l);
}

void getGoodY(char* goodY, unsigned long long* tabX, unsigned long long* lowSumPol, int v)
{
	for (int i = 0 ; i < nbtest ; i++){
		unsigned long long Wi = lowSumPol[nbiter + i] % (1 << known_low);
	    for (int j = 0 ; j < k ; j++){
	        unsigned long long Xij = tabX[i*k + j]; //unrotate(X[i], j);
	        unsigned long long goodYi1 = (((Xij % (1 << known_low)) ^ Wi) << known_up) ^ (j ^ (Xij >> (k - known_up)));
	        unsigned long long goodYi2 = (goodYi1 - 1) % (1 << (known_low + known_up));
	        setbit(goodY, goodYi1 + i * (1<<(known_low + known_up)), v);
	        setbit(goodY, goodYi2 + i * (1<<(known_low + known_up)), v);
	    }
	}
}


void getTabTmp(unsigned long long* tabTmp, unsigned long long* X, unsigned long long* lowSumPol, unsigned long long* sumPolY)
{
	for (int i = 0 ; i < k ; i++){
		for (int j = 0 ; j < nbiter ; j++){
    	    unsigned long long uX = unrotate(X[j], i);
			tabTmp[i * nbiter + j] = (((lowSumPol[j] % (1 << known_low)) ^ (uX % (1 << known_low))) << known_up) + (i ^ (uX >> (k - known_up))) - sumPolY[j];
		}
	}
}

static inline int checkY(const char* goodY, int i, unsigned long long Y)
{
	int idx = Y + i * (1 << (known_up + known_low));
	// idx = idx / 4;
	int j = idx / 8;
	int l = idx % 8;
	return (goodY[j] >> l) & 1;
}

static inline bool confirm(unsigned long long Y0, unsigned long long DS640, const unsigned long long* sumPolTest, const char* goodY)
{
	/**** Confirmation du DS640 ****/
    for (int i = 0 ; i < nbtest ; i++) {
        unsigned long long tmp2 = ((unsigned long long) polA[i + nbiter]) * DS640 + sumPolTest[i]; //ATTENTION cast pcg128_t
        unsigned long long Yi1 = (Y0 + (tmp2 >> (k - known_low - known_up))) % (1 << (known_low + known_up)); //avec ou sans retenue OK!
        if (!(checkY(goodY, i, Yi1)))
            return 0;
    }
    return 1;
}

static inline long long crazy_round(double x)
{
    union { double d; long long l; } magic; 
    magic.d = x + 6755399441055744.0; 
    magic.l <<= 13; 
    magic.l >>= 13;
    return magic.l;
}

bool solve_isgood(const char* goodY, const int* rot, const unsigned long long* tabTmp, const unsigned long long* sumPolY, const unsigned long long* sumPolTest)
{
    /**** Recherche du DS640 ****/
    unsigned long long tmp[nbiter];
    for (int i = 0; i < nbiter; i++) //Y
        tmp[i] = tabTmp[i + nbiter * rot[i]];  
    
    unsigned long long Y0 = tmp[0] + sumPolY[0];
    
    double tmp3[nbiter - 1];
    for (int i = 0; i < nbiter - 1; i++)  //DY
        tmp3[i] = (tmp[i+1] - tmp[i]) % (1 << (known_low + known_up));
    
    double u[nbiter-1];
    for (int i=0 ; i<nbiter-1 ; i++) {
        u[i] = 0.0;
        for (int j=0 ; j<nbiter-1 ; j++)
            u[i] += invG[i * (nbiter-1) + j] * tmp3[j];
    }

    unsigned long long DS640 = 0;
    for(int i = 0 ; i < nbiter-1 ; i++) {
        DS640 += Greduite[i] * crazy_round(u[i]);
    }
  
 	return confirm(Y0, DS640, sumPolTest, goodY);
}


void solve(unsigned long long* DS640, unsigned long long* Y0, char* goodY, int* rot, unsigned long long* tabTmp, unsigned long long* sumPolY, unsigned long long* sumPolTest)
{
    unsigned long long tmp[nbiter];

    /**** Recherche du DS640 ****/

    for (int i = 0; i < nbiter; i++) //Y
        tmp[i] = tabTmp[i + nbiter * rot[i]];  
    
    *Y0 = (tmp[0] + sumPolY[0]) % (1 << (known_low + known_up));
    
    unsigned long long tmp3[nbiter - 1];
    for(int i = 0; i < nbiter - 1; i++)  //DY
        tmp3[i] = (tmp[i+1] - tmp[i]) % (1 << (known_low + known_up));
    
    double u[nbiter-1];
    for (int i=0 ; i<nbiter-1 ; i++) {
        u[i] = 0.0;
        for (int j=0 ; j<nbiter-1 ; j++)
            u[i]+= invG[i * (nbiter-1) + j] * tmp3[j];
    }

    *DS640 = 0;
    for(int i = 0 ; i < nbiter-1 ; i++) {
    	(*DS640) += Greduite[i] * crazy_round(u[i]);
    }
  
    assert(confirm(*Y0, *DS640, sumPolTest, goodY));
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

