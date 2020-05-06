#include <assert.h>
#include <time.h>
#include <sys/time.h>

#include "fonctions.h"

pcg128_t a;
pcg128_t a_inv;
pcg128_t powA[nboutput];
pcg128_t polA[nboutput];


u64 Greduite[16] =
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

    a_inv = (((pcg128_t) 566787436162029664) << k) + ((pcg128_t) 11001107174925446285);

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

static inline u64 unrotate(u64 Xi, int i)
{
    return (Xi >> (k-i)) | (Xi << i);
}

char * setupGoodY()
{
	char* goodY = malloc((1<<(known_low + known_up)) * sizeof(char) * nbtest / 8);
	for (u64 y = 0 ; y < nbtest * (1<< (known_low + known_up)) / 8 ; y++)
	    goodY[y] = 0;
	return goodY;
}

static inline void setbit(char *goodY, int i, u64 Y, int v)
{
    int idx = Y + i * (1 << (known_up + known_low));
	int j = idx / 8;
	int l = idx % 8;
	if (v == 1)
		goodY[j] |= (1 << l);
	else
		goodY[j] &= ~(1 << l);
}

void getGoodY(char* goodY, const u64 (*X)[nboutput], const u64* lowSumPol, int v)
{
	for (int i = 0 ; i < nbtest ; i++){
		u64 Wi = lowSumPol[nbiter + i] % (1 << known_low);
	    for (int j = 0 ; j < k ; j++){
	        u64 Xij = unrotate((*X)[i + nbiter], j);
	        u64 goodYi1 = (((Xij % (1 << known_low)) ^ Wi) << known_up) ^ (j ^ (Xij >> (k - known_up)));
	        u64 goodYi2 = (goodYi1 - 1) % (1 << (known_low + known_up));
	        setbit(goodY, i, goodYi1, v);
	        setbit(goodY, i, goodYi2, v);
	    }
	}
}


static inline int checkY(const char* goodY, int i, u64 Y)
{
	// Y = Y(1 << (known_up + known_low - 1));
	int idx = Y + i * (1 << (known_up + known_low));
	// idx = idx / 4;
	int j = idx / 8;
	int l = idx % 8;
	return (goodY[j] >> l) & 1;
}

/* THIS is also time-critical */
static inline bool confirm(u64 Y0, u64 DS640, const struct task_t *task)
{
	/**** Confirmation du DS640 ****/
    // <= 8 ADD, 4 MUL, 8 SHIFT, 8 AND
    for (int i = 0 ; i < nbtest ; i++) {
        u64 tmp2 = ((u64) polA[i + nbiter]) * DS640 + task->sumPolTest[i]; //ATTENTION cast pcg128_t
        u64 Yi1 = (Y0 + (tmp2 >> (k - known_low - known_up))) % (1 << (known_low + known_up)); //avec ou sans retenue OK!
        if (!(checkY(task->goodY, i, Yi1)))
            return 0;
    }
    return 1;
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

/* THIS is the time-critical function */
bool solve_isgood(const struct task_t *task)
{
    /**** Recherche du DS640 ****/
    u64 tmp[nbiter];
    for (int i = 0; i < nbiter; i++) //Y
        tmp[i] = task->tabTmp[i + nbiter * task->rot[i]];  
    
    // ADD
    u64 Y0 = tmp[0] + task->sumPolY[0];
    
    // 4 SUB, 4 AND, 4 u64->double
    double tmp3[nbiter - 1];
    for (int i = 0; i < nbiter - 1; i++)  //DY
        tmp3[i] = (tmp[i+1] - tmp[i]) % (1 << (known_low + known_up));
    
    // 20 ADD (double), 16 MUL (double)
    double u[nbiter-1];
    for (int i = 0; i < nbiter - 1; i++) {
        u[i] = 0.0;
        for (int j = 0; j<nbiter - 1; j++)
            u[i] += invG[i * (nbiter-1) + j] * tmp3[j];
        u[i] += 6755399441055744.0;
    }

    // 8 SHIFT, 4 MUL, 4 ADD
    u64 DS640 = 0;
    for (int i = 0; i < nbiter-1; i++) {
        DS640 += Greduite[i] * light_crazy_round(u[i]);
    }
  
    // 8 ops dans la plupart des cas
 	return confirm(Y0, DS640, task);
}


void solve(const struct task_t *task, u64* DS640, u64* Y0)
{
    u64 tmp[nbiter];

    /**** Recherche du DS640 ****/

    for (int i = 0; i < nbiter; i++) //Y
        tmp[i] = task->tabTmp[i + nbiter * task->rot[i]];  
    
    *Y0 = (tmp[0] + task->sumPolY[0]) % (1 << (known_low + known_up));
    
    u64 tmp3[nbiter - 1];
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
  
    assert(confirm(*Y0, *DS640, task));
}


void pcg(pcg128_t *S, u64* X, pcg128_t S0, pcg128_t* c, int n)
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

void init_task(struct task_t *t)
{
    t->goodY = malloc((1<<(known_low + known_up)) * sizeof(char) * nbtest / 8);
    for (u64 y = 0 ; y < nbtest * (1<< (known_low + known_up)) / 8 ; y++)
        t->goodY[y] = 0;
}

void prepare_task(const u64 (*X)[nboutput], u64 W0, u64 WC, struct task_t *t)
{
    for (int i = 0 ; i < nbiter ; i++) {
        t->lowSumPol[i]  = (W0 * ((u64) powA[i]) + WC * ((u64) polA[i]));
        t->sumPolY[i] = (polA[i] * WC + powA[i] * W0) >> (k - known_up);
    }
    for (int i = 0 ; i < nbtest ; i++) {
        t->lowSumPol[nbiter + i] = (W0 * ((u64) powA[i + nbiter]) + WC * ((u64) polA[i + nbiter]));
        t->sumPolTest[i] = W0 * ((u64) (powA[i + nbiter] >> known_low) - 1) + WC * ((u64) (polA[i + nbiter] >> known_low) - 1);
    }

    getGoodY(t->goodY, X, t->lowSumPol, 1);

    // getTabTmp(t->tabTmp, X, t->lowSumPol, t->sumPolY);
    for (int i = 0 ; i < k ; i++) {
        for (int j = 0 ; j < nbiter ; j++) {
            u64 uX = unrotate((*X)[j], i);
            t->tabTmp[i * nbiter + j] = (((t->lowSumPol[j] % (1 << known_low)) ^ (uX % (1 << known_low))) << known_up) + (i ^ (uX >> (k - known_up))) - t->sumPolY[j];
        }
    }
    
    for (int i = 0 ; i < nbiter ; i++)
        t->rot[i] = 0;
}

void finish_task(const u64 (*X)[nboutput], struct task_t *t)
{
    getGoodY(t->goodY, X, t->lowSumPol, 0);
}