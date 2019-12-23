#include <assert.h>
#include <math.h>

#include "fonctions.h"

int testValid (FILE* f, int n) 
{
    int rot[nboutput];
    pcg128_t vraiS[nboutput];
    u64 X[nboutput];
    pcg128_t seeds[2];
    int cpt = 0;
    char* goodY = setupGoodY();

    for (int i = 0 ; i < n ; i++) {
        if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
            perror("Something went wrong when reading /dev/urandom");
            exit(EXIT_FAILURE);
        }
        
        pcg(vraiS, X, seeds[0], seeds+1, nboutput);
        
        for(int i = 0 ; i < nboutput ; i++)
            rot[i] = (int) (vraiS[i] >> (2 * k - known_up));

        u64 W0 = (u64) (vraiS[0] % (1<<known_low));
        u64 WC = (u64) (seeds[1] % (1<<known_low));
        
        /**** Polynômes en WC et W0 utilisés dans la résolution ****/
        u64 lowSumPol[nbiter + nbtest];
        u64 sumPolY[nbiter];
        u64 sumPolTest[nbtest];
        for(int i = 0 ; i < nbiter ; i++){
            lowSumPol[i]  = (W0 * ((u64) powA[i]) + WC * ((u64) polA[i]));
            sumPolY[i] = (polA[i] * WC + powA[i] * W0) >> (k - known_up);
        }
        for(int i = 0 ; i < nbtest ; i++){
            lowSumPol[nbiter + i] = (W0 * ((u64) powA[i + nbiter]) + WC * ((u64) polA[i + nbiter]));
            sumPolTest[i] = W0 * ((u64) (powA[i + nbiter] >> known_low) - 1) + WC * ((u64) (polA[i + nbiter] >> known_low) - 1);
        }

        u64 DS640;
        u64 Y0;
        /*u64 uX[nbiter];
        unrotateX(uX, X, rot);
        
        u64 Y[nbiter]; //utilisé dans testDS640
        getY(Y, W0, WC, rot, uX);
        FindDS64(DS64, Y, uX, rot, lowSumPol, sumPolY);
        cpt += testDS640(DS64[0], X, Y[0], sumPolTest, lowSumPol);*/

        u64 tabTmp[k * nbiter];
        getTabTmp(tabTmp, X, lowSumPol, sumPolY);    

        getGoodY(goodY, X, lowSumPol, 1);
        int a = solve_isgood(goodY, rot, tabTmp, sumPolY, sumPolTest); 
        cpt += a;
        if (a)
            solve(&DS640, &Y0, goodY, rot, tabTmp, sumPolY, sumPolTest);
        getGoodY(goodY, X, lowSumPol, 0);
    }
    return cpt;
}



int main()
{    
    /*  INITIALISATION DES PARAMETRES  */   
    init_var_globales();
    
    
    FILE *f = fopen("/dev/urandom", "r");
    if (f == NULL) {
        perror("Something went wrong when opening /dev/urandom");
        exit(EXIT_FAILURE);
    }
    
    printf("1..1\n");
    static const int nbtests = 1000000;
    int successes = testValid(f, nbtests);
    if (successes != nbtests)
        printf("not ok 1 - #success = %d / %d\n", successes, nbtests);
    else
        printf("ok 1 - all %d tests OK\n", nbtests);

    exit(0);
}
