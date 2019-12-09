#include <assert.h>
#include <math.h>

#include "fonctions.h"

int testValid (FILE* f, int n) 
{
    int rot[nboutput];
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    pcg128_t seeds[2];
    int cpt = 0;

    for (int i = 0 ; i < n ; i++) {
        if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
            perror("Something went wrong when reading /dev/urandom");
            exit(EXIT_FAILURE);
        }
        
        pcg(vraiS, X, seeds[0], seeds+1, nboutput);
        
        for(int i = 0 ; i < nboutput ; i++)
            rot[i] = (int) (vraiS[i] >> (2 * k - known_up));

        unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
        unsigned long long WC = (unsigned long long) (seeds[1] % (1<<known_low));
        
        /**** Polynômes en WC et W0 utilisés dans la résolution ****/
        unsigned long long lowSumPol[nbiter + nbtest];
        unsigned long long sumPolY[nbiter];
        unsigned long long sumPolTest[nbtest];
        for(int i = 0 ; i < nbiter ; i++){
            lowSumPol[i]  = (W0 * ((unsigned long long) powA[i]) + WC * ((unsigned long long) polA[i]));
            sumPolY[i] = (polA[i] * WC + powA[i] * W0) >> (k - known_up);
        }
        for(int i = 0 ; i < nbtest ; i++){
            lowSumPol[nbiter + i] = (W0 * ((unsigned long long) powA[i + nbiter]) + WC * ((unsigned long long) polA[i + nbiter]));
            sumPolTest[i] = W0 * ((unsigned long long) (powA[i + nbiter] >> known_low) - 1) + WC * ((unsigned long long) (polA[i + nbiter] >> known_low) - 1);
        }

        unsigned long long DS640;
        unsigned long long Y0;
        /*unsigned long long uX[nbiter];
        unrotateX(uX, X, rot);
        
        unsigned long long Y[nbiter]; //utilisé dans testDS640
        getY(Y, W0, WC, rot, uX);
        FindDS64(DS64, Y, uX, rot, lowSumPol, sumPolY);
        cpt += testDS640(DS64[0], X, Y[0], sumPolTest, lowSumPol);*/
        cpt += solve(&DS640, &Y0, X, rot, lowSumPol, sumPolY, sumPolTest);
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