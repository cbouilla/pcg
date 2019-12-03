#include <assert.h>
#include <math.h>

#include "fonctions.h"

int testValid (FILE* f, int nbtest) 
{
    int rot[nboutput];
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    pcg128_t seeds[2];
    int cpt = 0;

    for (int i = 0 ; i < nbtest ; i++) {
        if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
            perror("Something went wrong when reading /dev/urandom");
            exit(EXIT_FAILURE);
        }
        
        pcg(vraiS, X, seeds[0], seeds+1, nboutput);
        
        for(int i = 0 ; i < nboutput ; i++)
            rot[i] = (int) (vraiS[i] >> (2 * k - known_up));

        unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
        unsigned long long WC = (unsigned long long) (seeds[1] % (1<<known_low));
        unsigned long long uX[nbiter];
        unrotateX(uX, X, rot);
        
        unsigned long long Y[nbiter]; //utilisé dans testDS640
        getY(Y, W0, WC, rot, uX);
        
        unsigned long long DS64[nboutput];
        FindDS64(DS64, uX, rot, W0, WC);
        cpt += testDS640(DS64[0], X, Y[0], W0, WC, 3);
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
    /*if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
        perror("Something went wrong when reading /dev/urandom");
        exit(EXIT_FAILURE);
    }

    printf("Seed[0] : %016" PRIx64 " %016" PRIx64 "\n", (uint64_t) (seeds[0] >> 64), (uint64_t) seeds[0]);
    printf("Seed[1] : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seeds[1] >> 64), (uint64_t) seeds[1]);
    */
    
    //pcg128_t S0 = (((pcg128_t) 5995207026785010249) << k) + ((pcg128_t) 179350442155841024);
    //pcg128_t c = ((((pcg128_t) 6364136223846793005) << k) + 1442695040888963407) >> 1;
    //printVal(seeds[0], seeds[1]);
    
    printf("1..1\n");
    static const int nbtests = 1000000;
    int successes = testValid(f, nbtests);
    if (successes != nbtests)
        printf("not ok 1 - #success = %d / %d\n", successes, nbtests);
    else
        printf("ok 1 - all %d tests OK\n", nbtests);

    exit(0);
}
