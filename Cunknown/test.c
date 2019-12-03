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
    for(int i = 0 ; i < nbtest ; i++){
        if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
            perror("Something went wrong when reading /dev/urandom");
            exit(EXIT_FAILURE);
        }
        
        pcg(vraiS, X, seeds[0], seeds+1, nboutput);
        
        for(int i = 0 ; i < nboutput ; i++){
            rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
        }
        unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
        unsigned long long WC = (unsigned long long) (seeds[1] % (1<<known_low));
        unsigned long long uX[nbiter];
        unrotateX(uX, X, rot);
        
        unsigned long long Y[nbiter];//utilisé dans testDS640
        getY(Y, W0, WC, rot, uX);
        
        unsigned long long DS64[nboutput];
        FindDS64(DS64, uX, rot, W0, WC);
        cpt += testDS640(DS64[0], X, Y[0], W0, WC, 3);
    }
    return cpt;
}


int testFonctions()
{
    assert(known_low == 11);

    int rot[nboutput];
    pcg128_t S0 = (((pcg128_t) 5995207026785010249u) << k) + ((pcg128_t) 179350442155841024u);
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    pcg128_t c = ((((pcg128_t) 6364136223846793005u) << k) + 1442695040888963407u) >> 1;
    pcg(vraiS, X, S0, &c, nboutput);
    if(vraiS[2] != (((pcg128_t) 1792771836637573954u) << k) + ((pcg128_t) 11139816115278170276u)){
        printf("erreur sur pcg\n");
        return 0;
    }
    for(int i = 0 ; i < nboutput ; i++){
        rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
    }
    
    unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
    unsigned long long WC = (unsigned long long) (c % (1<<known_low));
 
    printf("1..5\n");
    printf("# nbiter = %d\n", nbiter);
    printf("# known_low = %d\n", known_low);

    /* test unrotate */
    unsigned long long uX[nbiter];
    unrotateX(uX, X, rot);
    if (uX[0] != 15007519919903780682u) {
        printf("not ok 1 - erreur sur unrotateX\n");
    } else {
        printf("ok 1 - unrotate\n");
    }

    /* test getY */
    unsigned long long Y[nbiter];
    printf("# W0 = %llu\n", W0);
    printf("# Wc = %llu\n", WC);
    printf("# rot = [");
    for (int i = 0 ; i < nboutput ; i++)
        printf("%d; ", rot[i]);
    printf("]\n");
    printf("# uX = [");
    for (int i = 0 ; i < nbiter ; i++)
        printf("%llu; ", uX[i]);
    printf("]\n");

    getY(Y, W0, WC, rot, uX);
    if(Y[3] != 129714){
        printf("not ok 2 - erreur sur getY. Y[3] == %llu / attendu : 129714\n", Y[3]);
    } else {
        printf("ok 2 - getY\n");
    }

    /* test getYprim */
    unsigned long long Yprim[nbiter];
    getYprim(Yprim, Y, W0, WC);
    if(Yprim[1] != 93486){
        printf("not ok 3 - erreur sur getY\n");
    }  else {
        printf("ok 3 - getYprim\n");
    }

    /* test getDY */
    unsigned long long DY[nbiter-1];
    getDY(DY, Yprim);
    if(DY[0] != 14609){
        printf("not ok 4 - erreur sur getDY\n");
    }  else {
        printf("ok 4 - getDY\n");
    }

    /* test FindDS64 */
    unsigned long long DS64[nbiter - 1];
    FindDS64(DS64, uX, rot, W0, WC);
    if(DS64[0] != 2055999906439120392u){
        printf("not ok 5 - erreur sur FindDS64\n");
    }  else {
        printf("ok 5 - FindDS64\n");
    }
    /*int rot2[nboutput * 64];
    int nbrot[nboutput];
    if(!FindRot(rot2, nbrot, DS64[0], X, Y[0], W0, WC, nboutput)){
        printf("erreur sur FindRot\n");
        return 0;
    }
    for(int i = 0 ;i < nboutput ; i++){
        int test = 0;
        for(int j = 0 ; j < nbrot[i] ; j++)
            if(rot2[i * k + j] == rot[i])
                test = 1;
        if(test == 0){
            printf("erreur sur FindRot\n");
            return 0;
        }
    }*/
    return 1;
}



void printVal(pcg128_t S0, pcg128_t c){
    int rot[nboutput];
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    pcg(vraiS, X, S0, &c, nboutput);
    printf("setseq : %llu %llu\n", (unsigned long long) (vraiS[2]>>64), (unsigned long long) vraiS[2]);
    
    //printf("%llu %llu\n", (unsigned long long) ((vraiS[1] - vraiS[0])>>64), (unsigned long long) (vraiS[1] - vraiS[0]));
    
    //done
    for(int i = 0 ; i < nboutput ; i++){
        rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
    }
    
    unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
    unsigned long long WC = (unsigned long long) (c % (1<<known_low));
    printf("W0 : %llu\n", W0);
    printf("WC : %llu\n", WC);
    
    unsigned long long uX[nbiter];
    unrotateX(uX, X, rot);
    printf("uX\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", uX[i]);
    printf("\n");
    /*FILE *f;
    f = fopen("result.txt","w");
    
    fprintf(f,"W0 : %llu\n", W0);*/
    unsigned long long Y[nbiter];
    getY(Y, W0, WC, rot, uX);
    printf("Y :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", Y[i]);
    printf("\n");
    
    
    /***** Tests de vérification des sous-fonctions *****/
    unsigned long long Yprim[nbiter];
    getYprim(Yprim, Y, W0, WC);
    printf("Yprim :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", Yprim[i]);
    printf("\n");
    
    unsigned long long DY[nbiter];
    getDY(DY, Yprim);
    printf("DY :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%llu ", DY[i]);
    printf("\n");
    
    unsigned long long DS64[nbiter - 1];
    FindDS64(DS64, uX, rot, W0, WC);
    printf("DS64 :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%llu ", DS64[i]);
    printf("\n");
    
    printf("testDS640 : %d\n", testDS640(DS64[0], X, Y[0], W0, WC, 3)); 
    
}

void print128(pcg128_t x)
{
    printf("%016" PRIx64 " %016" PRIx64, (uint64_t) (x >> 64), (uint64_t) x);
}


int main()
{    
    /*  INITIALISATION DES PARAMETRES  */   
    init_var_globales();
    
    for (int i = 0 ; i < nboutput ; i++) {
        printf("# polA[%2d] =", i);
        print128(polA[i]);
        printf("\n");
    }

    for (int i = 0 ; i < nboutput ; i++) {
        printf("# powA[%2d] =", i);
        print128(powA[i]);
        printf("\n");
    }

    testFonctions();
    exit(0);

#if 0
    pcg128_t seeds[2];
    
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
    
    printf("test de validité de taille 1000 : %d\n", testValid(f, 1000));
    return 0;
#endif
}
