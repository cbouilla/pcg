#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>


int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    testFonctions();
    
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
}
