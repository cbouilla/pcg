#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include "entropy.h"//version un peu modifiée pour ne PAS passer par pcg_variant.h parce que c'est déjà assez chiant comme ça !

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    pcg128_t seeds[2];
    entropy_getbytes((void*)seeds, sizeof(seeds)); // proposé et utilisé dans le code de pcg
    
    testValidite();
    //pcg128_t S0 = (((pcg128_t) 5995207026785010249) << k) + ((pcg128_t) 179350442155841024);
    //pcg128_t c = ((((pcg128_t) 6364136223846793005) << k) + 1442695040888963407) >> 1;
    //printVal(S0, c);
    for(int i = 0 ; i < 1000 ; i++)
        if(rand() > ((unsigned long long) 1) << 31)
            printf("Hey !\n");
    return(0);
}

+    // Seed the full-blown generator with external entropy
+    pcg128_t seeds[2];
+    
+    FILE *f = fopen("/dev/urandom", "r");
+    if (f == NULL) {
+        perror("Something went wrong when opening /dev/urandom");
+        exit(EXIT_FAILURE);
+    }
+    if (fread(seeds, sizeof(seeds), 1, f) != 1)  {
+        perror("Something went wrong when reading /dev/urandom");
+        exit(EXIT_FAILURE);
+    }
+
+    printf("Seed[0] : %016" PRIx64 " %016" PRIx64 "\n", (uint64_t) (seeds[0] >> 64), (uint64_t) seeds[0]);
+    printf("Seed[0] : %016" PRIx64 " %016" PRIx64 "\n\n", (uint64_t) (seeds[1] >> 64), (uint64_t) seeds[1]);
+
