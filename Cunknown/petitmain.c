#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    
    
    int rot[nboutput];
    pcg128_t S0 = (((pcg128_t) 5995207026785010249) << k) + ((pcg128_t) 179350442155841024);
    pcg128_t vraiS[nboutput];
    unsigned long long X[nboutput];
    pcg(vraiS, X, S0, c, nboutput);
    printf("%016llx %016llx\n", (unsigned long long) (vraiS[0]>>64), (unsigned long long) vraiS[0]);
    
    printf("%llu %llu\n", (unsigned long long) ((vraiS[1] - vraiS[0])>>64), (unsigned long long) (vraiS[1] - vraiS[0]));
    for(int i = 0 ; i < nboutput ; i++){
        rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
    }
    
    unsigned long long W0 = (unsigned long long) (vraiS[0] % (1<<known_low));
    unsigned long long WC = (unsigned long long) (c % (1<<known_low));
    
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
        printf("%016llx ", DS64[i]);
    printf("\n");
    
    pcg128_t vraiDS64[nboutput - 1];
    getDS64(vraiDS64, vraiS, W0, WC);
    printf("vraiDS64 :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%016llx\n", (unsigned long long) (vraiDS64[i]));
    printf("\n");
    
    pcg128_t DSmod[nbiter - 1];
    getDSmod(DSmod, DS64, W0, WC);
    printf("DSmod :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%016llx %016llx \n", (unsigned long long) (DSmod[i]>>64), (unsigned long long) DSmod[i]);
    printf("\n");
    
    /***** Vrai Test *****/
    int rot2[nboutput * 64];
    int nbrot[nboutput];
    //FindRoti(roti,DS64[0], X[1], 1, Y[0], W0, WC);
    if(FindRot(rot2, nbrot, DS64[0], X, Y[0], W0, WC, nboutput)){
            printf("rot :\n");
        for(int i = 0 ; i < nboutput ; i++)
            printf("%d ",rot2[k * i]);
        printf("\n");
    }
    else printf("On a pas trouvé :'(\n");
    int rotDS[nboutput-1];
    FindRotDS(rotDS, rot2, DS64[0], W0, WC);
    unsigned long long upDS[nboutput - 1];
    pcg128_t DS[nboutput - 1];
    FindUpDS(upDS, rotDS);
    FindDS(DS, upDS, DS64[0], W0, WC);
    for(int i = 0 ; i < nboutput - 1 ; i++)
        printf("%016llx %016llx \n", (unsigned long long) (DS[i]>>64), (unsigned long long) DS[i]);
    return(0);
}