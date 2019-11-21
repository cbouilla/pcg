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
    pcg(vraiS, X, S0, nboutput);
    printf("%016llx %016llx\n", (unsigned long long) (vraiS[0]>>64), (unsigned long long) vraiS[0]);
    
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
        printf("%llu ", DS64[i]);
    printf("\n");
    
    /***** Vrai Test *****/
    int rot2[10 * 64];
    int nbrot[10];
    //FindRoti(roti,DS64[0], X[1], 1, Y[0], W0, WC);
    if(FindRot(rot2, nbrot, DS64[0], X, Y[0], W0, WC, 10)){
            printf("rot :\n");
        for(int i = 0 ; i < 10 ; i++)
            printf("%d ",rot2[k * i]);
        printf("\n");
    }
    else printf("On a pas trouvé :'(\n");
    
    return(0);
}