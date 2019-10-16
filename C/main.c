#include "fonctions.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(){
    
    /*  INITIALISATION DES PARAMETRES  */
    
    init_var_globales();
        
    /********** Calculs/Tests plus ou moins à la con ***********/
    pcg128_t vraiS[nbiter];
    unsigned long long X[nbiter]; //= {15511475948252377726llu, 6234043309952295463llu, 9845631557281214488llu};
    pcg128_t S0 = (((pcg128_t) 5995207026785010249) << k) + ((pcg128_t) 179350442155841024);
    pcg(vraiS, X, S0, nbiter);
    printf("W0 : %llu, rot = %d, %d, %d \n", (unsigned long long) (vraiS[0]% (1<<known_low)), (int) (vraiS[0] >> (2*k - known_up)),(int) (vraiS[1] >> (2*k - known_up)),(int) (vraiS[2] >> (2*k - known_up)));
    unsigned long long W0;
    int rot[nbiter];
    unsigned long long sumPol[nbiter];
    unsigned long long sumPolY[nbiter];
    pcg128_t S[nbiter];
    pcg128_t polW[nbiter];
    
    FILE *f;
    f = fopen("result.txt","w");
    
    double t1 = wtime();
    int done = 0;
    
    for (W0 = 209817 ; W0 < (1<<known_low) ; W0++){//(1<<known_low)//W0=209818
        if((W0 % (1<<10)) == 0){
            printf("W0 : %llu (%.1f / s)\n", W0, done / (wtime() - t1));
            fprintf(f, "W0 : %llu\n", W0);
            fflush(f);
        }
        getPolW(polW, W0);
        getSumPol(sumPol,sumPolY, polW);
        for(int r = 0 ; r < 1<<(3*known_up) ; r++){
            //modif de rot :
            rot[0]=(rot[0] + 1) % k;
            int i = 0;
            while(rot[i] == 0 && i < nbiter){
                i++;
                rot[i]=(rot[i] + 1) %k;
            }
            unrotate(X, rot); //a optimiser (en meme temps que la modif de rot)
            //decaler X à rajouter avec vraie sortie de pcg
            
            if(solve(S, X, rot,sumPol,sumPolY)){
                fprintf(f,"S :\n");
                printf("S trouvé !!\n");
                for(int i = 0 ; i < nbiter ; i++)
                    fprintf(f,"%016llx %016llx\n", (unsigned long long) (S[i]>>64), (unsigned long long) S[i]);
                fprintf(f,"temps pour trouver la solution = %f\n", wtime() - t1 );
                fflush(f);
            }
            rotate(X,rot);
        }
        done++;;
    }
    fprintf(f,"temps total = %f\n", wtime() - t1);
    return(0);
}