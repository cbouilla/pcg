#include "fonctions.h"

////////////////Fonctions pour la récupération de S//////////////


/* Y = S[k-known_up:k+known_low] */
void getY(unsigned long long *Y, unsigned long long W0, unsigned long long WC, int* rot, unsigned long long* uX){
    for(int i = 0 ; i < nbiter ; i++){
        Y[i] = ((((unsigned long long) ((polA[i] * WC + powA[i] * W0) % (1 << known_low))) ^ (uX[i] % (1 << known_low))) << known_up) + (rot[i] ^ (uX[i] >> (k - known_up)));
    }
}

/* Y' = Y - composante en WC et en W0 */
void getYprim(unsigned long long *Yprim, unsigned long long *Y, unsigned long long W0, unsigned long long WC){
    for(int i = 0 ; i < nbiter ; i++)
        Yprim[i] = (unsigned long long) (Y[i] - ((polA[i] * WC + powA[i] * W0) >> (k - known_up))) % (1<<(known_up + known_low));
}

/* DY = différence sur les Y' */
void getDY(unsigned long long *DY, unsigned long long* Yprim){
    for(int i = 0 ; i < nbiter - 1 ; i++)
        DY[i] = (Yprim[i+1] - Yprim[i]) % (1<<(known_low + known_up));
}

/* DS64 = différence sur S'[known_low:known_low+k], avec S' = S - composante en WC, W0 */
void FindDS64(unsigned long long* DS64, unsigned long long* Y0, unsigned long long* uX,int* rot, unsigned long long* lowSumPol, unsigned long long* sumPolY){
    unsigned long long tmp[nbiter];
    for(int i = 0 ; i < nbiter ; i++){//Y
        tmp[i] = (((lowSumPol[i] % (1 << known_low)) ^ (uX[i] % (1 << known_low))) << known_up) + (rot[i] ^ (uX[i] >> (k - known_up)));
    }
    Y0[0] = tmp[0];
    for(int i = 0 ; i < nbiter ; i++)//Yprim
        tmp[i] = (tmp[i] - sumPolY[i]) % (1<<(known_up + known_low));

    for(int i = 0 ; i < nbiter - 1 ; i++){ //DY
        tmp[i] = (tmp[i+1] - tmp[i]) % (1<<(known_low + known_up));
        tmp[i] = tmp[i] << (k - known_up - known_low);
    }

    float u[nbiter-1];
    prodMatVecFFU(u, invG, tmp, nbiter-1);
    for(int i = 0 ; i < nbiter-1 ; i++)
        tmp[i] = (unsigned long long) llroundf(u[i]);

    *DS64 = 0;
    for(int i = 0 ; i < nbiter-1 ; i++)
        (*DS64) += Greduite[i] * tmp[i];

    //prodMatVecU(DS64, Greduite, tmp, nbiter-1);
}


/* vérifie si DS64 est cohérent avec les X_i */
int testDS640(unsigned long long DS640,  unsigned long long* X, unsigned long long Y0, unsigned long long* sumPolTest, unsigned long long* lowSumPol){
    for(int i = nbiter ; i < nbtest + nbiter ; i++){
        unsigned long long Xi = X[i];
        unsigned long long tmp = polA[i] * DS640; //ATTENTION cast pcg128_t
        tmp += sumPolTest[i - nbiter];//ATTENTION cast pcg128_t
        unsigned long long Yi1 = (Y0 + (tmp >> (k - known_low - known_up))) % (1<< (known_low + known_up)); //avec ou sans retenue OK!
        unsigned long long Yi2 = Yi1 + 1;
        unsigned long long Wi = lowSumPol[i] % (1<< known_low);//ATTENTION cast pcg128_t
        int test1, test2,test = 0;
        for(int j = 0 ; j < k ; j++){
            test1 = (((Xi ^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) && ((j ^ (Xi >> (k - known_up))) == Yi1 % (1 << known_up));
            test2 = (((Xi ^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) && ((j ^ (Xi >> (k - known_up))) == Yi2 % (1 << known_up));
            if (test1 || test2){
                test = 1;
            }
            Xi = unrotate1(Xi);
        }
        if(!test) {
            //printf("erreur a i = %d\n", i);
            return 0;
        }
    }
    return 1;
    
}



void getDS64(pcg128_t* DS64, pcg128_t* S,unsigned long long W0,unsigned long long WC){
    for(int i = 0 ; i < nboutput - 1 ; i++)  
        DS64[i] = (S[i+1] - S[i] - powA[i] * WC - (powA[i+1] - powA[i]) * W0) % (((pcg128_t) 1) << (k + known_low));
}

void getDSmod(pcg128_t* DSmod, unsigned long long* DS64, unsigned long long W0, unsigned long long WC){
    for(int i = 0 ; i < nbiter - 1 ; i++)  
        DSmod[i] = ((((pcg128_t) DS64[i]) << known_low) + powA[i] * WC + (powA[i+1] - powA[i]) * W0) % (((pcg128_t) 1) << (k + known_low));
}

int FindRoti(int *roti, unsigned long long DS640, unsigned long long Xi,int  i,unsigned long long Y0,unsigned long long W0,unsigned long long WC){//OK !
    unsigned long long tmp = polA[i] * DS640;
    tmp += W0 * ((unsigned long long) (powA[i] >> known_low) - 1) + WC * ((unsigned long long) (polA[i] >> known_low) - 1);
    unsigned long long Yi1 = (Y0 + (tmp >> (k - known_low - known_up))) % (1<< (known_low + known_up)); //avec ou sans retenue OK!
    unsigned long long Yi2 = Yi1 + 1;
    unsigned long long Wi = (W0 * ((unsigned long long) powA[i]) + WC * ((unsigned long long) polA[i])) % (1<< known_low);
    int test1, test2,cptRot = 0;
    for(int i = 0 ; i < k ; i++){
        test1 = (((Xi ^ (Yi1 >> known_up)) % (1 << known_low)) == Wi) && ((i ^ (Xi >> (k - known_up))) == Yi1 % (1 << known_up));
        test2 = (((Xi ^ (Yi2 >> known_up)) % (1 << known_low)) == Wi) && ((i ^ (Xi >> (k - known_up))) == Yi2 % (1 << known_up));
        if (test1 || test2){
            roti[cptRot] = i;
            cptRot++;
        }
        Xi = unrotate1(Xi);
    }
    return cptRot;
}

int FindRot(int *rot, int* nbrot, unsigned long long DS640,unsigned long long* X,unsigned long long Y0,unsigned long long W0,unsigned long long WC,int n){
    for(int i = 0 ; i < n ; i++){
        nbrot[i] = FindRoti(rot + i * k, DS640, X[i], i, Y0, W0, WC);
        if(nbrot[i] == 0)
            return 0;
    }
    return 1;
}

void FindRotDS(int* rotDS, int* rot, unsigned long long DS640, unsigned long long W0, unsigned long long WC){//cohérent
    int tmp;
    for(int i = 0 ; i < nboutput-1 ; i++){
        tmp = rot[i+1] - rot[i];
        rotDS[i] = tmp - ((int) ((((powA[i] * DS640) << known_low) + WC * powA[i] + W0 * (powA[i+1] - powA[i])) >> (2 * k - known_up)));
        //printf("rotDS[%d] = %d\n", i, rotDS[i]);
    }
}

/*void FindUpDS(unsigned long long* upDS, int* rotDS){
    unsigned long long tmp[nboutput - 1];
    for(int i = 0 ; i < nbiter ; i++)
        tmp[i] = ((unsigned long long) rotDS[i]) << (k - known_low - known_up);
    float tmp2[nboutput - 1];
    prodMatVecFFU(tmp2, invG2, tmp, nboutput - 1);
    for(int i = 0 ; i < nbiter ; i++)
        tmp[i] = (unsigned long long) llroundf(tmp2[i]);
    prodMatVecU(upDS, Greduite2, tmp, nboutput - 1);
}

void FindDS(pcg128_t* DS, unsigned long long* upDS, unsigned long long DS640, unsigned long long W0,unsigned long long WC){
    for(int i = 0 ; i < nboutput - 1 ; i++)
        DS[i] = (((pcg128_t) upDS[i]) << (k + known_low)) + ((powA[i] * DS640) << known_low) + WC * powA[i] + W0 * (powA[i+1] - powA[i]);
}*/


/*void pcgone(pcg128_t *S, unsigned long long* X, pcg128_t S0, int n) //a besoin de pcg_oneseq
{
    struct pcg_state_128 rng;
    pcg_oneseq_128_srandom_r(&rng, S0);
    for (int i = 0 ; i < n ; i++) {
        X[i] = pcg_output_xsl_rr_128_64(rng.state);
        S[i] = rng.state;
        pcg_oneseq_128_step_r(&rng);
    }
}*/  