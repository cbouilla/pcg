#include <assert.h>
#include <math.h>

#include "fonctions.h"

int testFonctions()
{
    assert(known_low == 11);

    u64 X[nboutput];
    pcg128_t vraiS[nboutput];

    pcg128_t S0 = (((pcg128_t) 5995207026785010249u) << k) + ((pcg128_t) 179350442155841024u);
    pcg128_t c = ((((pcg128_t) 6364136223846793005u) << k) + 1442695040888963407u) >> 1;
    pcg(vraiS, X, S0, &c, nboutput);
    
    if (vraiS[2] != (((pcg128_t) 1792771836637573954u) << k) + ((pcg128_t) 11139816115278170276u)) {
        printf("erreur sur pcg\n");
        return 0;
    }
    
    u64 W0 = (u64) (vraiS[0] % (1<<known_low));
    u64 WC = (u64) (c % (1<<known_low));

    struct task_t task;
    init_task(&task);
    prepare_task(X, W0, WC, &task);
    for (int i = 0 ; i < nbiter ; i++)
        task.rot[i] = (int) (vraiS[i] >> (2 * k - known_up));

    printf("1..7\n");
    printf("# nbiter = %d\n", nbiter);
    printf("# known_low = %d\n", known_low);

    /* test unrotate */
    u64 uX[nbiter];
    unrotateX(uX, X, task.rot);
    if (uX[0] != 15007519919903780682u) {
        printf("not ok 1 - erreur sur unrotateX\n");
    } else {
        printf("ok 1 - unrotate\n");
    }

    /* test getY */
    u64 Y[nbiter];
    printf("# W0 = %llu\n", W0);
    printf("# Wc = %llu\n", WC);

    getY(Y, W0, WC, task.rot, uX);
    if(Y[3] != 129714){
        printf("not ok 2 - erreur sur getY. Y[3] == %llu / attendu : 129714\n", Y[3]);
    } else {
        printf("ok 2 - getY\n");
    }

    /* test getYprim */
    u64 Yprim[nbiter];
    getYprim(Yprim, Y, W0, WC);
    if(Yprim[1] != 93486){
        printf("not ok 3 - erreur sur getY\n");
    }  else {
        printf("ok 3 - getYprim\n");
    }

    /* test getDY */
    u64 DY[nbiter-1];
    getDY(DY, Yprim);
    if(DY[0] != 14609){
        printf("not ok 4 - erreur sur getDY\n");
    }  else {
        printf("ok 4 - getDY\n");
    }

    /* test FindDS64 */
    u64 DS64[nbiter - 1];
    FindDS64(DS64, Y, uX, task.rot, task.lowSumPol, task.sumPolY);
    if (DS64[0] != 2055999906439120392u) {
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

    /*test getGoodY*/
    
    // u64 uXnbiter1 = unrotate(X[nbiter + 1], rot[nbiter + 1]);
    // u64 Ynbiter1 = ((((u64) ((polA[nbiter + 1] * WC + powA[nbiter + 1] * W0) % (1 << known_low))) ^ (uXnbiter1 % (1 << known_low))) << known_up) + (rot[nbiter + 1] ^ (uXnbiter1 >> (k - known_up)));
    // if (!checkY(goodY, 1, Ynbiter1)) {
    // //if(!goodY[Ynbiter1 + (1<<(known_low + known_up))]){
    //     printf("not ok 6 - erreur sur getGoodY\n");
    // }  else {
    //     printf("ok 6 - getGoodY\n");
    // }
    
    if (!solve_isgood(&task)) {
        printf("not ok 6 - erreur sur solve_isgood\n");
    }  else {
        printf("ok 6 - solve_isgood\n");
    }

    u64 DS640, Y0;
    solve(&task, &DS640, &Y0);
    printf("ok 7 - solve");

    return 1;
}



void printVal(pcg128_t S0, pcg128_t c) {
    int rot[nboutput];
    pcg128_t vraiS[nboutput];
    u64 X[nboutput];
    pcg(vraiS, X, S0, &c, nboutput);
    printf("setseq : %llu %llu\n", (u64) (vraiS[2]>>64), (u64) vraiS[2]);
    
    //printf("%llu %llu\n", (u64) ((vraiS[1] - vraiS[0])>>64), (u64) (vraiS[1] - vraiS[0]));
    
    //done
    for(int i = 0 ; i < nboutput ; i++){
        rot[i] = (int) (vraiS[i] >> (2 * k - known_up));
    }
    
    u64 W0 = (u64) (vraiS[0] % (1<<known_low));
    u64 WC = (u64) (c % (1<<known_low));
    printf("W0 : %llu\n", W0);
    printf("WC : %llu\n", WC);
    
    u64 uX[nbiter];
    unrotateX(uX, X, rot);
    printf("uX\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", uX[i]);
    printf("\n");
    /*FILE *f;
    f = fopen("result.txt","w");
    
    fprintf(f,"W0 : %llu\n", W0);*/
    u64 Y[nbiter];
    getY(Y, W0, WC, rot, uX);
    printf("Y :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", Y[i]);
    printf("\n");
    
    
    /***** Tests de vérification des sous-fonctions *****/
    u64 Yprim[nbiter];
    getYprim(Yprim, Y, W0, WC);
    printf("Yprim :\n");
    for(int i = 0 ; i < nbiter ; i++)
        printf("%llu ", Yprim[i]);
    printf("\n");
    
    u64 DY[nbiter];
    getDY(DY, Yprim);
    printf("DY :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%llu ", DY[i]);
    printf("\n");

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
    
    u64 DS64[nbiter - 1];
    FindDS64(DS64, Y, uX, rot, lowSumPol,sumPolY);
    printf("DS64 :\n");
    for(int i = 0 ; i < nbiter - 1 ; i++)
        printf("%llu ", DS64[i]);
    printf("\n");
    
    printf("testDS640 : %d\n", testDS640(DS64[0], X, Y[0], sumPolTest, lowSumPol)); 
    
}

void print128(pcg128_t x)
{
    printf("%016" PRIx64 " %016" PRIx64, (uint64_t) (x >> 64), (uint64_t) x);
}


int main()
{    
    /*  INITIALISATION DES PARAMETRES  */   
    init_var_globales();
    testFonctions();
    exit(0);
}
