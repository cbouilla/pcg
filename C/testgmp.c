#include <gmp.h>

int main(){
    mpz_t m;
    mpz_init_set_ui(m,2549297995355413924); 
    mpz_mul_2exp(m,m,64);
    mpz_add_ui(m,m,4865540595714422341);
    gmp_printf("%Zd\n",m);
    
    return 0;
}