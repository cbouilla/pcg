/***** Version oneseq (incrément par defaut) *****/
#include <inttypes.h>
#include "pcg_setseq.h"

typedef __uint128_t pcg128_t;

#define PCG_128BIT_CONSTANT(high,low) \
            ((((pcg128_t)high) << 64) + low)
#define PCG_HAS_128BIT_OPS 1
            
#define PCG_DEFAULT_MULTIPLIER_128 \
        PCG_128BIT_CONSTANT(2549297995355413924ULL,4865540595714422341ULL)
#define PCG_DEFAULT_INCREMENT_128  \
        PCG_128BIT_CONSTANT(6364136223846793005ULL,1442695040888963407ULL)

// #define pcg64s_random_r                 pcg_oneseq_128_xsl_rr_64_random_r

struct pcg_state_128 {
    pcg128_t state;
};


void pcg_oneseq_128_srandom_r(struct pcg_state_128* rng, const pcg128_t initstate);

uint64_t pcg_oneseq_128_xsl_rr_64_random_r(struct pcg_state_128* rng);

/*****    Dans setseq.h    *****/
/*static inline uint64_t pcg_rotr_64(uint64_t value, unsigned int rot)
{
#if 0 && PCG_USE_INLINE_ASM && __clang__ && __x86_64__
    // For whatever reason, clang actually *does* generator rotq by
    // itself, so we don't need this code.
    asm ("rorq   %%cl, %0" : "=r" (value) : "0" (value), "c" (rot));
    return value;
#else
    return (value >> rot) | (value << ((- rot) & 63));
#endif
}



static inline uint64_t pcg_output_xsl_rr_128_64(pcg128_t state)
{
    return pcg_rotr_64(((uint64_t)(state >> 64u)) ^ (uint64_t)state, state >> 122u);
}*/

static inline void pcg_oneseq_128_step_r(struct pcg_state_128* rng)
{
    rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_128 + PCG_DEFAULT_INCREMENT_128;
}