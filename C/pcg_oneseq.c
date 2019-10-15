#include "pcg_oneseq.h"

uint64_t pcg_rotr_64(uint64_t value, unsigned int rot)
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

uint64_t pcg_output_xsl_rr_128_64(pcg128_t state)
{
    return pcg_rotr_64(((uint64_t)(state >> 64u)) ^ (uint64_t)state,
                       state >> 122u);
}

void pcg_oneseq_128_step_r(struct pcg_state_128* rng)
{
    rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_128
                 + PCG_DEFAULT_INCREMENT_128;
}

void pcg_oneseq_128_srandom_r(struct pcg_state_128* rng,
                                     pcg128_t initstate)
{
    rng->state = 0U;
    pcg_oneseq_128_step_r(rng);
    rng->state += initstate;
    pcg_oneseq_128_step_r(rng);
}

uint64_t pcg_oneseq_128_xsl_rr_64_random_r(struct pcg_state_128* rng)
{
    pcg_oneseq_128_step_r(rng);
    return pcg_output_xsl_rr_128_64(rng->state);
}