#include "pcg_oneseq.h"


void pcg_oneseq_128_srandom_r(struct pcg_state_128* rng, const pcg128_t initstate)
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