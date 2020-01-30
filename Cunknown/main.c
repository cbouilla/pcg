#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <time.h>

#include <omp.h>
#include <mpi.h>

#include "fonctions.h"

static const bool VERBOSE = true;
static const char *CHKPT_FILENAME = "checkpoint.bin";
static const char *ALT_CHKPT_FILENAME = "checkpoint.bin.tmp";

enum chkpt_status {GOOD_CHECKPOINT, NO_CHECKPOINT, BAD_CHECKPOINT};

struct checkpoint_t {
    int nranks;
    int known_bits;
    u64 X[nbiter];
    u64 done;
    time_t when;
};


enum chkpt_status load_chkpt(int size, const u64 *X, u64 *done)
{
    *done = 0;

    /* try to load checkpoint file */
    FILE *f = fopen(CHKPT_FILENAME, "r");
    if (f == NULL) {
        perror("Cannot open checkpoint file");
        return NO_CHECKPOINT;
    }
    struct checkpoint_t chkpt;
    size_t check = fread(&chkpt, sizeof(chkpt), 1, f);
    fclose(f);
    if (check != 1) {
        perror("Cannot read checkpoint from file");
        return NO_CHECKPOINT;
    }

    /* verify checkpoint */
    if (size != chkpt.nranks) {
        printf("Communicator size mismatch. Now=%d, in checkpoint=%d.\n", size, chkpt.nranks);
        return BAD_CHECKPOINT;
    }
    if (known_low != chkpt.known_bits) {
        printf("Guessed bits mismatch. Now=%d, in checkpoint=%d.\n", known_low, chkpt.known_bits);
        return BAD_CHECKPOINT;
    }
    for (int i = 0; i < nbiter; i++)
        if (X[i] != chkpt.X[i]) {
           printf("X[%d] mismatch. Now=%llx, in checkpoint=%llx.\n", i, X[i], chkpt.X[i]);
            return BAD_CHECKPOINT;
        }

    /* checkpoint is fine */
    struct tm *tmp;
    tmp = localtime(&chkpt.when);
    if (tmp == NULL)
        err(1, "localtime");
    char outstr[255];
    if (strftime(outstr, sizeof(outstr), "%Y-%m-%d %H:%M:%S", tmp) == 0)
        errx(1, "strftime returned 0");
    printf("Correct checkpoint loaded from %s. Time = %s\n", CHKPT_FILENAME, outstr);
    printf("Tasks done per MPI rank: %lld.\n", chkpt.done);

    *done = chkpt.done;
    return GOOD_CHECKPOINT;
}


void save_chkpt(const u64 *X, u64 done)
{
    struct checkpoint_t chkpt;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* prepare checkpoint data */
    chkpt.nranks = size;
    chkpt.known_bits = known_low;
    for (int i = 0; i < nbiter; i++)
        chkpt.X[i] = X[i];
    chkpt.done = done;
    chkpt.when = time(NULL);

    /* try to open alternate checkpoint file */
    FILE *f = fopen(ALT_CHKPT_FILENAME, "w");
    if (f == NULL) {
        perror("WARNING ! Cannot open temporary checkpoint file");
        return;
    }
    size_t check = fwrite(&chkpt, sizeof(chkpt), 1, f);
    fclose(f);
    if (check != 1) {
        perror("WARNING ! Cannot write temporary checkpoint file");
        return;
    }

    /* writing the new checkpoint was successful: we erase an eventual old one. */
    if (rename(ALT_CHKPT_FILENAME, CHKPT_FILENAME) != 0)
        perror("WARNING ! Cannot rename tmp checkpoint file");
}



/* invoked at the beginning. Sets the range for the current MPI rank. */
void restart(const u64 *X, u64 *range_start, u64 *range_end, u64 *done)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    /* default */
    u64 n_tasks = 1 << (2 * known_low - 1);
    u64 tasks_per_rank = n_tasks / size;
    if (n_tasks % size != 0)
        tasks_per_rank += 1;

    *range_start = rank * tasks_per_rank;
    *range_end = (rank + 1) * tasks_per_rank;

    if (rank == 0) {
        enum chkpt_status status = load_chkpt(size, X, done);
        switch (status) {
        case BAD_CHECKPOINT:
            printf("BAD CHECKPOINT. Refusing to start. Please clean up the mess\n");
            exit(EXIT_FAILURE);
        case NO_CHECKPOINT:
            printf("COLD START.\n");
            break;
        case GOOD_CHECKPOINT:
            printf("WARM START.\n");
            break;
        }
    }

    MPI_Bcast(done, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    *range_start += *done;

    /* clip */
    if (*range_end > n_tasks)
        *range_end = n_tasks;
    

    if (VERBOSE)
        printf("MPI rank %d : [%llx:%llx]\n", rank, *range_start, *range_end);
}


/* checkpoints the current MPI rank */
void checkpoint(const u64 *X, u64 done)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* synchronize everybody */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        save_chkpt(X, done);
}


void result_found(const u64 *X, u64 W0, u64 WC, u64 r)
{
    /* we print it just in case something goes wrong... */
    printf("solution found : W_0 = %04llx / W_c = %04llx / r = %08llx\n", W0, WC, r);
    
    u64 c = (W0 << (known_low - 1)) + WC;
    char filename[255];
    sprintf(filename, "results/solution-%08llx.txt", c);
    FILE *f = fopen(filename, "a");
    if (f == NULL)
        err(1, "cannot open solution file");
    for (int i = 0; i < nbiter; i++)
        fprintf(f, "X[%d] = %llx\n", i, X[i]);
    fprintf(f, "W_0 = %04llx / W_c = %04llx / r = %08llx\n", W0, WC, r);
    fprintf(f, "==============================================\n");
    fclose(f);
}


void do_task(u64 current, struct task_t *task, const u64 *X)
{
    u64 W0 = current >> (known_low - 1);
    u64 WC = 1 + 2 * (current % (1 << (known_low - 1)));
    printf("Doing task %lld / %lld\n", W0, WC);
    
    prepare_task(X, W0, WC, task);
    for (u64 r = 0; r < 1 << (nbiter * known_up); r++) {
        task->rot[0] = (task->rot[0] + 1) % k;
        int i = 0;
        while (task->rot[i] == 0 && i < nbiter) {
            i++;
            task->rot[i] = (task->rot[i] + 1) % k;
        }
        if (solve_isgood(task))
            result_found(X, W0, WC, r);
    }
    finish_task(X, task);
}


int main(int argc, char **argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
        errx(1, "MPI Thread support not sufficient");

    /*  INITIALISATION DES PARAMETRES  */
    init_var_globales();
        
    u64 X[9];
    X[0] = 0x47a42ee112e8afb9;
    X[1] = 0xf5e7948dbc0c7e26;
    X[2] = 0x91724bdca45a78a4;
    X[3] = 0x1be0e7e5b398b248;
    X[4] = 0x6f8b727451e185a8;
    X[5] = 0x976d59bba78ef4e2;
    X[6] = 0xc588c4c6c9052cba;
    X[7] = 0x9cc0fc58615e1b87;
    X[8] = 0xec7c5d6ee9992147;

    u64 range_start, range_end, done;
    restart(X, &range_start, &range_end, &done);

    double t1 = wtime();

    /* init all tasks */
    int T = omp_get_max_threads();
    struct task_t task[T];
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        init_task(&task[tid]);
    }
    
    /* DEBUG */
    // W_0 = 7984
    // W_c = 7673
    range_start = (7984 << (known_low - 1)) + (7673 - 1) / 2;

    while (range_start < range_end) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            if (range_start + tid < range_end) {
                do_task(range_start + tid, &task[tid], X);
            }
        }
        range_start += T;
        done += T;
        checkpoint(X, done);
    }

    printf("temps total = %f\n", wtime() - t1);
    MPI_Finalize();
    return(0);
}
