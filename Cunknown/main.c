#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <time.h>

#include <mpi.h>

#include "fonctions.h"

static const bool VERBOSE = false;

enum chkpt_status {GOOD_CHECKPOINT, NO_CHECKPOINT, BAD_CHECKPOINT};
struct checkpoint_t {
    u64 known_bits;
    u64 X[nbiter];
};
u64 X[12];
enum task_status_t {READY, SENT, DONE};
int tag_hello = 0;
int tag_task = 1;
int tag_kill = 2;
int tag_result = 3;
const char * journal_filename = "journal.log";


enum chkpt_status load_journal(enum task_status_t *tasks)
{
	/* try to load checkpoint file */
	FILE *f = fopen(journal_filename, "r");
	if (f == NULL) {
		perror("Cannot open checkpoint file");
		return NO_CHECKPOINT;
	}
	struct checkpoint_t chkpt;
	size_t check = fread(&chkpt, sizeof(chkpt), 1, f);

	if (check != 1) {
		perror("Cannot read checkpoint from file");
		return NO_CHECKPOINT;
	}

	/* verify checkpoint */
	if (known_low != chkpt.known_bits) {
		printf("Guessed bits mismatch. Now=%d, in checkpoint=%lld.\n", known_low, chkpt.known_bits);
		return BAD_CHECKPOINT;
	}
	for (int i = 0; i < nbiter; i++)
		if (X[i] != chkpt.X[i]) {
			printf("X[%d] mismatch. Now=%llx, in checkpoint=%llx.\n", i, X[i], chkpt.X[i]);
		   	return BAD_CHECKPOINT;
		}

	/* header is fine */
	printf("Correct checkpoint loaded from %s.\n", journal_filename);

	/* load tasks already done */
	int found = 0;
	while (1) {
		int i;
		size_t x = fread(&i, sizeof(int), 1, f);
		if (x != 1) {
			if (feof(f))
				break;
			else
				err(1, "reading journal file");
		}
		tasks[i] = DONE;			
		found++;
	}
	printf("%s: %d tasks done\n", journal_filename, found);
	fclose(f);    

	return GOOD_CHECKPOINT;
}


void save_journal_header()
{
	struct checkpoint_t chkpt;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/* prepare checkpoint data */
	chkpt.known_bits = known_low;
	for (int i = 0; i < nbiter; i++)
		chkpt.X[i] = X[i];

	/* try to open journal file */
	FILE *f = fopen(journal_filename, "w");
	if (f == NULL) {
		perror("WARNING ! Cannot open journal file for writing");
		return;
	}
	size_t check = fwrite(&chkpt, sizeof(chkpt), 1, f);
	fclose(f);
	if (check != 1) {
		perror("WARNING ! Cannot write temporary checkpoint file");
		return;
	}
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
    double start = MPI_Wtime();
    u64 W0 = current >> (known_low - 1);
    u64 WC = 1 + 2 * (current % (1 << (known_low - 1)));
    
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
    double end = MPI_Wtime();

    printf("Done task %llx (W_0=%04llx / W_c=%04llx) [%.1fs]\n", current, W0, WC, end - start);
}


/* 
 * This runs an (embarassingly parallel) exhaustive search. 
 * The search space is split into "tasks". Tasks are numbered.
 * Tasks have the same running time
 * All tasks must be done to solve the problem.
 * Tasks are distributed to MPI ranks sequentially (block-wise).
 * Each MPI rank run tasks on all threads synchronously.
 * Inside an MPI rank, the state of each thread is initialized by init_task().
 * All threads run do_task() with a different task number synchronously.
 * After all threads have run a task, a checkpoint is taken.
 *
 * On jean-zay, tasks take between 41 and 45s / hyperthread. 
 * This mandates a master-slave design.
 *
 * With known_low=11, there are 2^21 tasks.
 * Let's conservatively assume that a task requires 44s.
 * Each jz node does 6545 tasks/hour ("heavy AVX512 mode" with maximum turbo).
 * Observed 6750 task/h on one node.
 * So we need 310 node-hours (we will be billed 12.4K cpu-hours).
 */
void master() 
{
	int n_tasks = 1 << (2 * known_low - 1);    
	enum task_status_t *tasks = malloc(sizeof(*tasks) * n_tasks);
	if (tasks == NULL)
		err(1, "cannot allocate task array!");

	/* All tasks must be done */
	for (int i = 0; i < n_tasks; i++)
		tasks[i] = READY;
	int task_ptr = 0;

	enum chkpt_status status = load_journal(tasks);
	switch (status) {
	case BAD_CHECKPOINT:
	    printf("BAD CHECKPOINT. Refusing to start. Please clean up the mess\n");
	    exit(EXIT_FAILURE);
	case NO_CHECKPOINT:
	    printf("COLD START.\n");
	    save_journal_header();
	    break;
	case GOOD_CHECKPOINT:
	    printf("WARM START.\n");
	    break;
	}

	/* DEBUG */
	// W_0 = 0x01ec
	// W_c = 0x040b
	// task_ptr = (0x01ec << (known_low - 1)) + (0x040b - 1) / 2;

	FILE *journal_file = fopen(journal_filename, "a");
	if (journal_file == NULL)
		err(1, "Impossible to open journal file %s for append\n", journal_filename);	
	double start = MPI_Wtime();
	
	int active_ranks;
	int done = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &active_ranks);
	active_ranks *= 2; /* because of the hack that runs two "threads" per MPI rank */

	while (active_ranks > 0) {
		MPI_Status status;
		int msg;
		MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
		if (status.MPI_TAG == tag_result) {
			if (VERBOSE)
				printf("RESULT(%d) <-- %d\n", msg, status.MPI_SOURCE);
			tasks[msg] = DONE;
			fwrite(&msg, sizeof(int), 1, journal_file);
			fflush(journal_file);
			done++;
		} else {
			if (VERBOSE)
				printf("HELLO <-- %d\n", status.MPI_SOURCE);
		}

		while (task_ptr < n_tasks && tasks[task_ptr] != READY)
			task_ptr++;
		if (task_ptr < n_tasks) {
			MPI_Bsend(&task_ptr, 1, MPI_INT, status.MPI_SOURCE, tag_task, MPI_COMM_WORLD);
			tasks[task_ptr] = SENT;
			if (VERBOSE)
				printf("TASK(%d) --> %d\n", task_ptr, status.MPI_SOURCE);
			task_ptr++;
		} else {
			MPI_Bsend(NULL, 0, MPI_INT, status.MPI_SOURCE, tag_kill, MPI_COMM_WORLD);
			active_ranks--;
			if (VERBOSE)
				printf("KILL --> %d\n", status.MPI_SOURCE);
		}

		if (done > 0) {
			double stop = MPI_Wtime();		
			double rate = done / (stop - start) * 3600;
			printf("%d tasks done in %.1fs -----> %.0f tasks / hour\n", done, stop - start, rate);
		}
	}
	fclose(journal_file);
	double stop = MPI_Wtime();
	printf("Finished in %.1fs\n", stop-start);
}

/* 
 * OK, so a bug in Intel's MPI library prevents us from using one MPI rank
 * per hyperthread, so here we cheat.
 */
void slave()
{
	int T = omp_get_max_threads();
	for (int i = 0; i < T; i++)
		MPI_Send(NULL, 0, MPI_INT, 0, tag_hello, MPI_COMM_WORLD);

	struct task_t task[T];
	#pragma omp parallel
	init_task(&task[omp_get_thread_num()]);

	bool stop = false;
	while (!stop) {
		int task_id[T];
		MPI_Status status;
		for (int i = 0; i < T; i++) {
			MPI_Recv(&task_id[i], 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (status.MPI_TAG == tag_kill) {
				task_id[i] = -1;
				stop = true;
			}
		}
		#pragma omp parallel
		{
			int t = omp_get_thread_num();
			if (task_id[t] >= 0)
				do_task(task_id[t], &task[t], X);
		}
		for (int i = 0; i < T; i++)
			MPI_Bsend(&task_id[i], 1, MPI_INT, 0, tag_result, MPI_COMM_WORLD);
	}
}


int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
	if (provided < MPI_THREAD_FUNNELED)
		errx(1, "MPI Thread support not sufficient");

	unsigned char Mpi_buffer[32768];
	MPI_Buffer_attach(Mpi_buffer, sizeof(Mpi_buffer));

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	init_var_globales();
	
	// X[ 0] = 0xe60f77ceac8f6cd9;  // W_0 = 00ed
	X[ 0] = 0xa21bd9d52bb064fa;  // W_0 = 01ec
	X[ 1] = 0xde3e78a8e09e5d29;  // W_0 = 00a7
	X[ 2] = 0xdbbe7755775e52ac;  // W_0 = 030e
	X[ 3] = 0xe553b06143b2e108;  // W_0 = 02d1
	X[ 4] = 0x4d99972ab762dae8;  // W_0 = 0460
	X[ 5] = 0x27af5ab26f846d1d;  // W_0 = 01eb
	X[ 6] = 0xd360f68fcc697b33;  // W_0 = 0262
	X[ 7] = 0x746b18d4b67a4e61;  // W_0 = 0475
	X[ 8] = 0x9be7410a83ddc14a;  // W_0 = 0594
	X[ 9] = 0x6c2b07bed5b9fe0b;  // W_0 = 04ef
	X[10] = 0x2fc807378a7bff2b;  // W_0 = 0276
    
	if (rank == 0)
    		master();
	else
    		slave();

	MPI_Finalize();
	return 0;
}