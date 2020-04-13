#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <time.h>

#include "fonctions.h"

void result_found(u64 W0, u64 WC, u64 r)
{
    /* we print it just in case something goes wrong... */
    printf("solution found : W_0 = %04llx / W_c = %04llx / r = %08llx\n", W0, WC, r);
}


void do_task(u64 current, struct task_t *task, const u64 (*X)[nboutput])
{
    u64 W0 = current >> (known_low - 1);
    u64 WC = 1 + 2 * (current % (1 << (known_low - 1)));
	printf("Doing task %llx (W_0=%04llx / W_c=%04llx)\n", current, W0, WC);
    
    prepare_task(X, W0, WC, task);
    for (u64 r = 0; r < 1 << (nbiter * known_up); r++) {
	task->rot[0] = (task->rot[0] + 1) % k;
	int i = 0;
	while (task->rot[i] == 0 && i < nbiter) {
	    i++;
	    task->rot[i] = (task->rot[i] + 1) % k;
	}
	if (solve_isgood(task))
	    result_found(W0, WC, r);
    }
    finish_task(X, task);
}




int main()
{
	init_var_globales();
	
	
	u64 X[nboutput];
	
	int W_c = 0x035b;
	//X[ 0] = 0x2a33f0ac9ba281a3;	 int W_0 = 0x0271;
	X[ 0] = 0xd4166f4c3e02d10a;	 int W_0 = 0x01d0;
	X[ 1] = 0x1d1ceb21e7737101;	 // W_0 = 006b
	X[ 2] = 0xf8b90f473a5426d3;	 // W_0 = 0232
	X[ 3] = 0xe3a3b7babb2ad9ca;	 // W_0 = 06d5
	X[ 4] = 0x0077f2c80987dd13;	 // W_0 = 00c4
	X[ 5] = 0xf8ddaf2431548a13;	 // W_0 = 002f
	X[ 6] = 0x80935e041bbab85a;	 // W_0 = 0206
	X[ 7] = 0xbe0fde3939201c50;	 // W_0 = 02f9
	X[ 8] = 0xe9604fdf6b2177b7;	 // W_0 = 0678

	/* DEBUG */
	u64 task_id = (W_0 << (known_low - 1)) + (W_c - 1) / 2;

	struct task_t task;
	init_task(&task);
	do_task(task_id, &task, &X);

	return 0;
}