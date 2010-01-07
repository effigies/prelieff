#include <stdio.h>
#include <stdlib.h>

#ifdef NO_MPI
void *malloc_dbg(int n, int x) {
	printf("%d Allocating %d bytes... ", n, x);
	void *tmp = malloc(x);
	if (tmp) {
		printf("SUCCESS\n");
	} else {
		printf("FAILURE\n");
		exit(1);
	}
	return tmp;
}
#else
void *malloc_dbg(int n, int x) { return malloc(x); }
#endif
