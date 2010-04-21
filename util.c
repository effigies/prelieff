#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

#if 0

/* #ifdef NO_MPI */
void *malloc_dbg (int n, int x)
{
	printf ("%d Allocating %d bytes... ", n, x);
	void *tmp = malloc (x);
	if (tmp) {
		printf ("SUCCESS\n");
		fflush (stdout);
	} else {
		printf ("FAILURE\n");
		fflush (stdout);
		exit (1);
	}
	return tmp;
}
#else
void *malloc_dbg (int n, int x)
{
	return malloc (x);
}
#endif

int *remove_int (int *array, int len, int elem)
{
	int skip;
	int *ret = calloc (len - 1, sizeof (int));

	for (skip = 0; array[skip] != elem && skip < len; skip++);

	if (skip == len) {
		free (ret);
		return NULL;
	}

	memcpy (ret, array, skip * sizeof (int));
	memcpy (&ret[skip], &array[skip + 1],
		(len - skip - 1) * sizeof (int));
	return ret;
}
