#include "index_sort.h"
#include <stdlib.h>

static double *m_sortArray;

int compindex (const void *p1, const void *p2)
{
	double x1 = m_sortArray[*(int *) p1];
	double x2 = m_sortArray[*(int *) p2];
	return (x1 < x2) ? 1 : ((x1 == x2) ? 0 : -1);
}

void index_sort (int *index, double *x, int size)
{
	int i;
	for (i = 0; i < size; i++) {
		index[i] = i;
	}

	m_sortArray = x;
	qsort (index, size, sizeof (int), compindex);
}
