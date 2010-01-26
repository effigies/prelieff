#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "arff.h"
#include "prelieff.h"
#include "java.h"
#include "index_sort.h"
#include "util.h"
#ifndef NO_MPI
#include "mpi.h"
#endif

int main (int argc, char **argv)
{
	arff_info_t *info;
	int me, selected, version;
	double *weights;
	int *indexes;
	int i, j;
	FILE *outfile;

	if (argc < 6) {
		printf ("Usage: prelieff <arff file> <output file> <class name> <number of attributes to select> <algorithm version>\n");
		return 1;
	}
#ifdef NO_MPI
	me = 0;
#else
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &me);
#endif

	outfile = fopen(argv[2], "w");
	if (outfile == NULL) {
		fprintf (stderr, "Could not open file for writing: %s\n", argv[2]);
		return 1;
	}

	info = read_arff (argv[1], argv[3]);
	selected = atoi (argv[4]);
	version = atoi (argv[5]);

	if (info == NULL) {
		fprintf (stderr, "%s, line %i\n", get_last_error (), get_lineno ());
		return 1;
	}

	weights = (double *) malloc_dbg (19, sizeof (double) * info->num_attributes);

	resetOptions ();
	setSampleSize (-1);
	setNumNeighbours (10);
	setWeightByDistance (true);
	setSigma (2);
	setVersion (version);

	buildEvaluator (info, weights);
	if (me == 0) {
		indexes =
			(int *) malloc_dbg (20, sizeof (int) * info->num_attributes);
		index_sort (indexes, weights, info->num_attributes);

		for (j = 0; j < selected; j++) {
			i = info->num_attributes - j - 1;
			fprintf (outfile, "%-30s%-10.3f\n",
				info->attributes[indexes[i]]->name,
				evaluateAttribute (indexes[i]));
		}
		fprintf (outfile, "\n");
		free (indexes);
	}

	release_read_info (info);
	free (weights);

#ifndef NO_MPI
	MPI_Finalize ();
#endif

	return 0;
}
