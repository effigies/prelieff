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
	arff_info_t output;
	int me, selected, version;
	double *weights;
	int *indexes;
	int i, j, k;

	if (argc < 5) {
		printf ("Usage: prelieff <arff file> <class name> <number of attributes to select> <algorithm version>\n");
		return 1;
	}
#ifdef NO_MPI
	me = 0;
#else
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &me);
#endif

	info = read_arff (argv[1], argv[2]);
	selected = atoi (argv[3]);
	version = atoi (argv[4]);

	if (info == NULL) {
		printf ("%s, line %i\n", get_last_error (), get_lineno ());
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

		printf ("%% arff file autogenerated by prelieff\n\n");
		printf ("%% selected attrubutes:\n");
		printf ("%%   %-30s%-10s\n", "name", "weight");
		for (j = 0; j < selected; j++) {
			i = info->num_attributes - j - 1;
			printf ("%%   %-30s%-10.3f\n",
				info->attributes[indexes[i]]->name,
				evaluateAttribute (indexes[i]));
		}
		printf ("\n");

		memcpy (&output, info, sizeof (arff_info_t));
		output.num_attributes = selected + 1;
		output.attributes =
			(attr_info_t **) malloc_dbg (21, sizeof (attr_info_t) *
						 (selected + 1));
		for (j = 0; j < selected; j++) {
			i = info->num_attributes - j - 1;
			output.attributes[j] = info->attributes[indexes[i]];
		}
		output.attributes[selected] =
			info->attributes[info->class_index];
		output.instances =
			(instance_t **) malloc_dbg (22, sizeof (attr_info_t) *
						info->num_instances);
		for (j = 0; j < info->num_instances; j++) {
			output.instances[j] =
				(instance_t *) malloc_dbg (23, sizeof (instance_t));
			output.instances[j]->data =
				(data_t *) malloc_dbg (24, sizeof (data_t) *
						   (selected + 1));
			for (k = 0; k < selected; k++) {
				i = info->num_attributes - k - 1;
				output.instances[j]->data[k].ival =
					info->instances[j]->data[indexes[i]].
					ival;
			}
			output.instances[j]->data[selected].ival =
				info->instances[j]->data[info->class_index].
				ival;
		}
		write_arff (&output, stdout);
		free (output.attributes);
		for (j = 0; j < info->num_instances; j++) {
			free (output.instances[j]->data);
			free (output.instances[j]);
		}
		free (output.instances);
		free (indexes);
	}

	release_read_info (info);
	free (weights);

#ifndef NO_MPI
	MPI_Finalize ();
#endif

	return 0;
}