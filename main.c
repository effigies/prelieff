#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "arff.h"
#include "prelieff.h"
#include "java.h"
#include "index_sort.h"
#include "util.h"
#ifndef NO_MPI
#include "mpi.h"
#endif

const char *argp_program_version = "prelieff 0.2";
const char *argp_program_bug_address = "<chris-johnson@utulsa.edu>";

#ifdef NO_MPI
static char doc[] = "prelieff - Relief-F (Compiled without MPI)";
#else
static char doc[] = "prelieff - Parallel Relief-F with MPI";
#endif

static char args_doc[] = "ARFF_FILE RANK_FILE";

static struct argp_option options[] = {
	{"algorithm", 'a', "INT", 0,
	 "Algorithm version (0 = P (Default); 1 = G)"},
	{"class", 'c', "NAME", 0,
	 "Class attribute name (Default: \"Class\")"},
	{"difference", 'd', "INT", 0,
	 "Difference metric (0 = genotype (Default); 1 = allele-sharing)"},
	{"arff", 'r', "FILE", 0, "Output ARFF File (Default: none)"},
	{"prune", 'p', "NUM", 0,
	 "Number (or percent) of attributes to prune (Default: 0)"},
	{0}
};

struct arguments {
	char *args[2];
	int algorithm, difference;
	char *class;
	char *prune;
	char *arff_out;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key) {
	case 'a':
		arguments->algorithm = atoi (arg);
		break;
	case 'c':
		arguments->class = arg;
		break;
	case 'd':
		arguments->difference = atoi (arg);
		break;
	case 'r':
		arguments->arff_out = arg;
		break;
	case 'p':
		/* Maintain string form until we know how many attributes there are, in case this is a percentage. */
		arguments->prune = arg;
		break;

	case ARGP_KEY_ARG:
		if (state->arg_num >= 2)
			argp_usage (state);

		arguments->args[state->arg_num] = arg;
		break;

	case ARGP_KEY_END:
		if (state->arg_num < 2)
			argp_usage (state);

		break;

	default:
		return ARGP_ERR_UNKNOWN;
	}

	return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

int main (int argc, char **argv)
{
	arff_info_t *info;
	int me;
	double *weights;
	int *indices;
	int i;
	FILE *outfile;
	FILE *arfffile;
	int prune = 0;

	/* Argument parsing */
	struct arguments arguments;

	arguments.algorithm = 0;	// Default to the P algorithm
	arguments.difference = 0;	// Default to genotype difference metric
	arguments.class = "Class";	// Default class name
	arguments.prune = "0";	// Prune 0 attributes by default
	arguments.arff_out = NULL;	// Do not write a new ARFF file by default

	if (argp_parse (&argp, argc, argv, 0, 0, &arguments))
		return 1;
	/* End argument parsing */


#ifdef NO_MPI
	me = 0;
#else
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &me);
#endif

	/* Verify we can write to our output files before doing work */
	outfile = fopen (arguments.args[1], "w");
	if (outfile == NULL) {
		fprintf (stderr, "Could not open file for writing: %s\n",
			 arguments.args[1]);
		return 1;
	}

	/* Only check on the ARFF file if we're writing it */
	if (arguments.arff_out != NULL) {
		arfffile = fopen (arguments.arff_out, "w");
		if (outfile == NULL) {
			fprintf (stderr,
				 "Could not open file for writing: %s\n",
				 arguments.args[1]);
			return 1;
		}
	}

	info = read_arff (arguments.args[0], arguments.class);

	if (info == NULL) {
		fprintf (stderr, "%s, line %i\n", get_last_error (),
			 get_lineno ());
		return 1;
	}

	/* Now that we know how many attributes there are, we can determine how
	 * many to prune, if we're pruning a percentage.
	 */
	if (arguments.prune[strlen (arguments.prune) - 1] == '%') {
		prune = (int) ((atof (arguments.prune) *
				info->num_attributes) / 100);
	} else {
		prune = atoi (arguments.prune);
	}

	/* Don't let's be silly. */
	if (prune >= info->num_attributes) {
		fprintf (stderr,
			 "Attempting to prune entire file. Not bothering to run Relief-F.\n");
		return 1;
	}

	weights = calloc (info->num_attributes, sizeof (double));

	resetOptions ();
	setSampleSize (-1);
	setNumNeighbours (10);
	setWeightByDistance (true);
	setSigma (2);
	setVersion (arguments.algorithm);
	setDifference (arguments.difference);

	buildEvaluator (info, weights);
	if (me == 0) {
		/* The number of attributes retained includes neither the class
		 * attribute nor the pruned attributes.
		 */
		int retained = info->num_attributes - 1 - prune;

		/* Rank the attributes, but remove the class attribute. */
		int *tmp = calloc (info->num_attributes, sizeof (int));
		index_sort (tmp, weights, info->num_attributes);
		indices =
			remove_int (tmp, info->num_attributes,
				    info->class_index);

		/* In case we've gone crazy, fall back to original behavior */
		if (indices == NULL)
			indices = tmp;
		else
			free (tmp);

		/* We generate a two-columned CSV */
		for (i = 0; i < retained; i++) {
			fprintf (outfile, "%s,%.3f\n",
				 info->attributes[indices[i]]->name,
				 evaluateAttribute (indices[i]));
		}

		/* Automatically generate an ARFF file, if requested.
		 * This is primarily useful if we are pruning, for instance, for
		 * iterated Relief-F, or filtering by another ARFF-accepting program.
		 */
		if (arfffile != NULL) {
			arff_info_t output;

			/* For the most part, we keep the same structures */
			memcpy (&output, info, sizeof (arff_info_t));

			/* The +1 is for the class attribute, which we ignored, earlier. */
			output.num_attributes = retained + 1;
			output.attributes =
				calloc (retained + 1,
							 sizeof
							 (attr_info_t));

			if (output.attributes == NULL) {
				fprintf (stderr,
					 "Could not allocate %d bytes for output.attributes.",
					 (retained +
					  1) * (sizeof (attr_info_t)));

				return 1;
			}

			/* We want to copy the retained attributes.
			 *
			 * Since these are more or less immutable, we don't need to
			 * deeply copy these data.
			 */
			for (i = 0; i < retained; i++)
				output.attributes[i] =
					info->attributes[indices[i]];

			output.attributes[retained] =
				info->attributes[info->class_index];

			output.instances =
				calloc (info->num_instances,
							sizeof (instance_t
								*));

			if (output.instances == NULL) {
				fprintf (stderr,
					 "Could not allocate %d bytes for output.instances.",
					 (info->num_instances) *
					 (sizeof (instance_t *)));

				return 1;
			}

			/* Populate the table, sorting attributes by rank */
			for (i = 0; i < info->num_instances; i++) {
				int j;

				output.instances[i] =
					calloc (1,
							       sizeof
							       (instance_t));

				if (output.instances[i] == NULL) {
					fprintf (stderr,
						 "Could not allocate %d bytes for output.instances[%d].",
						 (sizeof (instance_t)), i);

					return 1;
				}

				output.instances[i]->data =
					calloc (retained + 1,
							   sizeof (data_t));

				if (output.instances[i] == NULL) {
					fprintf (stderr,
						 "Could not allocate %d bytes for output.instances[%d].data.",
						 (retained +
						  1) * (sizeof (instance_t)),
						 i);

					return 1;
				}


				for (j = 0; j < retained; j++)
					output.instances[i]->data[j].ival =
						info->instances[i]->
						data[indices[j]].ival;

				output.instances[i]->data[retained].ival =
					info->instances[i]->data[info->
								 class_index].
					ival;
			}

			write_arff (&output, arfffile);

			/* Cleanup */
			free (output.attributes);
			for (i = 0; i < info->num_instances; i++) {
				free (output.instances[i]->data);
				free (output.instances[i]);
			}
			free (output.instances);
		}

		free (indices);
	}

	release_read_info (info);
	free (weights);

#ifndef NO_MPI
	MPI_Finalize ();
#endif

	return 0;
}
