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
static char doc[] = "prelieff - Relief-F";
#else
static char doc[] = "prelieff - Parallel Relief-F with MPI";
#endif

static char args_doc[] = "ARFF_FILE RANK_FILE";

static struct argp_option options[] = {
	{"algorithm",	'a',	"INT",	0,	"Algorithm version (0 = P (Default); 1 = G)"},
	{"class",		'c',	"NAME",	0,	"Class attribute name (Default: \"Class\")"},
	{"difference",	'd',	"INT",	0,	"Difference metric (0 = genotype (Default); 1 = allele-sharing)"},
	{"arff",		'r',	"FILE",	0,	"Output ARFF File (Default: none)"},
	{"prune",		'p',	"NUM",	0,	"Number (or percent) of attributes to prune (Default: 0)"},
	{ 0 }
};

struct arguments {
	char *args[2];
	int algorithm, difference;
	char *class;
	char *prune;
	char *arff_out;
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state) {
	struct arguments *arguments = state->input;

	switch (key) {
		case 'a':
			arguments->algorithm = atoi( arg );
			break;
		case 'c':
			arguments->class = arg;
			break;
		case 'd':
			arguments->difference = atoi( arg );
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
				argp_usage( state );

			arguments->args[state->arg_num] = arg;
			break;

		case ARGP_KEY_END:
			if (state->arg_num < 2)
				argp_usage( state );

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
	int *indexes;
	int i;
	FILE *outfile;

	struct arguments arguments;

	arguments.algorithm = 0;
	arguments.difference = 0;
	arguments.class = "Class";
	arguments.prune = "0";
	arguments.arff_out = "";

	if (argp_parse( &argp, argc, argv, 0, 0, &arguments))
		return 1;

#ifdef NO_MPI
	me = 0;
#else
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &me);
#endif

	outfile = fopen(arguments.args[1], "w");
	if (outfile == NULL) {
		fprintf (stderr, "Could not open file for writing: %s\n", arguments.args[1]);
		return 1;
	}

	info = read_arff (arguments.args[0], arguments.class);

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
	setVersion (arguments.algorithm);
	setDifference (arguments.difference);

	buildEvaluator (info, weights);
	if (me == 0) {
		int prune = 0;

		if(arguments.prune[strlen(arguments.prune + 1) - 1] == '%') {
			prune = (int)((atof(arguments.prune) * info->num_attributes) / 100);
		} else {
			prune = atoi(arguments.prune);
		}

		indexes =
			(int *) malloc_dbg (20, sizeof (int) * info->num_attributes);
		index_sort (indexes, weights, info->num_attributes);

		for (i = 0; i <= info->num_attributes - prune; i++) {
			fprintf (outfile, "%s,%.3f\n",
				info->attributes[indexes[i]]->name,
				evaluateAttribute (indexes[i]));
		}
		free (indexes);
	}

	release_read_info (info);
	free (weights);

#ifndef NO_MPI
	MPI_Finalize ();
#endif

	return 0;
}
