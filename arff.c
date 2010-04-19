#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "arff.h"
#include "util.h"

#define SUCCESS(x)     ((x) == 0)
#define FAIL(x)        ((x) != 0)
#ifdef TRY
#undef TRY
#endif
#define TRY(r,x)       { if (FAIL(r = (x))) return r; }
#define LOWERCASE(x)   ((((x) >= 'A') && ((x) <='Z')) ? ((x) - 'A' + 'a') : (x))
#define MIN(x,y)       (((x) < (y)) ? (x) : (y))

//#define ARFF_DEBUG_MAIN
//#define ARFF_DEBUG_MSG

#ifdef ARFF_DEBUG_MSG
#define DEBUGMSG(x)    printf x
#else
#define DEBUGMSG(x)
#endif

/* structures, enums, typedefs */
typedef enum {
	TOKEN_RELATION,
	TOKEN_ATTRIBUTE,
	TOKEN_DATA,
	TOKEN_NUMERIC,
	TOKEN_STRING,
	TOKEN_NEWLINE,
	TOKEN_COMMA,
	TOKEN_LEFTBRACE,
	TOKEN_RIGHTBRACE,
	TOKEN_TOKEN,
	TOKEN_NUM_TYPES,
	TOKEN_ANY
} token_t;

typedef struct {
	char *string;
	token_t type;
} special_token_t;

typedef enum {
	PARSE_STATE_NEWLINE = 0,
	PARSE_STATE_GET_NEWLINE,
	PARSE_STATE_RELATION1,
	PARSE_STATE_ATTRIBUTE1,
	PARSE_STATE_ATTRIBUTE2,
	PARSE_STATE_ATTRIBUTE_NOM1,
	PARSE_STATE_ATTRIBUTE_NOM2,
	PARSE_STATE_DATA1,
	PARSE_STATE_DATA2,
	PARSE_STATE_DATA3,
	PARSE_STATE_NUM,
	PARSE_STATE_ERROR
} parse_state_t;

typedef parse_state_t (*parse_action_t) (token_t, char *, arff_info_t *);

typedef struct {
	token_t token;
	parse_action_t action;
	parse_state_t nextstate;
} parse_table_t;

/* global data */
static special_token_t special[] = {
	{"@relation", TOKEN_RELATION},
	{"@attribute", TOKEN_ATTRIBUTE},
	{"@data", TOKEN_DATA},
	{"numeric", TOKEN_NUMERIC},
	{"real", TOKEN_NUMERIC},
	{NULL, 0}
};

static parse_state_t state;
static attr_info_t *first_attr;
static instance_t *first_inst;
static attr_info_t *curr_attr;
static instance_t *curr_instance;
static int curr_data;
static char *class_name;

static char error_string[200];
static int line_no;

static char *token_names[] = {
	"RELATION",
	"ATTRIBUTE",
	"DATA",
	"NUMERIC",
	"STRING",
	"newline",
	",",
	"{",
	"}",
	"TOKEN"
};

#ifdef ARFF_DEBUG_MSG
static char *state_names[] = {
	"PARSE_STATE_NEWLINE",
	"PARSE_STATE_GET_NEWLINE",
	"PARSE_STATE_RELATION1",
	"PARSE_STATE_ATTRIBUTE1",
	"PARSE_STATE_ATTRIBUTE2",
	"PARSE_STATE_ATTRIBUTE_NOM1",
	"PARSE_STATE_ATTRIBUTE_NOM2",
	"PARSE_STATE_DATA1",
	"PARSE_STATE_DATA2",
	"PARSE_STATE_DATA3",
	"PARSE_STATE_NUM",
	"PARSE_STATE_ERROR"
};
#endif

/* local function prototypes */
static arff_info_t *init (char *class_attribute_name);
static int parse (token_t type, char *tok, arff_info_t * info);
static void parse_end (arff_info_t * info);
static int lex (char *buf, arff_info_t * info);
static char *read_token (char **pbuf, token_t * delimiter);
static int stricmp (const char *x, const char *y);
static attr_info_t *add_attribute (arff_info_t * info, char *name);
static void add_nominal_class (nom_info_t * info, char *name);
static instance_t *add_instance (arff_info_t * info);

/* parser tables */
parse_state_t parse_relation1 (token_t token, char *name, arff_info_t * info);
parse_state_t parse_attribute1 (token_t token, char *name,
				arff_info_t * info);
parse_state_t parse_attribute2_numeric (token_t token, char *name,
					arff_info_t * info);
parse_state_t parse_attribute2_nominal (token_t token, char *name,
					arff_info_t * info);
parse_state_t parse_attribute_nom1 (token_t token, char *name,
				    arff_info_t * info);
parse_state_t parse_data (token_t token, char *name, arff_info_t * info);
parse_state_t parse_data_end (token_t token, char *name, arff_info_t * info);

static parse_table_t parse_newline_tab[] = {
	{TOKEN_RELATION, NULL, PARSE_STATE_RELATION1},
	{TOKEN_ATTRIBUTE, NULL, PARSE_STATE_ATTRIBUTE1},
	{TOKEN_DATA, NULL, PARSE_STATE_DATA1},
	{TOKEN_NEWLINE, NULL, PARSE_STATE_NEWLINE},
	{TOKEN_NUM_TYPES, NULL, 0}
};

static parse_table_t parse_get_newline_tab[] = {
	{TOKEN_NEWLINE, NULL, PARSE_STATE_NEWLINE},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_relation1_tab[] = {
	{TOKEN_ANY, parse_relation1, 0},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_attribute1_tab[] = {
	{TOKEN_ANY, parse_attribute1},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_attribute2_tab[] = {
	{TOKEN_NUMERIC, parse_attribute2_numeric},
	{TOKEN_LEFTBRACE, parse_attribute2_nominal},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_attribute_nom1_tab[] = {
	{TOKEN_TOKEN, parse_attribute_nom1},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_attribute_nom2_tab[] = {
	{TOKEN_COMMA, NULL, PARSE_STATE_ATTRIBUTE_NOM1},
	{TOKEN_RIGHTBRACE, NULL, PARSE_STATE_GET_NEWLINE},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_data1_tab[] = {
	{TOKEN_TOKEN, parse_data},
	{TOKEN_NEWLINE, NULL, PARSE_STATE_DATA1},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_data2_tab[] = {
	{TOKEN_COMMA, NULL, PARSE_STATE_DATA3},
	{TOKEN_NEWLINE, parse_data_end},
	{TOKEN_NUM_TYPES}
};

static parse_table_t parse_data3_tab[] = {
	{TOKEN_TOKEN, parse_data},
	{TOKEN_NUM_TYPES}
};

static parse_table_t *parse_table[PARSE_STATE_NUM] = {
	parse_newline_tab,
	parse_get_newline_tab,
	parse_relation1_tab,
	parse_attribute1_tab,
	parse_attribute2_tab,
	parse_attribute_nom1_tab,
	parse_attribute_nom2_tab,
	parse_data1_tab,
	parse_data2_tab,
	parse_data3_tab,
};

/* functions */
arff_info_t *read_arff (char *filename, char *class_attribute_name)
{
	arff_info_t *info = init (class_attribute_name);
	FILE *in;
	int size;
	char *buf;

	if (!(in = fopen (filename, "r"))) {
		sprintf (error_string, "unable to open file: %s", filename);
	} else {
		fseek (in, 0, SEEK_END);
		size = ftell (in);
		buf = (char *) malloc_dbg (25, size + 1);
		fseek (in, 0, SEEK_SET);
		size = fread (buf, 1, size, in);

		if (FAIL (lex (buf, info))) {
			release_read_info (info);
			info = NULL;
		}

		free (buf);
	}
	return info;
}

void write_arff (arff_info_t * info, FILE * out)
{
	char buf[1024];
	int i, j, k;
	nom_val_t *nom_val;
	double val;

	sprintf (buf, "@RELATION %s\n\n", info->relation_name);
	fwrite (buf, strlen (buf), sizeof (char), out);

	for (i = 0; i < info->num_attributes; i++) {
		sprintf (buf, "@ATTRIBUTE %s ", info->attributes[i]->name);
		fwrite (buf, strlen (buf), sizeof (char), out);

		switch (info->attributes[i]->type) {
		case ATTR_NUMERIC:
			fwrite ("REAL\n", 5, sizeof (char), out);
			break;
		case ATTR_NOMINAL:
			fwrite ("{", 1, sizeof (char), out);
			for (nom_val = info->attributes[i]->nom_info->first;
			     nom_val != NULL; nom_val = nom_val->next) {
				fwrite (nom_val->name, strlen (nom_val->name),
					sizeof (char), out);
				if (nom_val->next != NULL)
					fwrite (",", 1, sizeof (char), out);
			}
			fwrite ("}\n", 2, sizeof (char), out);
			break;
		default:
			break;
		}
	}

	fwrite ("\n@DATA\n", 7, sizeof (char), out);
	for (i = 0; i < info->num_instances; i++) {
		for (j = 0; j < info->num_attributes; j++) {
			switch (info->attributes[j]->type) {
			case ATTR_NUMERIC:
				val = info->instances[i]->data[j].fval;
				if (((double) (int) val) == val)
					sprintf (buf, "%i", (int) val);
				else
					sprintf (buf, "%f", val);
				break;
			case ATTR_NOMINAL:
				for (k = 0, nom_val =
				     info->attributes[j]->nom_info->first;
				     k < info->instances[i]->data[j].ival;
				     k++) {
					nom_val = nom_val->next;
				}
				strcpy (buf, nom_val->name);
				break;
			default:
				buf[0] = 0;
				break;
			}

			fwrite (buf, strlen (buf), sizeof (char), out);
			if ((j + 1) != info->num_attributes)
				fwrite (",", 1, sizeof (char), out);
		}
		fwrite ("\n", 1, sizeof (char), out);
	}

	fflush (out);
}

void release_read_info (arff_info_t * info)
{
	int i;
	nom_val_t *nom_val, *next_nom_val;

	free (info->relation_name);
	if (info->attributes != NULL) {
		for (i = 0; i < info->num_attributes; i++) {
			free (info->attributes[i]->name);
			if (info->attributes[i]->type == ATTR_NOMINAL) {
				for (nom_val =
				     info->attributes[i]->nom_info->first;
				     nom_val != NULL;
				     nom_val = next_nom_val) {
					free (nom_val->name);
					next_nom_val = nom_val->next;
					free (nom_val);
				}
				free (info->attributes[i]->nom_info);
			}
			free (info->attributes[i]);
		}
		free (info->attributes);
	}
	if (info->instances != NULL) {
		for (i = 0; i < info->num_instances; i++) {
			free (info->instances[i]->data);
			free (info->instances[i]);
		}
		free (info->instances);
	}
	free (info);
}

char *get_last_error ()
{
	return error_string;
}

int get_lineno ()
{
	return line_no;
}

arff_info_t *init (char *class_attribute_name)
{
	arff_info_t *info =
		(arff_info_t *) malloc_dbg (26, sizeof (arff_info_t));
	info->relation_name = NULL;
	info->num_attributes = 0;
	info->attributes = NULL;
	info->num_instances = 0;
	info->instances = NULL;
	info->class_index = -1;

	state = PARSE_STATE_NEWLINE;
	first_attr = NULL;
	first_inst = NULL;
	curr_instance = NULL;
	strcpy (error_string, "");
	class_name = class_attribute_name;
	line_no = 1;
	return info;
}

int lex (char *buf, arff_info_t * info)
{
	char *tok;
	char c;
	int i, r;
	token_t delimiter;

	for (;;) {
		if ((c = *buf++) == '\0')
			break;

		switch (c) {
		case ' ':
		case '\t':
		case '\r':
			// ignore
			break;
		case '%':
			/* a comment */
			do {
				c = *buf++;
			} while ((c != '\n') && (c != '\0'));
			buf--;
			break;
		case '\n':
			/* a new line */
			TRY (r, parse (TOKEN_NEWLINE, NULL, info));
			break;
		case ',':
			/* a comma */
			TRY (r, parse (TOKEN_COMMA, NULL, info));
			break;
		case '{':
			/* a left brace */
			TRY (r, parse (TOKEN_LEFTBRACE, NULL, info));
			break;
		case '}':
			/* a right brace */
			TRY (r, parse (TOKEN_RIGHTBRACE, NULL, info));
			break;
		default:
			/* something else */
			buf--;
			tok = read_token (&buf, &delimiter);

			/* check if it is a special token */
			for (i = 0; special[i].string != NULL; i++) {
				if (!stricmp (tok, special[i].string)) {
					TRY (r,
					     parse (special[i].type, tok,
						    info));
					break;
				}
			}
			if (special[i].string == NULL)
				TRY (r, parse (TOKEN_TOKEN, tok, info));
			if (delimiter < TOKEN_NUM_TYPES)
				TRY (r, parse (delimiter, NULL, info));
			break;
		}
	}

	parse_end (info);

	return 0;
}

char *read_token (char **pbuf, token_t * pdelimiter)
{
	char *tok = *pbuf;
	char c;

	do {
		c = *(*pbuf)++;
	} while ((c != ' ') && (c != '\t') && (c != '\n') && (c != '\r') &&
		 (c != '\0') && (c != ',') && (c != '{') && (c != '}'));
	*((*pbuf) - 1) = '\0';
	switch (c) {
	case '\n':
		*pdelimiter = TOKEN_NEWLINE;
		break;
	case ',':
		*pdelimiter = TOKEN_COMMA;
		break;
	case '{':
		*pdelimiter = TOKEN_LEFTBRACE;
		break;
	case '}':
		*pdelimiter = TOKEN_RIGHTBRACE;
		break;
	default:
		*pdelimiter = TOKEN_NUM_TYPES;
		break;
	}
	return tok;
}

int stricmp (const char *x, const char *y)
{
	char c1, c2;
	while ((*x != '\0') && (*y != '\0')) {
		c1 = *x++;
		c2 = *y++;
		if (LOWERCASE (c1) != LOWERCASE (c2))
			return 1;
	}
	return *x == *y ? 0 : 1;
}

int parse (token_t token, char *name, arff_info_t * info)
{
	DEBUGMSG (("Received token type %s, value = %s\n",
		   token_names[token], name == NULL ? "NULL" : name));
	DEBUGMSG (("  state = %s\n", state_names[state]));

	int i;
	parse_state_t nextstate = PARSE_STATE_ERROR;
	parse_table_t *table = parse_table[state];

	if (token == TOKEN_NEWLINE)
		line_no++;

	for (i = 0; table[i].token != TOKEN_NUM_TYPES; i++) {
		if ((table[i].token == token)
		    || (table[i].token == TOKEN_ANY)) {
			if (table[i].action == NULL) {
				nextstate = table[i].nextstate;
			} else {
				nextstate =
					table[i].action (token, name, info);
			}
			break;
		}
	}
	if (table[i].token == TOKEN_NUM_TYPES) {
		sprintf (error_string, "unexpected token: %s",
			 name == NULL ? token_names[token] : name);
	}
	DEBUGMSG (("  next state = %s\n", state_names[nextstate]));
	state = nextstate;

	return state == PARSE_STATE_ERROR ? 1 : 0;
}


parse_state_t parse_relation1 (token_t token, char *name, arff_info_t * info)
{
	parse_state_t r = PARSE_STATE_ERROR;
	if (name == NULL) {
		sprintf (error_string, "unexpected token: %s",
			 name == NULL ? token_names[token] : name);
	} else {
		info->relation_name =
			(char *) malloc_dbg (27, strlen (name) + 1);
		strcpy (info->relation_name, name);
		r = PARSE_STATE_GET_NEWLINE;
		DEBUGMSG (("  set relation name = %s\n", name));
	}
	return r;
}

parse_state_t parse_attribute1 (token_t token, char *name, arff_info_t * info)
{
	parse_state_t r = PARSE_STATE_ERROR;
	if (name == NULL) {
		sprintf (error_string, "unexpected token: %s",
			 name == NULL ? token_names[token] : name);
	} else {
		curr_attr = add_attribute (info, name);
		r = PARSE_STATE_ATTRIBUTE2;
	}
	return r;
}

parse_state_t parse_attribute2_numeric (token_t token, char *name,
					arff_info_t * info)
{
	curr_attr->type = ATTR_NUMERIC;
	DEBUGMSG (("  set attribute type = ATTR_NUMERIC\n"));
	return PARSE_STATE_NEWLINE;
}

parse_state_t parse_attribute2_nominal (token_t token, char *name,
					arff_info_t * info)
{
	curr_attr->type = ATTR_NOMINAL;
	curr_attr->nom_info =
		(nom_info_t *) malloc_dbg (28, sizeof (nom_info_t));
	curr_attr->nom_info->num_classes = 0;
	curr_attr->nom_info->first = NULL;
	DEBUGMSG (("  set attribute type = ATTR_NOMINAL\n"));
	return PARSE_STATE_ATTRIBUTE_NOM1;
}

parse_state_t parse_attribute_nom1 (token_t token, char *name,
				    arff_info_t * info)
{
	add_nominal_class (curr_attr->nom_info, name);
	return PARSE_STATE_ATTRIBUTE_NOM2;
}

parse_state_t parse_data (token_t token, char *name, arff_info_t * info)
{
	nom_val_t *x;
	parse_state_t r = PARSE_STATE_DATA2;

	if (curr_instance == NULL) {
		curr_instance = add_instance (info);
		curr_data = 0;
		curr_attr = first_attr;
	}

	if (curr_attr == NULL) {
		sprintf (error_string, "too many data values given");
		r = PARSE_STATE_ERROR;
	} else {
		switch (curr_attr->type) {
		case ATTR_NUMERIC:
			curr_instance->data[curr_data++].fval = atof (name);
			DEBUGMSG (("  add numeric data, val = %f\n",
				   atof (name)));
			break;

		case ATTR_NOMINAL:
			for (x = curr_attr->nom_info->first; x != NULL;
			     x = x->next) {
				if (!strcmp (name, x->name)) {
					curr_instance->data[curr_data++].
						ival = x->val;
					DEBUGMSG (("  add nominal data, val = %s, index = %i\n", x->name, x->val));
					break;
				}
			}
			if (x == NULL) {
				sprintf (error_string,
					 "unknown nominal class value: %s",
					 name);
				r = PARSE_STATE_ERROR;
			}
			break;
		default:
			break;
		}

		curr_attr = curr_attr->next;
	}

	return r;
}

parse_state_t parse_data_end (token_t token, char *name, arff_info_t * info)
{
	parse_state_t r = PARSE_STATE_DATA1;
	if (curr_data != info->num_attributes) {
		sprintf (error_string, "not enough data values given");
		r = PARSE_STATE_ERROR;
	}
	curr_instance = NULL;
	return r;
}

void parse_end (arff_info_t * info)
{
	int i;
	attr_info_t *attr_info;
	instance_t *instance;
	info->attributes =
		(attr_info_t **) malloc_dbg (29, sizeof (attr_info_t *) *
					     info->num_attributes);
	for (i = 0, attr_info = first_attr; attr_info != NULL;
	     attr_info = attr_info->next, i++) {
		info->attributes[i] = attr_info;
	}
	info->instances =
		(instance_t **) malloc_dbg (30, sizeof (instance_t *) *
					    info->num_instances);
	for (i = 0, instance = first_inst; instance != NULL;
	     instance = instance->next, i++) {
		info->instances[i] = instance;
	}
}

attr_info_t *add_attribute (arff_info_t * info, char *name)
{
	attr_info_t *r =
		(attr_info_t *) malloc_dbg (31, sizeof (attr_info_t)), *x;

	r->name = (char *) malloc_dbg (32, strlen (name) + 1);
	strcpy (r->name, name);
	r->type = ATTR_NUM_TYPES;
	r->nom_info = NULL;
	r->next = NULL;

	if (first_attr == NULL) {
		first_attr = r;
	} else {
		for (x = first_attr; x->next != NULL; x = x->next);
		x->next = r;
	}

	if (!stricmp (name, class_name)) {
		info->class_index = info->num_attributes;
	}

	info->num_attributes++;

	DEBUGMSG (("  create attribute %s, total = %i\n", name,
		   info->num_attributes));

	return r;
}

void add_nominal_class (nom_info_t * info, char *name)
{
	nom_val_t *r = (nom_val_t *) malloc_dbg (33, sizeof (nom_val_t)), *x;
	r->name = (char *) malloc_dbg (34, strlen (name) + 1);
	strcpy (r->name, name);
	r->val = info->num_classes++;
	r->next = NULL;

	if (info->first == NULL) {
		info->first = r;
	} else {
		for (x = info->first; x->next != NULL; x = x->next);
		x->next = r;
	}

	DEBUGMSG (("  create nominal class %s, index = %i\n", name, r->val));
}

instance_t *add_instance (arff_info_t * info)
{
	instance_t *r =
		(instance_t *) malloc_dbg (35, sizeof (instance_t)), *x;
	r->data =
		(data_t *) malloc_dbg (36,
				       sizeof (data_t) *
				       info->num_attributes);
	r->next = NULL;

	if (first_inst == NULL) {
		first_inst = r;
	} else {
		for (x = first_inst; x->next != NULL; x = x->next);
		x->next = r;
	}
	info->num_instances++;

	DEBUGMSG (("  create instance, total = %i\n", info->num_instances));

	return r;
}

#ifdef ARFF_DEBUG_MAIN
void print_arff_info (arff_info_t * info)
{
	nom_val_t *nom_val;
	int i, j, k;

	printf ("relation name = %s\n", info->relation_name);
	printf ("  num attributes = %i\n", info->num_attributes);
	printf ("  num instances = %i\n", info->num_instances);
	printf ("\n");
	printf ("  attributes:\n");
	for (i = 0; i < MIN (info->num_attributes, 30); i++) {
		printf ("    %s, ", info->attributes[i]->name);
		switch (info->attributes[i]->type) {
		case ATTR_NUMERIC:
			printf ("numeric");
			break;
		case ATTR_NOMINAL:
			printf ("{");
			for (nom_val = info->attributes[i]->nom_info->first;
			     nom_val != NULL; nom_val = nom_val->next) {
				printf ("%s%s", nom_val->name,
					nom_val->next == NULL ? "" : ", ");
			}
			printf ("}");
			break;
		}
		if (i == info->class_index) {
			printf (" *class attribute");
		}
		printf ("\n");
	}
	if (i < info->num_attributes)
		printf ("...\n");
	printf ("\n");
	printf ("  data:\n");
	for (j = 0; j < MIN (info->num_instances, 30); j++) {
		printf ("    ");
		for (i = 0; i < MIN (info->num_attributes, 16); i++) {
			switch (info->attributes[i]->type) {
			case ATTR_NUMERIC:
				printf ("%.1f",
					info->instances[j]->data[i].fval);
				break;
			case ATTR_NOMINAL:
				for (k = 0, nom_val =
				     info->attributes[i]->nom_info->first;
				     k < info->instances[j]->data[i].ival;
				     k++) {
					nom_val = nom_val->next;
				}
				printf ("%s(%i)", nom_val->name,
					info->instances[j]->data[i].ival);
				break;
			}
			if ((i + 1) < info->num_attributes)
				printf (", ");
		}
		if (i < info->num_attributes)
			printf ("...");
		printf ("\n");
	}
	if (j < info->num_attributes)
		printf ("...\n");
}

int main (int argc, char **argv)
{
	int i;
	arff_info_t *info;
	for (i = 1; i < (argc - 1); i += 2) {
		info = read_arff (argv[i], argv[i + 1]);
		if (info == NULL) {
			printf ("read_arff failed: %s, lineno=%i\n",
				get_last_error (), get_lineno ());
		} else {
			print_arff_info (info);
			release_read_info (info);
		}
	}
	return 0;
}
#endif
