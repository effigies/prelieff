#ifndef _ARFF_H
#define _ARFF_H

typedef struct _nom_val_t {
	char *name;
	int val;
	struct _nom_val_t *next;
} nom_val_t;

typedef struct {
	int num_classes;
	nom_val_t *first;
} nom_info_t;

typedef enum {
	ATTR_NUMERIC,
	ATTR_NOMINAL,
	ATTR_NUM_TYPES
} attr_t;

typedef struct _attr_info_t {
	char *name;
	attr_t type;
	nom_info_t *nom_info;
	struct _attr_info_t *next;
} attr_info_t;

typedef union {
	float fval;
	int ival;
	char *sval;
} data_t;

typedef struct _instance_t {
	data_t *data;
	struct _instance_t *next;
} instance_t;

typedef struct {
	char *relation_name;
	int num_attributes;
	attr_info_t **attributes;
	int num_instances;
	instance_t **instances;
	int class_index;
} arff_info_t;

arff_info_t *read_arff (char *filename, char *class_attribute_name);
void write_arff (arff_info_t * info, FILE * out);
void release_read_info (arff_info_t * info);
char *get_last_error ();
int get_lineno ();

#endif
