
/*
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 *    ReliefFAttributeEval.java
 *    Copyright (C) 1999 University of Waikato, Hamilton, New Zealand
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include "arff.h"
#include "java.h"
#include "index_sort.h"
#include "util.h"
#ifndef NO_MPI
#include "mpi.h"
#endif

#define SMALL       1e-6
#define EQ(a,b)     (((a-b)<SMALL) && ((b-a)<SMALL))

/** The training instances */
static arff_info_t *m_trainInstances;

/** The class index */
static int m_classIndex;

/** The number of attributes */
static int m_numAttribs;

/** The number of instances */
static int m_numInstances;

/** The attributes */
static attr_info_t **m_attributes;

/** The instances */
static instance_t **m_instances;

/** The number of classes if class is nominal */
static int m_numClasses;

/** Holds the weights that relief assigns to attributes */
static double *m_weights;
static double *m_finalWeights;

/** Prior class probabilities (discrete class case) */
static double *m_classProbs;

/** 
  * The number of instances to sample when estimating attributes
  * default == -1, use all instances
*/
static int m_sampleM;

/** The number of nearest hits/misses */
static int m_Knn;

/** k nearest scores + instance indexes for n classes */
static double ***m_karray;

/** Upper bound for numeric attributes */
static double *m_maxArray;

/** Lower bound for numeric attributes */
static double *m_minArray;

/** Keep track of the farthest instance for each class */
static double *m_worst;

/** Index in the m_karray of the farthest instance for each class */
static int *m_index;

/** Number of nearest neighbours stored of each class */
static int *m_stored;

/** Random number seed used for sampling instances */
static int m_seed;

/**
  *  used to (optionally) weight nearest neighbours by their distance
  *  from the instance in question. Each entry holds 
  *  exp(-((rank(r_i, i_j)/sigma)^2)) where rank(r_i,i_j) is the rank of
  *  instance i_j in a sequence of instances ordered by the distance
  *  from r_i. sigma is a user defined parameter, default=20
**/
static double *m_weightsByRank;
static int m_sigma;

/** Weight by distance rather than equal weights */
static boolean m_weightByDistance;

static double m_totalTime;

static double *tempDistClass;
static double *tempDistAtt;
static int *tempSortedClass;
static int **tempSortedAtt;
static double *distNormAtt;
static int *m_attributeRank;
static int m_numExcludedAttributes;

static int m_version = 0;	// version of the algorithm to use
int m_difference = 0;

void updateMinMax (instance_t * instance);
void findKHitMiss (int instNum);
void updateWeightsDiscreteClass (int instNum);

void setSigma (int s)
{
	m_sigma = s;
}

int getSigma ()
{
	return m_sigma;
}

void setNumNeighbours (int n)
{
	m_Knn = n;
}

int getNumNeighbours ()
{
	return m_Knn;
}

void setSeed (int s)
{
	m_seed = s;
}

int getSeed ()
{
	return m_seed;
}

void setSampleSize (int s)
{
	m_sampleM = s;
}

int getSampleSize ()
{
	return m_sampleM;
}

void setWeightByDistance (boolean b)
{
	m_weightByDistance = b;
}

boolean getWeightByDistance ()
{
	return m_weightByDistance;
}

double getTotalTime ()
{
	return m_totalTime;
}

void setVersion (int version)
{
	m_version = version;
}

void setDifference (int difference)
{
	m_difference = difference;
}

/**
  * Initializes a ReliefF attribute evaluator. 
  *
  * @param data set of instances serving as training data 
  * @throws Exception if the evaluator has not been 
  * generated successfully
  */
void buildEvaluator (arff_info_t * data, double *weights)
{

	int i, j, k, z, totalInstances;
	int num_nodes, my_rank;
	double t0, t1;
#ifdef PRINT_STATUS
	char buf[100];
#endif

	srand (time (NULL));

#ifdef NO_MPI
	t0 = (double) clock () / CLOCKS_PER_SEC;
	num_nodes = 1;
	my_rank = 0;
#else
	t0 = MPI_Wtime ();
	MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
#endif

	m_trainInstances = data;
	m_classIndex = m_trainInstances->class_index;
	m_numAttribs = m_trainInstances->num_attributes;
	m_numInstances = m_trainInstances->num_instances;
	m_attributes = m_trainInstances->attributes;
	m_instances = m_trainInstances->instances;

	m_numClasses = m_attributes[m_classIndex]->nom_info->num_classes;

	if (m_weightByDistance)	// set up the rank based weights
	{
		m_weightsByRank = (double *) malloc_dbg (1, sizeof (double) * m_Knn);

		for (i = 0; i < m_Knn; i++) {
			m_weightsByRank[i] =
				exp (-
				     ((i / (double) m_sigma) *
				      (i / (double) m_sigma)));
		}
	}
	// the final attribute weights
#ifdef NO_MPI
	m_weights = m_finalWeights = weights;
#else
	m_weights = (double *) malloc_dbg (2, sizeof (double) * m_numAttribs);
	m_finalWeights = weights;
#endif


	m_attributeRank = (int *) malloc_dbg (3, sizeof (int) * m_numAttribs);
	
	for (i = 0; i < m_numAttribs; i++) {
		m_attributeRank[i] = i;
	}
	m_numExcludedAttributes = 0;

	for (i = 0; i < m_numAttribs; i++) {
		m_weights[i] = m_finalWeights[i] = 0.0;
	}
	// num classes (1 for numeric class) knn neighbours, 
	// and 0 = distance, 1 = instance index
	m_karray = (double ***) malloc_dbg (4, sizeof (double **) * m_numClasses);
	for (i = 0; i < m_numClasses; i++) {
		m_karray[i] = (double **) malloc_dbg (5, sizeof (double *) * m_Knn);
		for (j = 0; j < m_Knn; j++) {
			m_karray[i][j] =
				(double *) malloc_dbg (6, sizeof (double) * 2);
		}
	}

	m_classProbs = (double *) malloc_dbg (7, sizeof (double) * m_numClasses);

	for (i = 0; i < m_numInstances; i++) {
		m_classProbs[m_instances[i]->data[m_classIndex].ival]++;
	}

	for (i = 0; i < m_numClasses; i++) {
		m_classProbs[i] /= m_numInstances;
	}

	m_worst = (double *) malloc_dbg (8, sizeof (double) * m_numClasses);
	m_index = (int *) malloc_dbg (9, sizeof (int) * m_numClasses);
	m_stored = (int *) malloc_dbg (10, sizeof (int) * m_numClasses);
	m_minArray = (double *) malloc_dbg (11, sizeof (double) * m_numAttribs);
	m_maxArray = (double *) malloc_dbg (12, sizeof (double) * m_numAttribs);

	tempDistClass = (double *) malloc_dbg (13, sizeof (double) * m_numAttribs);
	tempDistAtt = (double *) malloc_dbg (14, sizeof (double) * m_numAttribs);
	tempSortedClass = (int *) malloc_dbg (15, sizeof (int) * m_numAttribs);
	tempSortedAtt = (int **) malloc_dbg (16, sizeof (int *) * m_numClasses);
	distNormAtt = (double *) malloc_dbg (17, sizeof (double) * m_numClasses);

	for (i = 0; i < m_numClasses; i++) {
		if (i != m_classIndex)	// already done cl
		{
			tempSortedAtt[i] =
				(int *) malloc_dbg (18, sizeof (int) * m_numAttribs);
		}
	}

	for (i = 0; i < m_numAttribs; i++) {
		m_minArray[i] = m_maxArray[i] = DBL_MAX;
	}

	for (i = 0; i < m_numInstances; i++) {
		updateMinMax (m_instances[i]);
	}

	if ((m_sampleM > m_numInstances) || (m_sampleM < 0)) {
		totalInstances = m_numInstances;
	} else {
		totalInstances = m_sampleM;
	}

	// process each instance, updating attribute weights
	for (i = 0; i < totalInstances; i += num_nodes) {
#ifdef PRINT_STATUS
		sprintf (buf, "%05i:%05i", i, totalInstances);
		printf ("%s\n", buf);
		fflush (stdout);
#endif

		if (totalInstances == m_numInstances) {
			z = i + my_rank;
		} else {
			z = rand () % m_numInstances;
		}

		if (z < 0) {
			z *= -1;
		}
		// first clear the knn and worst index stuff for the classes
		for (j = 0; j < m_numClasses; j++) {
			m_index[j] = m_stored[j] = 0;

			for (k = 0; k < m_Knn; k++) {
				m_karray[j][k][0] = m_karray[j][k][1] = 0;
			}
		}

		if (z < totalInstances) {
			findKHitMiss (z);

			updateWeightsDiscreteClass (z);
		}

		if (m_version == 1) {
#ifndef NO_MPI
			MPI_Allreduce (m_weights, m_finalWeights,
				       m_numAttribs, MPI_DOUBLE, MPI_SUM,
				       MPI_COMM_WORLD);
			if (my_rank == 0) {
				memcpy (m_weights, m_finalWeights,
					m_numAttribs * sizeof (double));
			} else {
				memset (m_weights, 0,
					m_numAttribs * sizeof (double));
			}
#endif
			index_sort (m_attributeRank, m_finalWeights,
				    m_numAttribs);
			m_numExcludedAttributes++;
		}
	}

#ifndef NO_MPI
	if (m_version != 1) {
		MPI_Reduce (m_weights, m_finalWeights, m_numAttribs,
			    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	free (m_weights);
#endif

	// now scale weights by 1/m_numInstances (nominal class) or
	// calculate weights numeric class
	for (i = 0; i < m_numAttribs; i++) {
		if (i != m_classIndex) {
			m_finalWeights[i] *= (1.0 / (double) totalInstances);
		}
	}

	if (m_weightByDistance)
		free (m_weightsByRank);
	for (i = 0; i < m_numClasses; i++) {
		for (j = 0; j < m_Knn; j++)
			free (m_karray[i][j]);
		free (m_karray[i]);
	}
	free (m_karray);
	free (m_classProbs);
	free (m_worst);
	free (m_index);
	free (m_stored);
	free (m_minArray);
	free (m_maxArray);

	for (i = 0; i < m_numClasses; i++) {
		if (i != m_classIndex)	// already done cl
		{
			free (tempSortedAtt[i]);
		}
	}

	free (tempDistClass);
	free (tempDistAtt);
	free (tempSortedClass);
	free (tempSortedAtt);
	free (distNormAtt);
	free (m_attributeRank);

#ifdef NO_MPI
	t1 = (double) clock () / CLOCKS_PER_SEC;
#else
	t1 = MPI_Wtime ();
#endif

	m_totalTime = t1 - t0;
}


/**
  * Evaluates an individual attribute using ReliefF's instance based approach.
  * The actual work is done by buildEvaluator which evaluates all features.
  *
  * @param attribute the index of the attribute to be evaluated
  * @throws Exception if the attribute could not be evaluated
  */
double evaluateAttribute (int attribute)
{
	return m_finalWeights[attribute];
}


/**
  * Reset options to their default values
  */
void resetOptions ()
{
	m_trainInstances = null;
	m_sampleM = -1;
	m_Knn = 10;
	m_sigma = 2;
	m_weightByDistance = false;
	m_seed = 1;
}


/**
  * Normalizes a given value of a numeric attribute.
  *
  * @param x the value to be normalized
  * @param i the attribute's index
  * @return the normalized value
  */
double norm (double x, int i)
{
	if ((m_minArray[i] == DBL_MAX) || EQ (m_maxArray[i], m_minArray[i])) {
		return 0;
	} else {
		return (x - m_minArray[i]) / (m_maxArray[i] - m_minArray[i]);
	}
}


/**
  * Updates the minimum and maximum values for all the attributes
  * based on a new instance.
  *
  * @param instance the new instance
  */
void updateMinMax (instance_t * instance)
{
	int j;
	for (j = 0; j < m_numAttribs; j++) {
		if (m_attributes[j]->type == ATTR_NUMERIC) {
			if (m_minArray[j] == DBL_MAX) {
				m_minArray[j] = instance->data[j].fval;
				m_maxArray[j] = instance->data[j].fval;
			} else {
				if (instance->data[j].fval < m_minArray[j]) {
					m_minArray[j] =
						instance->data[j].fval;
				} else {
					if (instance->data[j].fval >
					    m_maxArray[j]) {
						m_maxArray[j] =
							instance->data[j].
							fval;
					}
				}
			}
		}
	}
}

/**
  * Computes the difference between two given attribute
  * values.
  */
double difference (int index, data_t *dat1, data_t *dat2)
{
	double x;
	switch (m_attributes[index]->type) {
	case ATTR_NOMINAL:
		if( m_difference == 0 ) {
			return (dat1[index].ival != dat2[index].ival);
		} else {
			// This option is for when your nominal values differ by their distance
			// in the value list. e.g. AA Aa aa
			return abs(dat1[index].ival - dat2[index].ival);
		}
	case ATTR_NUMERIC:
		// If attribute is numeric
		x = norm (dat1[index].fval, index) - norm (dat2[index].fval, index);
		return x >= 0 ? x : -x;
	default:
		return 0;
	}
}

/**
  * Calculates the distance between two instances
  *
  * @param first the first instance
  * @param second the second instance
  * @return the distance between the two given instances, between 0 and 1
  */
double distance (instance_t * first, instance_t * second)
{

	double distance = 0, diff;
	int i, a;

	for (i = 0; i < (m_numAttribs - m_numExcludedAttributes); i++) {
		a = m_attributeRank[i];
		if (a == m_classIndex) {
			continue;
		}
		diff = difference (a, first->data, second->data);
		//      distance += diff * diff;
		distance += diff;
	}

	//    return Math.sqrt(distance / m_NumAttributesUsed);
	return distance;
}

/**
  * update attribute weights given an instance when the class is discrete
  *
  * @param instNum the index of the instance to use when updating weights
  */
void updateWeightsDiscreteClass (int instNum)
{
	int i, j, k, l;
	int cl;
	double temp_diff, w_norm = 1.0;
	double distNormClass = 1.0;

	// store the indexes (sparse instances) of non-zero elements
	instance_t *inst = m_instances[instNum];

	// get the class of this instance
	cl = inst->data[m_classIndex].ival;

	// sort nearest neighbours and set up normalization variables
	if (m_weightByDistance) {
		// do class (hits) first
		// sort the distances

		for (j = 0, distNormClass = 0; j < m_stored[cl]; j++) {
			// copy the distances
			tempDistClass[j] = m_karray[cl][j][0];
			// sum normalizer
			distNormClass += m_weightsByRank[j];
		}

		index_sort (tempSortedClass, tempDistClass, m_stored[cl]);

		for (k = 0; k < m_numClasses; k++) {
			if (k != cl)	// already done cl
			{
				// sort the distances
				for (j = 0, distNormAtt[k] = 0;
				     j < m_stored[k]; j++) {
					// copy the distances
					tempDistAtt[j] = m_karray[k][j][0];
					// sum normalizer
					distNormAtt[k] += m_weightsByRank[j];
				}

				index_sort (tempSortedAtt[k], tempDistAtt,
					    m_stored[k]);
			}
		}
	}

	if (m_numClasses > 2) {
		// the amount of probability space left after removing the
		// probability of this instance's class value
		w_norm = (1.0 - m_classProbs[cl]);
	}
	// do the k nearest hits of the same class
	for (j = 0, temp_diff = 0.0; j < m_stored[cl]; j++) {
		instance_t *cmp;
		cmp = (m_weightByDistance)
			? m_instances[(int)
				      m_karray[cl][tempSortedClass[j]][1]]
			: m_instances[(int) m_karray[cl][j][1]];

		for (k = 0; k < m_numAttribs; k++) {
			if (k == m_classIndex) {
				continue;
			}

			i = k;
			temp_diff =
				difference (k, inst->data,
					    cmp->data);

			if (m_weightByDistance) {
				temp_diff *=
					(m_weightsByRank[j] / distNormClass);
			} else {
				if (m_stored[cl] > 0) {
					temp_diff /= (double) m_stored[cl];
				}
			}
			m_weights[i] -= temp_diff;

		}
	}


	// now do k nearest misses from each of the other classes
	temp_diff = 0.0;

	for (k = 0; k < m_numClasses; k++) {
		if (k != cl)	// already done cl
		{
			for (j = 0; j < m_stored[k]; j++) {
				instance_t *cmp;
				cmp = (m_weightByDistance)
					? m_instances[(int)
						      m_karray[k]
						      [tempSortedAtt[k][j]]
						      [1]]
					: m_instances[(int)
						      m_karray[k][j][1]];

				for (l = 0; l < m_numAttribs; l++) {
					if (l == m_classIndex) {
						continue;
					}
					i = l;
					temp_diff =
						difference (l,
							    inst->data,
							    cmp->data);

					if (m_weightByDistance) {
						temp_diff *=
							(m_weightsByRank[j] /
							 distNormAtt[k]);
					} else {
						if (m_stored[k] > 0) {
							temp_diff /= (double)
								m_stored[k];
						}
					}
					if (m_numClasses > 2) {
						m_weights[i] +=
							((m_classProbs[k] /
							  w_norm) *
							 temp_diff);
					} else {
						m_weights[i] += temp_diff;
					}
				}
			}
		}
	}

}

/**
  * Find the K nearest instances to supplied instance if the class is numeric,
  * or the K nearest Hits (same class) and Misses (K from each of the other
  * classes) if the class is discrete.
  *
  * @param instNum the index of the instance to find nearest neighbours of
  */
void findKHitMiss (int instNum)
{
	int i, j;
	int cl;
	double ww;
	double temp_diff = 0.0;
	instance_t *thisInst = m_instances[instNum];

	for (i = 0; i < m_numInstances; i++) {
		if (i != instNum) {
			instance_t *cmpInst = m_instances[i];
			temp_diff = distance (cmpInst, thisInst);

			// class of this training instance
			cl = cmpInst->data[m_classIndex].ival;

			// add this diff to the list for the class of this instance
			if (m_stored[cl] < m_Knn) {
				m_karray[cl][m_stored[cl]][0] = temp_diff;
				m_karray[cl][m_stored[cl]][1] = i;
				m_stored[cl]++;

				// note the worst diff for this class
				for (j = 0, ww = -1.0; j < m_stored[cl]; j++) {
					if (m_karray[cl][j][0] > ww) {
						ww = m_karray[cl][j][0];
						m_index[cl] = j;
					}
				}

				m_worst[cl] = ww;
			} else
				/* if we already have stored knn for this class then check to
				   see if this instance is better than the worst */
			{
				if (temp_diff < m_karray[cl][m_index[cl]][0]) {
					m_karray[cl][m_index[cl]][0] =
						temp_diff;
					m_karray[cl][m_index[cl]][1] = i;

					for (j = 0, ww = -1.0;
					     j < m_stored[cl]; j++) {
						if (m_karray[cl][j][0] > ww) {
							ww = m_karray[cl][j]
								[0];
							m_index[cl] = j;
						}
					}

					m_worst[cl] = ww;
				}
			}
		}
	}
}
