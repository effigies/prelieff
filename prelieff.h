#ifndef _PRELIEFF_H
#define _PRELIEFF_H

#include "arff.h"
#include "java.h"

void buildEvaluator (arff_info_t * data, double *weights);
double evaluateAttribute (int attribute);

void resetOptions ();
void setSigma (int s);
int getSigma ();
void setNumNeighbours (int n);
int getNumNeighbours ();
void setSeed (int s);
int getSeed ();
void setSampleSize (int s);
int getSampleSize ();
void setWeightByDistance (boolean b);
boolean getWeightByDistance ();
double getTotalTime ();
void setVersion (int version);
void setDifference (int);

#endif
