/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMATH_H
#define GRIDMATH_H

#include "master_header.h"

BEGIN_C_DECLORATION

void gridSubtractVector(double *v1, double *v2, double *result);
double gridDotProduct(double *v1, double *v2);
void gridCrossProduct(double *norm, double *edge1, double *edge2);
void gridVectorNormalize(double *norm);

END_C_DECLORATION

#endif /* GRIDMATH_H */
