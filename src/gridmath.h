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

#define gridSubtractVector(v1,v2,result) \
result[0] = v1[0] - v2[0]; result[1] = v1[1] - v2[1]; result[2] = v1[2] - v2[2];

#define gridDotProduct(v1, v2) ( v1[0]*v2[0] +  v1[1]*v2[1] + v1[2]*v2[2] )

void gridCrossProduct(double *norm, double *edge1, double *edge2);
double gridVectorLength(double *norm);
void gridVectorNormalize(double *norm);

END_C_DECLORATION

#endif /* GRIDMATH_H */
