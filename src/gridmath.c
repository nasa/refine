/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
#include <stdlib.h>
#include <math.h>
#include "gridmath.h"

void gridSubtractVector(double *v1, double *v2, double *result)
{
  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];
}

double gridDotProduct(double *v1, double *v2)
{
  return ( v1[0]*v2[0] +  v1[1]*v2[1] + v1[2]*v2[2] );
}

void gridCrossProduct(double *edge1, double *edge2, double *norm)
{
  norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; 
  norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; 
  norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 
}

void gridVectorNormalize(double *norm)
{
  double length;
  length = sqrt(gridDotProduct(norm,norm));
  if (length > 0 ) {
    norm[0] /= length;
    norm[1] /= length;
    norm[2] /= length;
  }
}

