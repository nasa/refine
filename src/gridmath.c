/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
#include <stdlib.h>
#include <math.h>
#include "gridmath.h"

double gridVectorLength(double *v)
{
  return sqrt( gridDotProduct(v, v) );
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

