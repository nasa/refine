
/* Computes geometric shape functions for elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:Mike.Park@NASA.Gov
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include "gridshape.h"

Grid* gridShapeJacobian1(Grid *grid, 
			 double *n0, double *n1,double *n2, double *n3, 
			 double *where, double *j )
{
  double gphi[3];
  j[0] = 0.0;
  j[1] = 0.0;
  j[2] = 0.0;
  j[3] = 0.0;
  j[4] = 0.0;
  j[5] = 0.0;
  j[6] = 0.0;
  j[7] = 0.0;
  j[8] = 0.0;

  gphi[0] = -1.0;
  gphi[1] = -1.0;
  gphi[2] = -1.0;
  j[0] += n0[0] * gphi[0];
  j[1] += n0[0] * gphi[1];
  j[2] += n0[0] * gphi[2];
  j[3] += n0[1] * gphi[0];
  j[4] += n0[1] * gphi[1];
  j[5] += n0[1] * gphi[2];
  j[6] += n0[2] * gphi[0];
  j[7] += n0[2] * gphi[1];
  j[8] += n0[2] * gphi[2];

  gphi[0] = 1.0;
  gphi[1] = 0.0;
  gphi[2] = 0.0;
  j[0] += n1[0] * gphi[0];
  j[1] += n1[0] * gphi[1];
  j[2] += n1[0] * gphi[2];
  j[3] += n1[1] * gphi[0];
  j[4] += n1[1] * gphi[1];
  j[5] += n1[1] * gphi[2];
  j[6] += n1[2] * gphi[0];
  j[7] += n1[2] * gphi[1];
  j[8] += n1[2] * gphi[2];

  gphi[0] = 0.0;
  gphi[1] = 1.0;
  gphi[2] = 0.0;
  j[0] += n2[0] * gphi[0];
  j[1] += n2[0] * gphi[1];
  j[2] += n2[0] * gphi[2];
  j[3] += n2[1] * gphi[0];
  j[4] += n2[1] * gphi[1];
  j[5] += n2[1] * gphi[2];
  j[6] += n2[2] * gphi[0];
  j[7] += n2[2] * gphi[1];
  j[8] += n2[2] * gphi[2];

  gphi[0] = 0.0;
  gphi[1] = 0.0;
  gphi[2] = 1.0;
  j[0] += n3[0] * gphi[0];
  j[1] += n3[0] * gphi[1];
  j[2] += n3[0] * gphi[2];
  j[3] += n3[1] * gphi[0];
  j[4] += n3[1] * gphi[1];
  j[5] += n3[1] * gphi[2];
  j[6] += n3[2] * gphi[0];
  j[7] += n3[2] * gphi[1];
  j[8] += n3[2] * gphi[2];

  return grid;
}
