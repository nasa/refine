
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
			 double *n0,double *n1,double *n2, double *n3, 
			 double *where, double *jacobian )
{
  jacobian[0] = 1.0;
  jacobian[1] = 0.0;
  jacobian[2] = 0.0;
  jacobian[3] = 0.0;
  jacobian[4] = 1.0;
  jacobian[5] = 0.0;
  jacobian[6] = 0.0;
  jacobian[7] = 0.0;
  jacobian[8] = 1.0;
  return grid;
}
