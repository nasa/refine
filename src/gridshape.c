
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

Grid* gridShapeJacobian2(Grid *grid, 
			 double *n0, double *n1, double *n2, double *n3,
			 double *e01, double *e02, double *e03,
			 double *e12, double *e13, double *e23,
			 double *w, double *j )
{
  double x, y, z;
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

  x = w[0];
  y = w[1];
  z = w[2];

  gphi[0] = -3.0+4.0*x+4.0*y+4.0*z;
  gphi[1] = -3.0+4.0*x+4.0*y+4.0*z;
  gphi[2] = -3.0+4.0*x+4.0*y+4.0*z;
  j[0] += n0[0] * gphi[0];
  j[1] += n0[0] * gphi[1];
  j[2] += n0[0] * gphi[2];
  j[3] += n0[1] * gphi[0];
  j[4] += n0[1] * gphi[1];
  j[5] += n0[1] * gphi[2];
  j[6] += n0[2] * gphi[0];
  j[7] += n0[2] * gphi[1];
  j[8] += n0[2] * gphi[2];

  gphi[0] = -1.0+4.0*x;
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
  gphi[1] = -1.0+4.0*y;
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
  gphi[2] = -1.0+4.0*y;
  j[0] += n3[0] * gphi[0];
  j[1] += n3[0] * gphi[1];
  j[2] += n3[0] * gphi[2];
  j[3] += n3[1] * gphi[0];
  j[4] += n3[1] * gphi[1];
  j[5] += n3[1] * gphi[2];
  j[6] += n3[2] * gphi[0];
  j[7] += n3[2] * gphi[1];
  j[8] += n3[2] * gphi[2];

  gphi[0] = 4.0-8.0*x-4.0*y-4.0*z;
  gphi[1] = -4.0*x;
  gphi[2] = -4.0*x;
  j[0] += e01[0] * gphi[0];
  j[1] += e01[0] * gphi[1];
  j[2] += e01[0] * gphi[2];
  j[3] += e01[1] * gphi[0];
  j[4] += e01[1] * gphi[1];
  j[5] += e01[1] * gphi[2];
  j[6] += e01[2] * gphi[0];
  j[7] += e01[2] * gphi[1];
  j[8] += e01[2] * gphi[2];

  gphi[0] = -4.0*y;
  gphi[1] = 4.0-4.0*x-8.0*y-4.0*z;
  gphi[2] = -4.0*y;
  j[0] += e02[0] * gphi[0];
  j[1] += e02[0] * gphi[1];
  j[2] += e02[0] * gphi[2];
  j[3] += e02[1] * gphi[0];
  j[4] += e02[1] * gphi[1];
  j[5] += e02[1] * gphi[2];
  j[6] += e02[2] * gphi[0];
  j[7] += e02[2] * gphi[1];
  j[8] += e02[2] * gphi[2];

  gphi[0] = -4.0*z;
  gphi[1] = -4.0*z;
  gphi[2] = 4.0-4.0*x-4.0*y-8.0*z;
  j[0] += e03[0] * gphi[0];
  j[1] += e03[0] * gphi[1];
  j[2] += e03[0] * gphi[2];
  j[3] += e03[1] * gphi[0];
  j[4] += e03[1] * gphi[1];
  j[5] += e03[1] * gphi[2];
  j[6] += e03[2] * gphi[0];
  j[7] += e03[2] * gphi[1];
  j[8] += e03[2] * gphi[2];

  gphi[0] = 4.0*y;
  gphi[1] = 4.0*x;
  gphi[2] = 0.0;
  j[0] += e12[0] * gphi[0];
  j[1] += e12[0] * gphi[1];
  j[2] += e12[0] * gphi[2];
  j[3] += e12[1] * gphi[0];
  j[4] += e12[1] * gphi[1];
  j[5] += e12[1] * gphi[2];
  j[6] += e12[2] * gphi[0];
  j[7] += e12[2] * gphi[1];
  j[8] += e12[2] * gphi[2];

  gphi[0] = 4.0*z;
  gphi[1] = 0.0;
  gphi[2] = 4.0*x;
  j[0] += e13[0] * gphi[0];
  j[1] += e13[0] * gphi[1];
  j[2] += e13[0] * gphi[2];
  j[3] += e13[1] * gphi[0];
  j[4] += e13[1] * gphi[1];
  j[5] += e13[1] * gphi[2];
  j[6] += e13[2] * gphi[0];
  j[7] += e13[2] * gphi[1];
  j[8] += e13[2] * gphi[2];

  gphi[0] = 0.0;
  gphi[1] = 4.0*z;
  gphi[2] = 4.0*y;
  j[0] += e23[0] * gphi[0];
  j[1] += e23[0] * gphi[1];
  j[2] += e23[0] * gphi[2];
  j[3] += e23[1] * gphi[0];
  j[4] += e23[1] * gphi[1];
  j[5] += e23[1] * gphi[2];
  j[6] += e23[2] * gphi[0];
  j[7] += e23[2] * gphi[1];
  j[8] += e23[2] * gphi[2];

  return grid;
}
