
/* Computes geometric shape functions for elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:Mike.Park@NASA.Gov 
 */
  
/* $Id$ */

#ifndef GRIDSHAPE_H
#define GRIDSHAPE_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridShapeJacobian1(Grid *, double *n0,double *n1,double *n2, double *n3, 
			 double *where, double *jacobian );
END_C_DECLORATION

#endif /* GRIDSHAPE_H */
