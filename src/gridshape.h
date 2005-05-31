
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

Grid *gridCurvedEdgeMidpoint(Grid *,int node0, int node1, double *midpoint);

double gridMinCellJacDet2(Grid *, int *nodes);

Grid *gridPlotMinDeterminateAtSurface(Grid *);

Grid *gridWriteTecplotCurvedGeom(Grid *, char *filename );

Grid *gridShapeJacobian1(Grid *,
			 double *n0, double *n1, double *n2, double *n3,
			 double *where, double *jacobian );
Grid *gridShapeJacobian2(Grid *,
			 double *n0, double *n1, double *n2, double *n3,
			 double *e01, double *e02, double *e03,
			 double *e12, double *e13, double *e23,
			 double *where, double *jacobian );

END_C_DECLORATION

#endif /* GRIDSHAPE_H */
