
/* Interp, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id: interp.h 2867 2005-07-18 17:32:34Z mikepark $ */

#ifndef INTERP_H
#define INTERP_H

#include <stdio.h>
#include "refine_defs.h"

BEGIN_C_DECLORATION

typedef struct Interp Interp;

#include "grid.h"

struct Interp {
  int function_id;
  int error_order;
  int order;
  int dim;
  double *w;
  double *f;
  Grid *grid;
};

Interp *interpCreate( Grid *grid, int function_id, int order, int error_order );
Interp *interpDirect( Grid *grid );

void interpFree( Interp * );

Interp *interpContinuousReconstruction( Interp *, int order, int error_order );

#define interpFunctionId(interp) ((interp)->function_id)
#define interpErrorOrder(interp) ((interp)->error_order)
#define interpOrder(interp) ((interp)->order)
#define interpDim(interp) ((interp)->dim)
#define interpGrid(interp) ((interp)->grid)

int interpNB(Interp *);
GridBool interpPhi( Interp *, double *bary, double *phi );
GridBool interpLoc2Row( Interp *, int cell, int *loc2row );

GridBool interpFunction( Interp *, double *xyz, 
			 double *func, double *weight );
GridBool interpFunctionInCell( Interp *, int cell, double *bary, 
			       double *func, double *weight );
GridBool interpMetric( Interp *, double *xyz, double *m );
GridBool interpError( Interp *, 
		      double *xyz0, double *xyz1, double *xyz2, double *xyz3, 
		      double *error );

GridBool interpError1D( Interp *interp,
			double *xyz0, double *xyz1, 
			double *error );

GridBool interpSplitImprovement1D( Interp *, 
				   double *xyz0, double *xyz1,
				   double *error_before, double *error_after );

GridBool interpSplitImprovement( Interp *, 
				 double *xyz0, double *xyz1,
				 double *xyz2, double *xyz3,
				 double *error_before, double *error_after );

GridBool interpTecplot( Interp *, char *filename );

double interpTotalError( Grid *grid );

END_C_DECLORATION

#endif /* INTERP_H */
