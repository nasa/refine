
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
  int order;
  double *f;
  Grid *grid;
};

Interp *interpCreate( Grid *grid, int function_id, int order );
void interpFree( Interp * );

#define interpFunctionId(interp) ((interp)->function_id)
#define interpOrder(interp) ((interp)->order)
#define interpGrid(interp) ((interp)->grid)

GridBool interpFunction( Interp *, double *xyz, double *func );
GridBool interpFunctionInCell( Interp *, int cell, double *bary, double *func );
GridBool interpMetric( Interp *, double *xyz, double *m );
GridBool interpError( Interp *, 
		      double *xyz0, double *xyz1, double *xyz2, double *xyz3, 
		      double *error );

GridBool interpTecplot( Interp *, char *filename );

double interpTotalError( Grid *grid );

END_C_DECLORATION

#endif /* INTERP_H */
