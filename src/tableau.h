
/* Tableau, an implementation of the tableau method for solving linear
 *   programming problems.
 *
 * Michael A. Park
 * Computational AeroSciences Branch
 * NASA Langley Research Center
 * Phone: (757) 864-6604
 * Email: Mike.Park@NASA.Gov
 */

/* $Id$ */

#ifndef TABLEAU_H
#define TABLEAU_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

typedef struct Tableau Tableau;

struct Tableau {
  int constraints;
  int dimension;
  
  double *constraint_matrix;
  double *constraint;
  double *cost;

  double *t;

  int *basis;
  int *in_basis;
};

Tableau *tableauCreate( int constraints, int dimension );
void tableauFree( Tableau * );

Tableau *tableauConstraintMatrix( Tableau *, double *constraint_matrix );
Tableau *tableauConstraint( Tableau *, double *constraint );
Tableau *tableauCost( Tableau *, double *cost );

int tableauConstraints( Tableau * );
int tableauDimension( Tableau * );

Tableau *tableauBasis( Tableau *, int *basis );

Tableau *tableauInit( Tableau * );

Tableau *tableauAuxillaryPivot( Tableau *, int *pivot_row, int *pivot_col );
Tableau *tableauLargestPivot( Tableau *, int *pivot_row, int *pivot_col );

Tableau *tableauSolve( Tableau * );
Tableau *tableauPivotAbout( Tableau *, int row, int column );

Tableau *tableauTableau( Tableau *, double *tableau );

Tableau *tableauShow( Tableau * );
Tableau *tableauShowTransposed( Tableau * );

END_C_DECLORATION

#endif /* RING_H */
