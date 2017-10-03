
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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
double tableauBound( Tableau * );

Tableau *tableauShow( Tableau * );
Tableau *tableauShowTransposed( Tableau * );

END_C_DECLORATION

#endif /* RING_H */
