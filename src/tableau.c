
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

#include <stdlib.h>
#include <stdio.h>
#ifndef __APPLE__       /* Not needed on Mac OS X */
#include <malloc.h>
#endif
#include "tableau.h"

Tableau* tableauCreate( int constraints, int dimension )
{
  int i;
  Tableau *tableau;

  tableau = (Tableau *)malloc( sizeof(Tableau) );

  tableau->constraints = constraints;
  tableau->dimension = dimension;

  tableau->basis = (int *)malloc( tableauConstraints(tableau) * sizeof(int) );
  for (i=0;i<tableauConstraints(tableau);i++) 
    tableau->basis[i] = tableauDimension( tableau ) + i ;

  return tableau;
}

void tableauFree( Tableau *tableau )
{
  if (NULL != tableau->basis) free(tableau->basis);
  free( tableau );
}

int tableauConstraints( Tableau *tableau )
{
  return tableau->constraints;
}

int tableauDimension( Tableau *tableau )
{
  return tableau->dimension;
}

Tableau *tableauBasis( Tableau *tableau, int *basis )
{
  int i;
  for (i=0;i<tableauConstraints(tableau);i++) basis[i] = tableau->basis[i];
  return tableau;
}
