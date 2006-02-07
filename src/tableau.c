
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

  tableau->constriant_matrix = (double *)malloc( tableauConstraints(tableau) *
						 tableauDimension( tableau ) * 
						 sizeof(double) );
  tableau->constriant        = (double *)malloc( tableauConstraints(tableau) *
						 sizeof(double) );
  tableau->cost              = (double *)malloc( tableauDimension(tableau) *
						 sizeof(double) );

  return tableau;
}

void tableauFree( Tableau *tableau )
{
  if (NULL != tableau->constriant_matrix) free(tableau->constriant_matrix);
  if (NULL != tableau->constriant) free(tableau->constriant);
  if (NULL != tableau->cost) free(tableau->cost);
  if (NULL != tableau->basis) free(tableau->basis);
  free( tableau );
}

Tableau *tableauConstraintMatrix( Tableau *tableau, double *constriant_matrix )
{
  int i, length;
  length = tableauConstraints( tableau ) * tableauDimension( tableau );
  for (i=0;i<length;i++) tableau->constriant_matrix[i] = constriant_matrix[i];
  return tableau;
}

Tableau *tableauConstraint( Tableau *tableau, double *constriant )
{
  int i, length;
  length = tableauConstraints( tableau );
  for (i=0;i<length;i++) tableau->constriant[i] = constriant[i];
  return tableau;
}

Tableau *tableauCost( Tableau *tableau, double *cost )
{
  int i, length;
  length = tableauDimension( tableau );
  for (i=0;i<length;i++) tableau->cost[i] = cost[i];
  return tableau;
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
