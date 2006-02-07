
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
  int i, length;
  Tableau *tableau;

  tableau = (Tableau *)malloc( sizeof(Tableau) );

  tableau->constraints = constraints;
  tableau->dimension = dimension;

  tableau->constriant_matrix = (double *)malloc( tableauConstraints(tableau) *
						 tableauDimension( tableau ) * 
						 sizeof(double) );
  tableau->constriant        = (double *)malloc( tableauConstraints(tableau) *
						 sizeof(double) );
  tableau->cost              = (double *)malloc( tableauDimension(tableau) *
						 sizeof(double) );

  length = (1+tableauConstraints(tableau)) * 
           (1+tableauDimension( tableau )+tableauConstraints(tableau));
  tableau->t = (double *)malloc(  length * sizeof(double) );
  for (i=0;i<length;i++) tableau->t[i] = 0.0 ;

  tableau->basis = (int *)malloc( tableauConstraints(tableau) * sizeof(int) );
  for (i=0;i<tableauConstraints(tableau);i++) 
    tableau->basis[i] = tableauDimension( tableau ) + i ;

  return tableau;
}

void tableauFree( Tableau *tableau )
{
  if (NULL != tableau->constriant_matrix) free(tableau->constriant_matrix);
  if (NULL != tableau->constriant)        free(tableau->constriant);
  if (NULL != tableau->cost)              free(tableau->cost);
  if (NULL != tableau->basis)             free(tableau->basis);
  if (NULL != tableau->t)                 free(tableau->t);
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

Tableau *tableauSolve( Tableau *tableau )
{

  int i, j;
  int m, n;
  int t_index, A_index;

  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );

  for (i=0;i<tableauConstraints( tableau );i++) {
    for (j=0;j<tableauDimension( tableau );j++) {
      t_index = (1+i)+(1+j)*m;
      A_index = i+j*tableauDimension( tableau );
      tableau->t[t_index] = tableau->constriant_matrix[A_index];
    }
  }

  tableauShow(tableau);

  return tableau;
}

Tableau *tableauShow( Tableau *tableau )
{
  int i, j;
  int m, n;

  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );

  printf("\n");
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      printf(" %8.4f",tableau->t[i+j*m]);
    }
    printf("\n");
  }
  return tableau;
}
