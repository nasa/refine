
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

  tableau->constraint_matrix = (double *)malloc( tableauConstraints(tableau) *
						 tableauDimension( tableau ) * 
						 sizeof(double) );
  tableau->constraint        = (double *)malloc( tableauConstraints(tableau) *
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
  if (NULL != tableau->constraint_matrix) free(tableau->constraint_matrix);
  if (NULL != tableau->constraint)        free(tableau->constraint);
  if (NULL != tableau->cost)              free(tableau->cost);
  if (NULL != tableau->basis)             free(tableau->basis);
  if (NULL != tableau->t)                 free(tableau->t);
  free( tableau );
}

Tableau *tableauConstraintMatrix( Tableau *tableau, double *constraint_matrix )
{
  int i, length;
  length = tableauConstraints( tableau ) * tableauDimension( tableau );
  for (i=0;i<length;i++) tableau->constraint_matrix[i] = constraint_matrix[i];
  return tableau;
}

Tableau *tableauConstraint( Tableau *tableau, double *constraint )
{
  int i, length;
  length = tableauConstraints( tableau );
  for (i=0;i<length;i++) tableau->constraint[i] = constraint[i];
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

Tableau *tableauInit( Tableau *tableau )
{

  int i, j;
  int m, n;
  int t_index, A_index;
  double M;

  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );

  /* zero out current entries */
  for (i=0;i<m*n;i++) tableau->t[i] = 0.0 ;

  /* copy the A constraint matrix into tableu */
  for (i=0;i<tableauConstraints( tableau );i++) {
    for (j=0;j<tableauDimension( tableau );j++) {
      t_index = (1+i)+(1+j)*m;
      A_index = i+j*tableauDimension( tableau );
      tableau->t[t_index] = tableau->constraint_matrix[A_index];
    }
  }

  /* add identity slack variables */
  for (i=0;i<tableauConstraints( tableau );i++) {
    j = i+tableauDimension( tableau );
    t_index = (1+i)+(1+j)*m;
    tableau->t[t_index] = 1.0;
  }

  /* initial bfs */
  for (i=0;i<tableauConstraints(tableau);i++) 
    tableau->basis[i] = tableauDimension( tableau ) + i ;

  for (i=0;i<tableauConstraints( tableau );i++) {
    t_index = (1+i);
    tableau->t[t_index] = tableau->constraint[i];
  }
  
  /* find a large M value to start phase 1 */ 
  M = 0.0;
  for (j=0;j<tableauDimension( tableau );j++) {
    M = M + ABS(tableau->cost[j]);
  }
  
  /* intial cost */ 
  tableau->t[0] = 0.0;
  for (i=0;i<tableauConstraints( tableau );i++) {
    t_index = (1+i);
    tableau->t[0] -= M*tableau->t[t_index];
  }

  /* Reduced costs */ 
  for (j=0;j<tableauDimension( tableau );j++) {
    tableau->t[m*(1+j)] = tableau->cost[j];
    for (i=0;i<tableauConstraints( tableau );i++) {
      t_index = (1+i)+m*(1+j);
      tableau->t[m*(1+j)] -= M*tableau->t[t_index];
    }
  }

  return tableau;
}

Tableau *tableauSolve( Tableau *tableau )
{
  
  if ( tableau != tableauInit( tableau ) ) {
    printf( "%s: %d: %s: tableauInit NULL\n",
	    __FILE__, __LINE__, "tableauSolve");
    return NULL;
  }

  tableauShow(tableau);

  return tableau;
}

Tableau *tableauTableau( Tableau *tableau, double *tab )
{
  int i;
  int m, n;
  int length;

  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );
  length = m*n;

  for (i=0;i<length;i++) {
    tab[i] = tableau->t[i];
  }
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
