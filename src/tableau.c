
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
#include <math.h>
#include <limits.h>
#include <values.h>
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
  length = 1+tableauDimension( tableau )+tableauConstraints(tableau);
  tableau->in_basis = (int *)malloc( length * sizeof(int) );
  for (i=0;i<length;i++)  tableau->in_basis[i] = EMPTY;
  for (i=0;i<tableauConstraints(tableau);i++) 
    tableau->in_basis[tableau->basis[i]+1] = i;
  

  return tableau;
}

void tableauFree( Tableau *tableau )
{
  if (NULL != tableau->constraint_matrix) free(tableau->constraint_matrix);
  if (NULL != tableau->constraint)        free(tableau->constraint);
  if (NULL != tableau->cost)              free(tableau->cost);
  if (NULL != tableau->basis)             free(tableau->basis);
  if (NULL != tableau->in_basis)          free(tableau->in_basis);
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
  for (j=0;j<n;j++) tableau->in_basis[j] = EMPTY;
  for (i=0;i<tableauConstraints(tableau);i++) 
    tableau->in_basis[tableau->basis[i]+1] = i;

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

Tableau *tableauLargestPivot( Tableau *tableau, int *pivot_row, int *pivot_col )
{
  int i,j;
  int m,n;

  double zero;

  double divisor, best_divisor;
  double reduced_cost;
  
  int best_row;
  double feasable_step_length, this_step_length;
  double pivot;

  *pivot_row = EMPTY;
  *pivot_col = EMPTY;

  zero = 1.0e-14;
  
  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );
  
  best_divisor = 0.0;

  for (j=1;j<n;j++) { /* test all basis, first col is solution */
    if (EMPTY == tableau->in_basis[j]) { /* skip active basis */
      reduced_cost = tableau->t[m*j];
      if ( 0 > reduced_cost ) {  /* try a negative reduced cost */
	best_row = EMPTY;
	feasable_step_length = DBL_MAX;
	for (i=1;i<m;i++) {
	  pivot = tableau->t[i+m*j];
	  if ( pivot > zero ) {
	    this_step_length = reduced_cost / pivot;
	    if ( this_step_length < feasable_step_length ) {
	      best_row = i;
	      feasable_step_length = this_step_length;
	      divisor = pivot;
	    }
	  }
	} /* end loop over rows */
	/* if this column has the best pivot so far, keep it */
	if ( EMPTY != best_row && ABS(divisor) > ABS(best_divisor) ) {
	  *pivot_row = best_row;
	  *pivot_col = j;
	  best_divisor = divisor;
	}
      }
    }
  }/* end loop over columns */
  if ( EMPTY == (*pivot_row) ) return NULL;
  return tableau;
}

Tableau *tableauSolve( Tableau *tableau )
{
  int row, column;
  if ( tableau != tableauInit( tableau ) ) {
    printf( "%s: %d: %s: tableauInit NULL\n",
	    __FILE__, __LINE__, "tableauSolve");
    return NULL;
  }

  while ( tableau == tableauLargestPivot( tableau, &row, &column ) ) {
    tableauPivotAbout(tableau, row, column);
  }

  return tableau;
}

Tableau *tableauPivotAbout( Tableau *tableau, int row, int column )
{
  int m, n;
  int i, j;
  double pivot;
  double factor;

  m = 1 + tableauConstraints( tableau );
  n = 1 + tableauDimension( tableau ) + tableauConstraints( tableau );

  if ( row < 1 || row >= m ||  column < 1 || column >= n ) {
    printf( "%s: %d: %s: requested row %d or column %d is outside %d m %d n\n",
	    __FILE__, __LINE__, "tableauPivotAbout",row,column,m,n);
    return NULL;
  }

  if (EMPTY != tableau->in_basis[column]) {
    printf( "%s: %d: %s: requested column %d is already active\n",
	    __FILE__, __LINE__, "tableauPivotAbout",column);
    return NULL;
  }

  tableau->in_basis[tableau->basis[row-1]+1] = EMPTY;
  tableau->basis[row-1] = column-1;
  tableau->in_basis[tableau->basis[row-1]+1] = row-1;

  /* normalize row */
  pivot = tableau->t[row+m*column];
  for (j=0;j<n;j++) {
    tableau->t[row+m*j] /= pivot;
  }

  for (i=0;i<m;i++) {
    if (i!=row) {
      factor = tableau->t[i+m*column];
      for (j=0;j<n;j++) {
	tableau->t[i+m*j] -= tableau->t[row+m*j]*factor;
      }
    }
  }

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
