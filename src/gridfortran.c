
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridfortran.h"
#include "grid.h"

static Grid *grid;

int gridcreate_( int *nnode, double *x, double *y, double *z )
{
  int node;
  grid = gridCreate( *nnode, 5000, 5000, 0);
  for ( node=0; node<*nnode; node++) gridAddNode(grid,x[node],y[node],z[node]);
  printf("populated grid object with %d nodes\n",gridNNode(grid));
}


