
/* Computes metrics from faces and tets 
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridinsert.h"
#include "gridmpi.h"

Grid *gridIdentityGlobal(Grid *grid, int offset )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodeGlobal(grid,node,node+offset)) return NULL;
  return grid;
}

Grid *gridSetAllLocal(Grid *grid )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodePart(grid,node,gridPartId(grid))) return NULL;
  return grid;
}

Grid *gridSetGhost(Grid *grid, int node )
{
  return gridSetNodePart(grid,node,EMPTY);
}

int gridParallelEdgeSplit(Grid *grid, Queue *queue, int node0, int node1 )
{
  double xyz0[3], xyz1[3];
  double newX, newY, newZ;
  int newnode;

  if ( gridNodeGhost(grid,node0) && gridNodeGhost(grid,node1) ) return EMPTY;

  if (grid != gridNodeXYZ(grid,node0,xyz0)) return EMPTY;
  if (grid != gridNodeXYZ(grid,node1,xyz1)) return EMPTY;

  newX = ( xyz0[0] + xyz1[0] ) * 0.5;
  newY = ( xyz0[1] + xyz1[1] ) * 0.5;
  newZ = ( xyz0[2] + xyz1[2] ) * 0.5;

  queueNewTransaction(queue);
  newnode = gridSplitEdgeAt( grid, queue, node0, node1, newX, newY, newZ );
  if (EMPTY == newnode) return EMPTY;
  gridSetNodePart(grid,newnode,gridPartId(grid));
  return newnode;
}
