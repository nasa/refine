
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

Grid *gridIdentityGlobal(Grid *grid )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodeGlobal(grid,node,node)) return NULL;
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

Grid *gridParallelEdgeSplit(Grid *grid, Queue *queue, int node1, int node2 )
{
  int newnode;
  if ( gridNodeGhost(grid,node1) && gridNodeGhost(grid,node2) ) return NULL;
  queueNewTransaction(queue);
  newnode = gridSplitEdge(grid,node1, node2);
  if (EMPTY == newnode) return NULL;
  gridSetNodePart(grid,newnode,gridPartId(grid));
  return grid;
}
