
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

Grid *gridIdentityNodeGlobal(Grid *grid, int offset )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodeGlobal(grid,node,node+offset)) return NULL;
  return grid;
}

Grid *gridIdentityCellGlobal(Grid *grid, int offset )
{
  int cell;
  for (cell = 0; cell < gridNCell(grid) ; cell++ )
    if (grid != gridSetCellGlobal(grid,cell,cell+offset)) return NULL;
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
  if (EMPTY == newnode) {
    printf("WARNING: %s: %d: global n cel counter may be hosed.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }
  gridSetNodePart(grid,newnode,gridPartId(grid));
  
  return newnode;
}

Grid *gridApplyQueue(Grid *grid, Queue *queue )
{
  int transaction;
  int removed, removedcell;
  int i, globalnodes[5], localnodes[4];
  int cell;

  removedcell = 0;
  for (transaction=0;transaction<queueTransactions(queue);transaction++){
    printf( "transaction %d has %d removed cells\n",
	    transaction,queueRemovedCells(queue,transaction));
    for (removed=0;removed<queueRemovedCells(queue,transaction);removed++) {
      queueRemovedCellNodes( queue, removedcell, globalnodes );
      printf( "global cell %d %d %d %d\n",
	      globalnodes[0],
	      globalnodes[1],
	      globalnodes[2],
	      globalnodes[3]);
      removedcell++;
      for(i=0;i<4;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
      printf( "local cell %d %d %d %d\n",
	      localnodes[0],
	      localnodes[1],
	      localnodes[2],
	      localnodes[3]);
      cell = gridFindCell(grid,localnodes);
      gridRemoveCell(grid,cell);
    }
  }

  return grid;
}
