
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
#include "gridmetric.h"
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

  gridMakeGem(grid, node0, node1 );
  if ( NULL == queue && !gridGemIsLocal(grid)) return EMPTY;
  if ( NULL != queue && gridGemIsLocal(grid)) return EMPTY;

  if (grid != gridNodeXYZ(grid,node0,xyz0)) return EMPTY;
  if (grid != gridNodeXYZ(grid,node1,xyz1)) return EMPTY;

  newX = ( xyz0[0] + xyz1[0] ) * 0.5;
  newY = ( xyz0[1] + xyz1[1] ) * 0.5;
  newZ = ( xyz0[2] + xyz1[2] ) * 0.5;

  if (NULL != queue) queueNewTransaction(queue);
  newnode = gridSplitEdgeAt( grid, queue, node0, node1, newX, newY, newZ );
  if (EMPTY == newnode) {
    printf("WARNING: %s: %d: global n cel counter may be hosed.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }
  
  return newnode;
}

Grid *gridApplyQueue(Grid *grid, Queue *gq )
{
  int transaction;
  int removed, removedcell, removedface;
  int i, globalnodes[5], localnodes[4];
  int cell, face;
  int added, addedcell, addedface;
  double xyz[12], uv[6];
  Queue *lq;

  lq = queueCreate();

  removedcell = 0;
  removedface = 0;
  addedcell = 0;
  addedface = 0;
  for (transaction=0;transaction<queueTransactions(gq);transaction++){
    for (removed=0;removed<queueRemovedCells(gq,transaction);removed++) {
      queueRemovedCellNodes( gq, removedcell, globalnodes );
      removedcell++;
      for(i=0;i<4;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
      cell = gridFindCell(grid,localnodes);
      if (grid == gridRemoveCell(grid,cell)) queueRemoveCell(lq,localnodes);
    }
    for (removed=0;removed<queueRemovedFaces(gq,transaction);removed++) {
      queueRemovedFaceNodes( gq, removedface, globalnodes );
      removedface++;
      for(i=0;i<3;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
      face = gridFindFace(grid,localnodes[0],localnodes[1],localnodes[2]);
      gridRemoveFace(grid,face);
    }

    for(added=0;added<queueAddedCells(gq,transaction);added++) {
      queueAddedCellNodes( gq, addedcell, globalnodes );
      queueAddedCellXYZs( gq, addedcell, xyz );
      addedcell++;
      for(i=0;i<4;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
      if ( gridNodeLocal(grid,localnodes[0]) ||
	   gridNodeLocal(grid,localnodes[1]) ||
	   gridNodeLocal(grid,localnodes[2]) ||
	   gridNodeLocal(grid,localnodes[3]) ) {
	for(i=0;i<4;i++) {
	  if ( EMPTY == localnodes[i] ) {
	    localnodes[i]=gridAddNode(grid,xyz[0+3*i],xyz[1+3*i],xyz[2+3*i]);
	    gridSetNodeGlobal(grid, localnodes[i], globalnodes[i]);
	    gridSetNodePart(grid, localnodes[i], EMPTY);
	  }
	}
	cell = gridAddCell(grid,
			   localnodes[0],localnodes[1],
			   localnodes[2],localnodes[3]);
	gridSetCellGlobal(grid,cell,globalnodes[4]);
      }
    }
    for(added=0;added<queueAddedFaces(gq,transaction);added++) {
      queueAddedFaceNodes( gq, addedface, globalnodes );
      queueAddedFaceUVs( gq, addedface, uv );
      addedface++;
      for(i=0;i<3;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
      localnodes[3] = globalnodes[3];
      if ( gridNodeLocal(grid,localnodes[0]) ||
	   gridNodeLocal(grid,localnodes[1]) ||
	   gridNodeLocal(grid,localnodes[2]) ) {
	for(i=0;i<3;i++)
	  if ( EMPTY == localnodes[i] ) printf("face insert error.");
	face = gridAddFaceUV(grid,
			     localnodes[0],uv[0],uv[1],
			     localnodes[1],uv[2],uv[3],
			     localnodes[2],uv[4],uv[5],
			     globalnodes[3]);
      }
    }
    queueNewTransaction(lq);
  }
  for (removed=0;removed<queueTotalRemovedCells(lq);removed++) {
    queueRemovedCellNodes( lq, removed, localnodes );
    for (i=0;i<4;i++) {
      if ( 0 == gridCellDegree(grid,localnodes[i]) ) {
	gridRemoveNode(grid,localnodes[i]);
      }
    }
  }
  

  queueFree(lq);
  return grid;
}

Grid *gridParallelAdaptWithOutCAD(Grid *grid, Queue *queue,
				  double minLength, double maxLength )
{
  int n0, n1, adaptnode, origNNode, newnode;
  int report, nnodeAdd, nnodeRemove;
  double ratio;

  origNNode = gridNNode(grid);
  adaptnode =0;
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (gridNNode(grid) > 100) report = gridNNode(grid)/10;

  for ( n0=0; adaptnode<origNNode; n0++ ) { 
    adaptnode++;
    if (adaptnode > 100 &&adaptnode/report*report == adaptnode )
      printf(" %6d adapt node %8d nnode %8d added %8d AR %15.10f\n",
	     gridPartId(grid),
	     adaptnode,gridNNode(grid),nnodeAdd,gridMinAR(grid));
    if ( gridValidNode( grid, n0) && 
	 !gridNodeFrozen( grid, n0 ) && 
	 gridNodeLocal( grid, n0 ) ) {
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      if ( !gridNodeFrozen( grid, n1 ) && ratio > maxLength ) {
	newnode = gridParallelEdgeSplit(grid, queue, n0, n1);
	if ( newnode != EMPTY ){
	  nnodeAdd++;
	}
      }
    }else{
      adaptnode++;
    }
  }
  return grid;
}
