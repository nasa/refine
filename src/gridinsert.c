
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <values.h>
#include "adj.h"
#include "gridStruct.h"
#include "gridmetric.h"
#include "gridinsert.h"

Grid *gridThrash(Grid *grid)
{
  int ncell, nodelimit, cellId, nodes[4];
  ncell = grid->ncell;
  nodelimit = grid->nnode*3/2;
  for (cellId=0;cellId<ncell && grid->nnode<nodelimit;cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      gridSplitEdge( grid, nodes[0], nodes[1] );
  
  return grid;
}

Grid *gridAdapt(Grid *grid)
{
  AdjIterator it;
  int i, n0, n1, adaptnode;
  double ratio;
  
  adaptnode = grid->nnode;

  for ( n0=0; n0<adaptnode ; n0++ ) { 
    if ( gridValidNode( grid, n0) ) {
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      // printf("n0 %d n1 %d ratio %f \n",n0,n1,ratio);
      if ( ratio > 2.0 ) gridSplitEdge(grid, n0, n1);
    }else{
      adaptnode++;
    }
  }
  return grid;
}

Grid *gridSplitEdge(Grid *grid, int n0, int n1 )
{
  int  igem, cell, nodes[4], inode, node, newnode, newnodes0[4], newnodes1[4];
  double newX, newY, newZ;
  int gap0, gap1, face0, face1, faceId0, faceId1;
  int edge, edgeId;
  double t0,t1, newT;

  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;

  newX = ( grid->xyz[0+3*n0] + grid->xyz[0+3*n1] ) * 0.5;
  newY = ( grid->xyz[1+3*n0] + grid->xyz[1+3*n1] ) * 0.5;
  newZ = ( grid->xyz[2+3*n0] + grid->xyz[2+3*n1] ) * 0.5;
  newnode = gridAddNode(grid, newX, newY, newZ );
  if ( newnode == EMPTY ) return NULL;
  grid->spacing[newnode] = 0.5*(grid->spacing[n0]+grid->spacing[n1]);

  for ( igem=0 ; igem<grid->ngem ; igem++ ){
    cell = grid->gem[igem];
    gridCell(grid, cell, nodes);
    gridRemoveCell(grid, cell);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      node = nodes[inode];
      newnodes0[inode]=node;
      newnodes1[inode]=node;
      if ( node == n0 ) newnodes0[inode] = newnode;
      if ( node == n1 ) newnodes1[inode] = newnode;
    }
    gridAddCell(grid, newnodes0[0], newnodes0[1], newnodes0[2], newnodes0[3] );
    gridAddCell(grid, newnodes1[0], newnodes1[1], newnodes1[2], newnodes1[3] );
    
  }

  //test face
  if ( grid->nequ != grid->ngem ){
    double n0Id0uv[2], n1Id0uv[2], n0Id1uv[2], n1Id1uv[2];
    double gap0uv[2], gap1uv[2], newId0uv[2], newId1uv[2]; 
    gap0 = grid->equ[0];
    gap1 = grid->equ[grid->ngem];
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    faceId0 = gridFaceId(grid, n0, n1, gap0 );
    faceId1 = gridFaceId(grid, n0, n1, gap1 );
    gridNodeUV(grid,n0,faceId0,n0Id0uv);
    gridNodeUV(grid,n1,faceId0,n1Id0uv);
    gridNodeUV(grid,n0,faceId1,n0Id1uv);
    gridNodeUV(grid,n1,faceId1,n1Id1uv);
    gridNodeUV(grid,gap0,faceId0,gap0uv);
    gridNodeUV(grid,gap1,faceId1,gap1uv);
    newId0uv[0] = 0.5 * (n0Id0uv[0]+n1Id0uv[0]);
    newId0uv[1] = 0.5 * (n0Id0uv[1]+n1Id0uv[1]);
    newId1uv[0] = 0.5 * (n0Id1uv[0]+n1Id1uv[0]);
    newId1uv[1] = 0.5 * (n0Id1uv[1]+n1Id1uv[1]);

    if ( faceId0 == EMPTY || faceId1 == EMPTY ) return NULL;

    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );
    gridAddFaceUV(grid, 
		  n0, n0Id0uv[0], n0Id0uv[1],
		  newnode, newId0uv[0],newId0uv[1],
		  gap0, gap0uv[0], gap0uv[1],
		  faceId0 );
    gridAddFaceUV(grid, 
		  n1, n1Id0uv[0], n1Id0uv[1], 
		  gap0, gap0uv[0], gap0uv[1], 
		  newnode, newId0uv[0], newId0uv[1], 
		  faceId0 );
    gridAddFaceUV(grid, 
		  n0, n0Id1uv[0], n0Id1uv[1], 
		  gap1, gap1uv[0], gap1uv[1], 
		  newnode, newId1uv[0], newId1uv[1], 
		  faceId1 );
    gridAddFaceUV(grid, 
		  n1, n1Id1uv[0], n1Id1uv[1], 
		  newnode, newId1uv[0], newId1uv[1], 
		  gap1, gap1uv[0], gap1uv[1], 
		  faceId1 );
    edge = gridFindEdge(grid,n0,n1);
    if ( edge != EMPTY ) {
      edgeId = gridEdgeId(grid,n0,n1);
      t0 = grid->edgeT[0+2*edge];
      t1 = grid->edgeT[1+2*edge];
      newT = 0.5 * (t0+t1);
      gridRemoveEdge(grid,edge);
      gridAddEdge(grid,n0,newnode,edgeId,t0,newT);
      gridAddEdge(grid,n1,newnode,edgeId,t1,newT);
    }
  }

  return grid;
}

Grid *gridCollapseEdge(Grid *grid, int n0, int n1 )
{
  int i, cell, face, face0, face1, faceId;
  double xyz0[3], xyz1[3], xyzAvg[3];
  double uv0[2], uv1[2], uvAvg[2];
  AdjIterator it;

  if ( gridGeometryEdge(grid, n0) ) return NULL;
  if ( gridGeometryEdge(grid, n1) ) return NULL;
  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;
  gridNodeXYZ( grid, n0, xyz0);
  gridNodeXYZ( grid, n1, xyz1);
  
  for (i=0 ; i<3 ; i++) xyzAvg[i] = 0.5 * ( xyz0[i] + xyz1[i] );
  for (i=0 ; i<3 ; i++) grid->xyz[i+3*n0] = xyzAvg[i];
  for (i=0 ; i<3 ; i++) grid->xyz[i+3*n1] = xyzAvg[i];
  
  if ( gridNegCellAroundNodeExceptGem( grid, n0 ) || 
       gridNegCellAroundNodeExceptGem( grid, n1 ) ) {
    for (i=0 ; i<3 ; i++) grid->xyz[i+3*n0] = xyz0[i];
    for (i=0 ; i<3 ; i++) grid->xyz[i+3*n1] = xyz1[i];
    return NULL;
  }

  for (i=0 ; i<grid->ngem ; i++) gridRemoveCell( grid, grid->gem[i] );

  it = adjFirst(grid->cellAdj, n1);
  while (adjValid(it)) {
    cell = adjItem(it);
    for ( i=0 ; i<4 ; i++ ) 
      if (grid->c2n[i+4*cell] == n1 ) 
	grid->c2n[i+4*cell] = n0;
    adjRemove( grid->cellAdj, n1, cell );
    adjRegister( grid->cellAdj, n0, cell );
    it = adjNext(it);
  }

  if ( grid->nequ != grid->ngem ) {
    faceId = grid->faceId[adjItem(adjFirst(grid->faceAdj,n0))];
    gridNodeUV(grid,n0,faceId,uv0);
    gridNodeUV(grid,n1,faceId,uv1);
    for (i=0 ; i<2 ; i++) uvAvg[i] = 0.5 * ( uv0[i] + uv1[i] );

    face0 = gridFindFace(grid, n0, n1, grid->equ[0] );
    face1 = gridFindFace(grid, n0, n1, grid->equ[grid->ngem] );
    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );

    it = adjFirst(grid->faceAdj, n1);
    while (adjValid(it)) {
      face = adjItem(it);
      for ( i=0 ; i<3 ; i++ ) 
	if (grid->f2n[i+3*face] == n1 ) 
	  grid->f2n[i+3*face] = n0;
      adjRemove( grid->faceAdj, n1, face );
      adjRegister( grid->faceAdj, n0, face );
      it = adjNext(it);
    }
    gridSetNodeUV(grid, n0, faceId, uvAvg[0], uvAvg[1]);
  }

  gridRemoveNode(grid, n1);

  return grid;
}

