
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
#include "gridmetric.h"
#include "gridinsert.h"
#include "gridcad.h"

Grid *gridThrash(Grid *grid)
{
  int ncell, nodelimit, cellId, nodes[4];
  ncell = gridNCell(grid);
  nodelimit = gridNNode(grid)*3/2;
  for (cellId=0;cellId<ncell && gridNNode(grid)<nodelimit;cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      gridSplitEdge( grid, nodes[0], nodes[1] );
  
  return grid;
}

Grid *gridRemoveAllNodes(Grid *grid )
{
  AdjIterator it;
  int n0, i, cell, nodes[4];
  bool nodeExists;
  double ratio;
 
  for ( n0=0; n0<gridMaxNode(grid); n0++ ) { 
    if ( gridValidNode(grid, n0) && !gridNodeFrozen(grid, n0) ) {
      nodeExists = TRUE;
      for ( it = adjFirst(gridCellAdj(grid),n0); 
	    nodeExists && adjValid(it); 
	    it = adjNext(it) ){
	cell = adjItem(it);
	gridCell( grid, cell, nodes);
	for (i=0;nodeExists && i<4;i++){
	  if (n0 != nodes[i]) {
	    ratio = 1.0;
	    if (!gridNodeFrozen(grid, nodes[i]))ratio=0.5;
	    nodeExists = (grid == gridCollapseEdge(grid, nodes[i], n0, 1.0));
	  }
	}
      }
    }
  }
  return grid;
}

Grid *gridAdapt(Grid *grid, double minLength, double maxLength )
{
  AdjIterator it;
  int i, n0, n1, adaptnode, origNNode, newnode;
  int report, nnodeAdd, nnodeRemove;
  double ratio;
  
  origNNode = gridNNode(grid);
  adaptnode =0;
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (gridNNode(grid) > 100) report = gridNNode(grid)/10;

  for ( n0=0; 
	adaptnode<origNNode && n0<gridMaxNode(grid); 
	n0++ ) { 
    if (adaptnode > 100 &&adaptnode/report*report == adaptnode )
      printf("adapt node %8d nnode %8d added %8d removed %8d\n",
	     adaptnode,gridNNode(grid),nnodeAdd,nnodeRemove);
    if ( gridValidNode( grid, n0) && !gridNodeFrozen( grid, n0 ) ) {
      adaptnode++;
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      if ( !gridNodeFrozen( grid, n1 ) && ratio > maxLength ) {
	newnode = gridSplitEdge(grid, n0, n1);
	if ( newnode != EMPTY ){
	  nnodeAdd++;
	  gridSwapNearNode( grid, newnode );
	  if (gridGeometryFace( grid, newnode ) ){
	    gridRobustProjectNode(grid, newnode);
	    gridSwapNearNode( grid, newnode );
	  }
	}
      }else{
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( !gridNodeFrozen( grid, n1 ) && ratio < minLength ) { 
	  if ( grid == gridCollapseEdge(grid, n0, n1, 0.5) ) {
	    nnodeRemove++;
	    gridSwapNearNode( grid, n0 );
	    if (  gridGeometryFace( grid, n0 ) ) {
	      gridRobustProjectNode(grid, n0);
	      gridSwapNearNode( grid, n0 );
	    }
	  }
	}
      }
    }else{
      adaptnode++;
    }
  }
  return grid;
}

Grid *gridAdaptSurface(Grid *grid, double minLength, double maxLength )
{
  AdjIterator it;
  int i, n0, n1, adaptnode, origNNode, newnode;
  int report, nnodeAdd, nnodeRemove;
  double ratio;
  
  origNNode = gridNNode(grid);
  adaptnode =0;
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (gridNNode(grid) > 100) report = gridNNode(grid)/10;

  for ( n0=0; 
	adaptnode<origNNode && n0<gridMaxNode(grid); 
	n0++ ) { 
    if (adaptnode > 100 &&adaptnode/report*report == adaptnode )
      printf("adapt node %8d nnode %8d added %8d removed %8d\n",
	     adaptnode,gridNNode(grid),nnodeAdd,nnodeRemove);
    if ( gridValidNode( grid, n0) && 
	 !gridNodeFrozen( grid, n0 ) && 
	 gridGeometryFace( grid, n0)) {
      adaptnode++;
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      if ( !gridNodeFrozen( grid, n1 ) &&
	   gridGeometryFace( grid, n1) &&
	   ratio > maxLength ) {
	newnode = gridSplitEdge(grid, n0, n1);
	if ( newnode != EMPTY ){
	  nnodeAdd++;
	  gridSwapNearNode( grid, newnode );
	  if (gridGeometryFace( grid, newnode ) ){
	    gridRobustProjectNode(grid, newnode);
	    gridSwapNearNode( grid, newnode );
	  }
	}
      }else{
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( !gridNodeFrozen( grid, n1 ) && ratio < minLength ) { 
	  if ( grid == gridCollapseEdge(grid, n0, n1, 0.5) ) {
	    nnodeRemove++;
	    gridSwapNearNode( grid, n0 );
	    if (  gridGeometryFace( grid, n0 ) ) {
	      gridRobustProjectNode(grid, n0);
	      gridSwapNearNode( grid, n0 );
	    }
	  }
	}
      }
    }else{
      adaptnode++;
    }
  }
  return grid;
}

Grid *gridAdaptWithOutCAD(Grid *grid, double minLength, double maxLength )
{
  AdjIterator it;
  int i, n0, n1, adaptnode, origNNode, newnode;
  int report, nnodeAdd, nnodeRemove;
  double ratio;
  
  origNNode = gridNNode(grid);
  adaptnode =0;
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (gridNNode(grid) > 100) report = gridNNode(grid)/10;

  for ( n0=0; 
	adaptnode<origNNode && n0<gridMaxNode(grid); 
	n0++ ) { 
    if (adaptnode > 100 &&adaptnode/report*report == adaptnode )
      printf("adapt node %8d nnode %8d added %8d removed %8d\n",
	     adaptnode,gridNNode(grid),nnodeAdd,nnodeRemove);
    if ( gridValidNode( grid, n0) && !gridNodeFrozen( grid, n0 ) ) {
      adaptnode++;
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      if ( !gridNodeFrozen( grid, n1 ) && ratio > maxLength ) {
	newnode = gridSplitEdge(grid, n0, n1);
	if ( newnode != EMPTY ){
	  nnodeAdd++;
	  gridSwapNearNode( grid, newnode );
	}
      }else{
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( !gridNodeFrozen( grid, n1 ) && ratio < minLength ) { 
	  if ( grid == gridCollapseEdge(grid, n0, n1, 0.5) ) {
	    nnodeRemove++;
	    gridSwapNearNode( grid, n0 );
	  }
	}
      }
    }else{
      adaptnode++;
    }
  }
  return grid;
}

int gridSplitEdge(Grid *grid, int n0, int n1)
{
  double xyz0[3], xyz1[3];
  double newX, newY, newZ;

  if (grid != gridNodeXYZ(grid,n0,xyz0)) return EMPTY;
  if (grid != gridNodeXYZ(grid,n1,xyz1)) return EMPTY;

  newX = ( xyz0[0] + xyz1[0] ) * 0.5;
  newY = ( xyz0[1] + xyz1[1] ) * 0.5;
  newZ = ( xyz0[2] + xyz1[2] ) * 0.5;
  
  return gridSplitEdgeAt(grid, NULL, n0, n1,
			 newX, newY, newZ );

}

int gridSplitEdgeAt(Grid *grid, Queue *queue, int n0, int n1,
		    double newX, double newY, double newZ )
{
  int i, igem, cell, nodes[4], inode, node;
  int newnode, newnodes0[4], newnodes1[4];
  int gap0, gap1, face0, face1, faceId0, faceId1;
  int edge, edgeId;
  double t0,t1, newT;
  double map, map0, map1;

  if ( NULL == gridEquator( grid, n0, n1) ) return EMPTY;

  newnode = gridAddNode(grid, newX, newY, newZ );
  if ( newnode == EMPTY ) return EMPTY;
  gridSetMapMatrixToAverageOfNodes(grid, newnode, n0, n1 );

  for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
    cell = gridGem(grid,igem);
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
  if ( !gridContinuousEquator(grid) ){
    double n0Id0uv[2], n1Id0uv[2], n0Id1uv[2], n1Id1uv[2];
    double gap0uv[2], gap1uv[2], newId0uv[2], newId1uv[2]; 
    gap0 = gridEqu(grid,0);
    gap1 = gridEqu(grid,gridNGem(grid));
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

    if ( faceId0 == EMPTY || faceId1 == EMPTY ) {
      return EMPTY;
    }

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
  }

  edge = gridFindEdge(grid,n0,n1);
  if ( edge != EMPTY ) {
    edgeId = gridEdgeId(grid,n0,n1);
    gridNodeT(grid,n0,edgeId,&t0);
    gridNodeT(grid,n1,edgeId,&t1);
    newT = 0.5 * (t0+t1);
    gridRemoveEdge(grid,edge);
    gridAddEdge(grid,n0,newnode,edgeId,t0,newT);
    gridAddEdge(grid,n1,newnode,edgeId,t1,newT);
  }
  
  if ( gridNegCellAroundNode(grid, newnode) ) {
    gridCollapseEdge(grid, n0, newnode, 0.0 );
    return EMPTY;
  }else{
    return newnode;
  }
}

int gridSplitEdgeIfNear(Grid *grid, int n0, int n1,
			double newX, double newY, double newZ)
{
  int i;
  int newnode;
  double oldXYZ[3], newXYZ[3], xyz0[3], xyz1[3];
  double edgeXYZ[3], edgeLength, edgeDir[3];
  double newEdge[3], edgePosition, radius, radiusVector[3];

  newXYZ[0] = newX;  newXYZ[1] = newY;  newXYZ[2] = newZ;

  gridNodeXYZ(grid, n0, xyz0);
  gridNodeXYZ(grid, n1, xyz1);
  for (i=0;i<3;i++) edgeXYZ[i] = xyz1[i]   - xyz0[i];
  edgeLength = sqrt ( edgeXYZ[0]*edgeXYZ[0] + 
		      edgeXYZ[1]*edgeXYZ[1] +
		      edgeXYZ[2]*edgeXYZ[2] );
  for (i=0;i<3;i++) edgeDir[i] = edgeXYZ[i]/edgeLength;
  for (i=0;i<3;i++) newEdge[i] = newXYZ[i] - xyz0[i];
  edgePosition =  ( newEdge[0]*edgeDir[0] + 
		    newEdge[1]*edgeDir[1] + 
		    newEdge[2]*edgeDir[2] ) / edgeLength;
  for (i=0;i<3;i++) radiusVector[i] = newEdge[i] -
		      edgePosition*edgeLength*edgeDir[i];
  radius = sqrt ( radiusVector[0]*radiusVector[0] + 
		  radiusVector[1]*radiusVector[1] + 
		  radiusVector[2]*radiusVector[2] );
      
  if ( edgePosition > -0.05 && edgePosition < 0.05 && 
       radius < 0.05*edgeLength) {
    newnode = n0;
    gridNodeXYZ(grid, newnode, oldXYZ );
    gridSetNodeXYZ(grid, newnode, newXYZ );
    if ( gridNegCellAroundNode(grid, newnode ) ){
      gridSetNodeXYZ(grid, newnode, oldXYZ );
    }else{
      return newnode;
    }
  }

  if ( edgePosition > 0.95 && edgePosition < 1.05 && 
       radius < 0.05*edgeLength) {
    newnode = n1;
    gridNodeXYZ(grid, newnode, oldXYZ );
    gridSetNodeXYZ(grid, newnode, newXYZ );
    if ( gridNegCellAroundNode(grid, newnode ) ){
      gridSetNodeXYZ(grid, newnode, oldXYZ );
    }else{
      return newnode;
    }
  }

  if ( edgePosition > 0.0 && edgePosition < 1.0 && 
       radius < 0.001*edgeLength) {
    newnode = gridSplitEdgeAt(grid, NULL, n0, n1, newX, newY, newZ);
    return newnode;
  }

  return EMPTY;
}

int gridSplitFaceAt(Grid *grid, int face,
		    double newX, double newY, double newZ )
{
  int newnode;
  int nodes[4], newnodes[4], faceId, cell;
  double U[3], V[3], avgU, avgV, newU[3], newV[3];
  int n, i;

  cell = gridFindCellWithFace(grid, face );
  if ( EMPTY == cell ) return EMPTY;
  
  /* set up nodes so that nodes[0]-nodes[2] is face 
     and nodes[0]-nodes[3] is the cell */

  if (grid != gridCell(grid, cell, nodes) ) return EMPTY;
  nodes[3] = nodes[0] + nodes[1] + nodes[2] + nodes[3];
  if (grid != gridFace(grid, face, nodes, &faceId ) ) return EMPTY;
  nodes[3] = nodes[3] -  nodes[0] - nodes[1] - nodes[2];

  //if (0.0>=gridVolume(grid, nodes ) ) return EMPTY;

  for (i=0;i<3;i++) {
    U[i] = gridNodeU(grid, nodes[i], faceId);
    V[i] = gridNodeV(grid, nodes[i], faceId);
  }
  avgU = (U[0] + U[1] + U[2]) / 3;
  avgV = (V[0] + V[1] + V[2]) / 3;

  newnode = gridAddNode(grid, newX, newY, newZ );
  if ( newnode == EMPTY ) return EMPTY;

  if ( grid != gridRemoveCell(grid, cell ) ) return EMPTY;
  if ( grid != gridRemoveFace(grid, face ) ) return EMPTY;

  for (i=0;i<3;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes[n];
    newnodes[i] = newnode; 
    for (n=0;n<3;n++) newU[n] = U[n];
    for (n=0;n<3;n++) newV[n] = V[n];
    newU[i] = avgU;
    newV[i] = avgV;
    if (EMPTY == gridAddFaceUV(grid, 
			       newnodes[0], newU[0], newV[0], 
			       newnodes[1], newU[1], newV[1],  
			       newnodes[2], newU[2], newV[2],  
			       faceId ) ) return EMPTY;
  }

  if (grid!=gridProjectNodeToFace(grid, newnode, faceId )) return EMPTY;

  for (i=0;i<3;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes[n];
    newnodes[i] = newnode; 
    if (EMPTY == gridAddCell(grid, newnodes[0], newnodes[1], 
			           newnodes[2], newnodes[3] ) ) return EMPTY;
  }

  if ( gridNegCellAroundNode(grid, newnode) ) {
    gridCollapseEdge(grid, nodes[0], newnode, 0.0 );
    return EMPTY;
  }else{
    return newnode;
  }
}

int gridSplitCellAt(Grid *grid, int cell,
		    double newX, double newY, double newZ )
{
  int newnode;
  int nodes[4], newnodes[4];
  int n, i;

  if (grid != gridCell(grid, cell, nodes) ) return EMPTY;

  //if (0.0>=gridVolume(grid, nodes ) ) return EMPTY;

  newnode = gridAddNode(grid, newX, newY, newZ );
  if ( newnode == EMPTY ) return EMPTY;

  if ( grid != gridRemoveCell(grid, cell ) ) return EMPTY;

  for (i=0;i<4;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes[n];
    newnodes[i] = newnode; 
    if (EMPTY == gridAddCell(grid, newnodes[0], newnodes[1], 
		     	           newnodes[2], newnodes[3] ) ) return EMPTY;
  }

  if ( gridNegCellAroundNode(grid, newnode) ) {
    gridCollapseEdge(grid, nodes[0], newnode, 0.0 );
    return EMPTY;
  }else{
    return newnode;
  }
}

int gridInsertInToGeomEdge(Grid *grid, double newX, double newY, double newZ)
{
  int i, edge, maxedge, edgeId, nodes[2];
  int newnode;

  newnode = EMPTY;
  edge = 0;
  maxedge = gridMaxEdge(grid);
  while ( EMPTY == newnode && edge < maxedge ) {
    if (grid == gridEdge(grid, edge, nodes, &edgeId) ){
      newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],newX,newY,newZ);
    }
    edge++;
  }

  if ( newnode != EMPTY ) gridSafeProjectNode(grid, newnode, 1.0);

  return newnode;
}

int gridInsertInToGeomFace(Grid *grid, double newX, double newY, double newZ)
{
  int foundFace;
  int i, face, maxface, faceId, nodes[3];
  double newxyz[3], xyz0[3], xyz1[3], xyz2[3];
  double edge0[3], edge1[3], edge2[3];
  double leg0[3], leg1[3], leg2[3];
  double norm[3], norm0[3], norm1[3], norm2[3];
  double normLength, unit[3], normDistance;
  int newnode;

  bool edgeSplit;

  newxyz[0] = newX;  newxyz[1] = newY;  newxyz[2] = newZ;

  foundFace = EMPTY;
  edgeSplit = FALSE;
  face = 0;
  maxface = gridMaxFace(grid);
  while ( EMPTY == foundFace && !edgeSplit && face < maxface ) {
    if (grid == gridFace(grid, face, nodes, &faceId) ){
      /* first try putting it on an edge */
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[2],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[2],nodes[0],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit) {
	/* then try putting it on the face */
	gridNodeXYZ(grid, nodes[0], xyz0);
	gridNodeXYZ(grid, nodes[1], xyz1);
	gridNodeXYZ(grid, nodes[2], xyz2);
	gridSubtractVector(xyz1, xyz0, edge0);
	gridSubtractVector(xyz2, xyz1, edge1);
	gridSubtractVector(xyz0, xyz2, edge2);
	gridSubtractVector(newxyz, xyz0, leg0);
	gridSubtractVector(newxyz, xyz1, leg1);
	gridSubtractVector(newxyz, xyz2, leg2);
	gridCrossProduct(edge0,edge1,norm);
	gridCrossProduct(edge0,leg0,norm0);
	gridCrossProduct(edge1,leg1,norm1);
	gridCrossProduct(edge2,leg2,norm2);
	normLength = sqrt( gridDotProduct( norm, norm ) );
	unit[0] = norm[0] / normLength;
	unit[1] = norm[1] / normLength;
	unit[2] = norm[2] / normLength;
	normDistance = gridDotProduct( unit, leg0 );
	if (FALSE) {
	  printf("f%d X%8.5f Y%8.5f Z%8.5f %6.3f 0 %6.3f 1 %6.3f 2 %6.3f\n",
		 face, newX, newY, newZ, normDistance,
		 gridDotProduct( norm, norm0 ),
		 gridDotProduct( norm, norm1 ),
		 gridDotProduct( norm, norm2 ) );
	}
	if ( gridDotProduct( norm, norm0 ) > 0.00 &&
	     gridDotProduct( norm, norm1 ) > 0.00 &&
	     gridDotProduct( norm, norm2 ) > 0.00 &&
	     ABS(normDistance) < 0.01*normLength ){
	  foundFace = face;
	}
      }
    }
    face++;
  }

  if ( edgeSplit ) return newnode;
  if ( EMPTY == foundFace ) return EMPTY;

  newnode = gridSplitFaceAt(grid, foundFace, newX, newY, newZ);
  if (EMPTY == newnode) return EMPTY;

  return newnode;
}

int gridInsertInToVolume(Grid *grid, double newX, double newY, double newZ)
{
  int cell, maxcell, nodes[4], foundCell, newnode;
  bool edgeSplit;
  double newxyz[3], xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  double edge0[3], edge1[3];
  double leg0[3], leg1[3];
  double norm[3], norm0[3], norm1[3], norm2[3], norm3[3];

  newxyz[0] = newX;  newxyz[1] = newY;  newxyz[2] = newZ;

  foundCell = EMPTY;
  edgeSplit = FALSE;
  cell = 0;
  maxcell = gridMaxCell(grid);
  while ( EMPTY == foundCell && !edgeSplit && cell < maxcell ) {
    if (grid == gridCell(grid, cell, nodes) ){
      /* first try putting it on an edge */
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[2],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[3],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[2],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[3],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[2],nodes[3],newX,newY,newZ);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit) {
	/* then try putting it in the cell -or splitting faces?-*/
	gridNodeXYZ(grid, nodes[0], xyz0);
	gridNodeXYZ(grid, nodes[1], xyz1);
	gridNodeXYZ(grid, nodes[2], xyz2);
	gridNodeXYZ(grid, nodes[3], xyz3);

	gridSubtractVector(xyz3, xyz1, edge0);
	gridSubtractVector(xyz2, xyz1, edge1);
	gridCrossProduct(edge0,edge1,norm0);

	gridSubtractVector(xyz2, xyz0, edge0);
	gridSubtractVector(xyz3, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm1);

	gridSubtractVector(xyz3, xyz0, edge0);
	gridSubtractVector(xyz1, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm2);

	gridSubtractVector(xyz1, xyz0, edge0);
	gridSubtractVector(xyz2, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm3);

	gridSubtractVector(newxyz, xyz0, leg0);
	gridSubtractVector(newxyz, xyz1, leg1);

	if (FALSE) {	
	  printf("c%d X%8.5f Y%8.5f Z%8.5f 0 %6.3f 1 %6.3f 2 %6.3f 3 %6.3f\n",
		 cell, newX, newY, newZ, 
		 gridDotProduct( leg1, norm0 ),
		 gridDotProduct( leg0, norm1 ),
		 gridDotProduct( leg0, norm2 ),
		 gridDotProduct( leg0, norm3 ) );
	  printf("xyz0 X%8.5f Y%8.5f Z%8.5f\n",xyz0[0],xyz0[1],xyz0[2]);
	  printf("xyz1 X%8.5f Y%8.5f Z%8.5f\n",xyz1[0],xyz1[1],xyz1[2]);
	  printf("xyz2 X%8.5f Y%8.5f Z%8.5f\n",xyz2[0],xyz2[1],xyz2[2]);
	  printf("xyz3 X%8.5f Y%8.5f Z%8.5f\n",xyz3[0],xyz3[1],xyz3[2]);
	  printf("norm0 X%8.5f Y%8.5f Z%8.5f\n",norm0[0],norm0[1],norm0[2]);
	}

	if ( gridDotProduct( leg1, norm0 ) > 0.00 &&
	     gridDotProduct( leg0, norm1 ) > 0.00 &&
	     gridDotProduct( leg0, norm2 ) > 0.00 &&
	     gridDotProduct( leg0, norm3 ) > 0.00 ){
	  foundCell = cell;
	} 
      }
    }
    cell++;
  }

  if ( edgeSplit ) return newnode;

  if (EMPTY == foundCell) return EMPTY;

  return gridSplitCellAt(grid,foundCell,newX,newY,newZ);
}

Grid *gridCollapseEdge(Grid *grid, int n0, int n1, double ratio )
{
  int i, cell, face, face0, face1, faceId;
  double xyz0[3], xyz1[3], xyzAvg[3];
  double uv0[2], uv1[2], uvAvg[2];
  AdjIterator it;
  bool volumeEdge;

  if ( gridGeometryNode(grid, n1) ) return NULL;
  if ( gridGeometryEdge(grid, n1) && EMPTY == gridFindEdge(grid, n0, n1) ) 
    return NULL;

  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;

  volumeEdge = gridContinuousEquator(grid);

  if ( volumeEdge && gridGeometryFace(grid, n0) && gridGeometryFace(grid, n1))
    return NULL;

  if ( NULL == gridNodeXYZ( grid, n0, xyz0) ) return NULL;
  if ( NULL == gridNodeXYZ( grid, n1, xyz1) ) return NULL;
  
  for (i=0 ; i<3 ; i++) {
    xyzAvg[i] = (1.0-ratio) * xyz0[i] + ratio * xyz1[i];
    if ( volumeEdge && gridGeometryFace(grid, n0) ) xyzAvg[i] = xyz0[i];
    if ( volumeEdge && gridGeometryFace(grid, n1) ) xyzAvg[i] = xyz1[i];
    if ( gridGeometryEdge(grid, n0) && !gridGeometryEdge(grid, n1) ) 
      xyzAvg[i] = xyz0[i];
    if ( gridGeometryNode(grid, n0) ) xyzAvg[i] = xyz0[i];
  }

  if ( NULL == gridSetNodeXYZ( grid, n0, xyzAvg) ) return NULL;
  if ( NULL == gridSetNodeXYZ( grid, n1, xyzAvg) ) return NULL;  

  if ( gridNegCellAroundNodeExceptGem( grid, n0 ) || 
       gridNegCellAroundNodeExceptGem( grid, n1 ) ) {
    gridSetNodeXYZ( grid, n0, xyz0);
    gridSetNodeXYZ( grid, n1, xyz1);
    return NULL;
  }

  for (i=0 ; i<gridNGem(grid) ; i++) gridRemoveCell( grid, gridGem(grid,i) );
  
  gridReconnectAllCell(grid, n1, n0 );

  if ( !volumeEdge ) {
    int gap0, gap1, faceId0, faceId1;
    double n0Id0uv[2], n1Id0uv[2], n0Id1uv[2], n1Id1uv[2];
    double newId0uv[2], newId1uv[2]; 
    int edge, edgeId;
    double t0, t1, newT;
    gap0 = gridEqu(grid,0);
    gap1 = gridEqu(grid,gridNGem(grid));
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    faceId0 = gridFaceId(grid, n0, n1, gap0 );
    faceId1 = gridFaceId(grid, n0, n1, gap1 );
    if ( faceId0 == EMPTY || faceId1 == EMPTY ) return NULL;
    gridNodeUV(grid,n0,faceId0,n0Id0uv);
    gridNodeUV(grid,n1,faceId0,n1Id0uv);
    gridNodeUV(grid,n0,faceId1,n0Id1uv);
    gridNodeUV(grid,n1,faceId1,n1Id1uv);
    if ( gridGeometryEdge(grid, n0) && !gridGeometryEdge(grid, n1) ||
	 gridGeometryNode(grid, n0) ) {
      newId0uv[0] = n0Id0uv[0];
      newId0uv[1] = n0Id0uv[1];
      newId1uv[0] = n0Id1uv[0];
      newId1uv[1] = n0Id1uv[1];
    }else{
      newId0uv[0] = (1.0-ratio) * n0Id0uv[0] + ratio * n1Id0uv[0];
      newId0uv[1] = (1.0-ratio) * n0Id0uv[1] + ratio * n1Id0uv[1];
      newId1uv[0] = (1.0-ratio) * n0Id1uv[0] + ratio * n1Id1uv[0];
      newId1uv[1] = (1.0-ratio) * n0Id1uv[1] + ratio * n1Id1uv[1];
    }

    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );

    gridReconnectAllFace(grid, n1, n0);
    gridSetNodeUV(grid, n0, faceId0, newId0uv[0], newId0uv[1]);
    gridSetNodeUV(grid, n0, faceId1, newId1uv[0], newId1uv[1]);

    edge = gridFindEdge(grid,n0,n1);
    if ( edge != EMPTY ) {
      edgeId = gridEdgeId(grid,n0,n1);
      gridNodeT(grid, n0, edgeId, &t0 );
      gridNodeT(grid, n1, edgeId, &t1 );
      newT =  (1.0-ratio) * t0 + ratio * t1;
      if (gridGeometryNode(grid, n0) ) newT = t0;
      gridRemoveEdge(grid,edge);

      gridReconnectAllEdge(grid, n1, n0 );
      gridSetNodeT(grid, n0, edgeId, newT);
    }
  }

  if ( volumeEdge && gridGeometryFace(grid, n1) ) 
    gridReconnectAllFace(grid,n1,n0);

  gridRemoveNode(grid, n1);

  return grid;
}

Grid *gridFreezeGoodNodes(Grid *grid, double goodAR, 
			  double minLength, double maxLength )
{
  int n0, n1;
  double ratio;
  double ar;

  for ( n0=0; n0<gridMaxNode(grid); n0++ ) { 
    if ( gridValidNode( grid, n0) && !gridNodeFrozen( grid, n0 ) ) {
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) 
	return NULL;
      if ( ratio < maxLength ) {
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( ratio > minLength ) { 
	  gridNodeAR(grid,n0,&ar);
	  if ( grid == gridSafeProjectNode(grid, n0, 1.0 ) &&
	       ar > goodAR ) gridFreezeNode(grid,n0);
	}
      }
    }
  }
  return grid;
}

Grid *gridVerifyEdgeExists(Grid *grid, int n0, int n1 )
{
  int i0, i1, nodes0[4], nodes1[4];
  bool gotIt;
  int removeNode;
  double ratio;
  AdjIterator it0, it1;

  gotIt=gridCellEdge( grid, n0, n1 );

  if( !gotIt && n0 != EMPTY && n1 != EMPTY ) {
    for ( it0 = adjFirst(gridCellAdj(grid),n0); 
	  adjValid(it0) && !gotIt; 
	  it0 = adjNext(it0) ) {
      gridCell( grid, adjItem(it0), nodes0 );
      for(i0=0;i0<4&&!gotIt;i0++){
	for ( it1 = adjFirst(gridCellAdj(grid),nodes0[i0]); 
	      adjValid(it1) && !gotIt; 
	      it1 = adjNext(it1) ) {
	  gridCell( grid, adjItem(it1), nodes1 );
	  for(i1=0;i1<4&&!gotIt;i1++){
	    if (nodes1[i1] == n1 && !gridNodeFrozen( grid, nodes0[i0] )) {
	      gridCollapseEdge(grid, n0, nodes0[i0], 0.0);
	      gotIt = gridCellEdge( grid, n0, n1 );
	      if (!gotIt) gridCollapseEdge(grid, n1, nodes0[i0], 0.0);
	      gotIt = gridCellEdge( grid, n0, n1 );
	    }	    
	  }     
	}
      }
    }
  }

  if (gotIt) return grid;
  return NULL;
}

Grid *gridVerifyFaceExists(Grid *grid, int n0, int n1, int n2 )
{
  bool gotIt;

  gotIt=gridCellFace( grid, n0, n1, n2 );

  if (gotIt) return grid;
  return NULL;
}

