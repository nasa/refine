
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
#include "gridmetric.h"
#include "gridswap.h"

#define COSTLIMIT (-0.5)

Grid *gridRemoveTwoFaceCell(Grid *grid, Queue *queue, int cell )
{
  int cellnodes[4], facenodes[3];
  int face, faces[4], faceIds[4];
  int facecount;
  int face0, face1;
  int faceId0, faceId1;
  int node;
  int degree[4];
  double uv[8];
  int newface0, newface1;

  int cell2face[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

  if (grid!=gridCell(grid,cell,cellnodes)) return NULL;

  if (NULL==queue && gridCellHasGhostNode(grid,cellnodes)) return NULL;
  if (NULL!=queue && !gridCellHasGhostNode(grid,cellnodes)) return NULL;

  facecount = 0;
  face0 = face1 = EMPTY;
  for(face=0;face<4;face++) {
    faces[face] = gridFindFace( grid,
				cellnodes[cell2face[face][0]],
				cellnodes[cell2face[face][1]],
				cellnodes[cell2face[face][2]]);
    if (EMPTY==faces[face]) {
      faceIds[face]=EMPTY;
    }else{
      gridFace(grid,faces[face],facenodes,&faceIds[face]);      
    }
    if (faces[face]!=EMPTY) {
      facecount++;
      if (face0!=EMPTY) face1 = face;
      if (face0==EMPTY) face0 = face;
    }
  }

  if (2==facecount) {
    faceId0 = faceIds[face0];
    faceId1 = faceIds[face1];
    if (faceId0==faceId1) {

      /* make sure I'm not missing a face off proc */
      for(node=0;node<4;node++) degree[node]=0;
      for(face=0;face<4;face++)
	if (EMPTY!=faces[face])
	  for(node=0;node<3;node++) degree[cell2face[face][node]]++;
      for(node=0;node<4;node++) 
	if (2==degree[node] && gridNodeGhost(grid,cellnodes[node]))
	  return NULL;

      gridRemoveCellAndQueue(grid,queue,cell);
      for(node=0;node<4;node++) 
	gridNodeUV(grid, cellnodes[node], faceId0, &uv[2*node]);
      
      gridRemoveFaceAndQueue(grid, queue, faces[face0] );
      gridRemoveFaceAndQueue(grid, queue, faces[face1] );

      newface0 = newface1 = EMPTY;
      for(face=0;face<4;face++) {
	if (faces[face]==EMPTY) {
	  if (newface0!=EMPTY) newface1 = face;
	  if (newface0==EMPTY) newface0 = face;
	}
      }

      /* add opposite face in left-handed for the removed tet, 
	 right-handed for rest of grid */
 
      facenodes[0] = cell2face[newface0][1];
      facenodes[1] = cell2face[newface0][0];
      facenodes[2] = cell2face[newface0][2];

      gridAddFaceUVAndQueue(grid, queue, 
			    cellnodes[facenodes[0]], 
			    uv[0+2*facenodes[0]],
			    uv[1+2*facenodes[0]],
			    cellnodes[facenodes[1]], 
			    uv[0+2*facenodes[1]],
			    uv[1+2*facenodes[1]],
			    cellnodes[facenodes[2]], 
			    uv[0+2*facenodes[2]],
			    uv[1+2*facenodes[2]],
			    faceId0 );

      facenodes[0] = cell2face[newface1][1];
      facenodes[1] = cell2face[newface1][0];
      facenodes[2] = cell2face[newface1][2];

      gridAddFaceUVAndQueue(grid, queue, 
			    cellnodes[facenodes[0]], 
			    uv[0+2*facenodes[0]],
			    uv[1+2*facenodes[0]],
			    cellnodes[facenodes[1]], 
			    uv[0+2*facenodes[1]],
			    uv[1+2*facenodes[1]],
			    cellnodes[facenodes[2]], 
			    uv[0+2*facenodes[2]],
			    uv[1+2*facenodes[2]],
			    faceId1 );

      return grid;
    }
  }
  return NULL;
}

Grid *gridSwapFace(Grid *grid, Queue *queue, int n0, int n1, int n2 )
{
  int cell0, cell1;
  int nodes0[4], nodes1[4], topnode, bottomnode;
  int tent[5], nodes[3][4];
  int i;
  double origcost, bestcost;

  if ( gridNodeFrozen( grid, n0 ) && 
       gridNodeFrozen( grid, n1 ) && 
       gridNodeFrozen( grid, n2 ) )return NULL;  

  cell0 = gridFindOtherCellWith3Nodes(grid, n0, n1, n2, EMPTY );
  if ( EMPTY == cell0 ) return NULL;
  cell1 = gridFindOtherCellWith3Nodes(grid, n0, n1, n2, cell0 );
  if ( EMPTY == cell1 ) return NULL;

  gridCell(grid, cell0, nodes0);
  topnode = nodes0[0] + nodes0[1] + nodes0[2] + nodes0[3] - n0 - n1 - n2; 

  tent[0] = topnode;
  if ( topnode == nodes0[1] ) {
    tent[1] = nodes0[0];
  }else{
    tent[1] = nodes0[1];
  }
  gridOrient(grid,nodes0,tent);

  gridCell(grid, cell1, nodes1);
  bottomnode = nodes1[0] + nodes1[1] + nodes1[2] + nodes1[3] - n0 - n1 - n2; 

  //set nodes[0-3] to topnode orientation
  tent[4] = bottomnode;

  nodes[0][0] = tent[0];
  nodes[0][1] = tent[1];
  nodes[0][2] = tent[2];
  nodes[0][3] = tent[4];
  nodes[1][0] = tent[0];
  nodes[1][1] = tent[2];
  nodes[1][2] = tent[3];
  nodes[1][3] = tent[4];
  nodes[2][0] = tent[0];
  nodes[2][1] = tent[3];
  nodes[2][2] = tent[1];
  nodes[2][3] = tent[4];

  origcost = MIN( gridAR( grid, nodes0 ), gridAR( grid, nodes1 ) ); 
  bestcost = MIN( gridAR( grid, nodes[0] ), gridAR( grid, nodes[1] ) ); 
  bestcost = MIN( gridAR( grid, nodes[2] ), bestcost ); 

  if ( bestcost > origcost && bestcost > COSTLIMIT ) {
    gridRemoveCell(grid,cell0);
    gridRemoveCell(grid,cell1);
    for ( i = 0 ; i < 3 ; i++ )
      gridAddCell( grid, nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );
    return grid;
  }
  return NULL;
}

Grid *gridSwapEdge3(Grid *grid, Queue *queue, int n0, int n1 );
Grid *gridSwapEdge4(Grid *grid, Queue *queue, int n0, int n1 );
Grid *gridSwapEdge5(Grid *grid, Queue *queue, int n0, int n1 );
Grid *gridSwapEdge6(Grid *grid, Queue *queue, int n0, int n1 );
Grid *gridSwapEdge7(Grid *grid, Queue *queue, int n0, int n1 );

Grid *gridSwapEdge(Grid *grid, Queue *queue, int n0, int n1 )
{
  int gap0, gap1, face0, face1, faceId0, faceId1, newFaceId0, newFaceId1;
  double origMR, newMR;
  Grid *swapStatus;

  if ( gridNodeFrozen( grid, n0 ) && gridNodeFrozen( grid, n1 ) )return NULL;  
  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;
  
  //test face
  gap0 = gap1 = EMPTY;
  face0 = face1 = EMPTY;
  faceId0 = EMPTY;
  if ( !gridContinuousEquator(grid) ){
    gap0 = gridEqu(grid,0);
    gap1 = gridEqu(grid,gridNGem(grid));
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    faceId0 = gridFaceId(grid, n0, n1, gap0 );
    faceId1 = gridFaceId(grid, n0, n1, gap1 );
    
    if ( faceId0 == EMPTY || faceId1 == EMPTY ) return NULL;
    if ( faceId0 != faceId1 ) return NULL;

    newFaceId0 = gridFaceId(grid, n0, gap0, gap1 );
    newFaceId1 = gridFaceId(grid, n1, gap0, gap1 );
    if ( newFaceId0 != EMPTY || newFaceId1 != EMPTY ) return NULL;

    if (gridCOST_FCN_EDGE_LENGTH != gridCostFunction(grid)) {
      // make sure that face MR can improve
      origMR = MIN( gridFaceMR(grid, n0, n1, gap0 ), 
		    gridFaceMR(grid, n0, n1, gap1 ) );
      newMR  = MIN( gridFaceMR(grid, n0, gap0, gap1 ),
		    gridFaceMR(grid, n1, gap0, gap1 ) );
      if ( origMR > newMR ) return NULL;
    }
  }

  switch (gridNEqu(grid)) {
  case 3: 
    swapStatus = gridSwapEdge3(grid, queue, n0, n1); break;
  case 4: 
    swapStatus = gridSwapEdge4(grid, queue, n0, n1); break;
  case 5:
    swapStatus = gridSwapEdge5(grid, queue, n0, n1); break;
  case 6:
    swapStatus = gridSwapEdge6(grid, queue, n0, n1); break;
  case 7:
    swapStatus = gridSwapEdge7(grid, queue, n0, n1); break;
  default:
    swapStatus = NULL; break;
  }

  if ( !gridContinuousEquator(grid) && swapStatus != NULL ) {
    double n0uv[2], n1uv[2], gap0uv[2], gap1uv[2]; 

    gridNodeUV(grid, n0,   faceId0, n0uv);
    gridNodeUV(grid, n1,   faceId0, n1uv);
    gridNodeUV(grid, gap0, faceId0, gap0uv);
    gridNodeUV(grid, gap1, faceId0, gap1uv);

    gridRemoveFaceAndQueue(grid, queue, face0 );
    gridRemoveFaceAndQueue(grid, queue, face1 );

    gridAddFaceUVAndQueue(grid, queue, 
			  n0,   n0uv[0],   n0uv[1], 
			  gap1, gap1uv[0], gap1uv[1], 
			  gap0, gap0uv[0], gap0uv[1], 
			  faceId0 );
    gridAddFaceUVAndQueue(grid, queue, 
			  n1,   n1uv[0],   n1uv[1], 
			  gap0, gap0uv[0], gap0uv[1], 
			  gap1, gap1uv[0], gap1uv[1], 
			  faceId0 );
  }

  return swapStatus;
}

Grid *gridSwapNearNode(Grid *grid, int node)
{
  int nswap, nodes[4];
  AdjIterator it;

  nswap = 0;
  it = adjFirst(gridCellAdj(grid),node);
  while ( adjValid(it) ){
    gridCell( grid, adjItem(it), nodes);
    if ( gridAR(grid, nodes) < 0.5 ) {
      if ( ( NULL != gridSwapFace( grid, NULL, nodes[1],nodes[2],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[2],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[1],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[1],nodes[2]) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[1] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[2], nodes[3] ) ) ) {
	it = adjFirst(gridCellAdj(grid),node);
	nswap++;
      }else{
	it = adjNext(it);
      } 
    }else{
      it = adjNext(it);
    }
    if (nswap>100) {
      printf("node %d swap out.",node);
      return NULL;
    }
  }

  return grid;
}

Grid *gridSwapNearNodeExceptBoundary(Grid *grid, int node)
{
  int nswap, nodes[4];
  AdjIterator it;

  nswap = 0;
  it = adjFirst(gridCellAdj(grid),node);
  while ( adjValid(it) ){
    gridCell( grid, adjItem(it), nodes);
    if ( gridAR(grid, nodes) < 0.5 ) {
      if ( ( !gridGeometryFace(grid, nodes[0]) && 
	     !gridGeometryFace(grid, nodes[1]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[1] ) ) ) ||
	   ( !gridGeometryFace(grid, nodes[0]) && 
	     !gridGeometryFace(grid, nodes[2]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[2] ) ) ) ||
	   ( !gridGeometryFace(grid, nodes[0]) && 
	     !gridGeometryFace(grid, nodes[3]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[3] ) ) ) ||
	   ( !gridGeometryFace(grid, nodes[1]) && 
	     !gridGeometryFace(grid, nodes[2]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[2] ) ) ) ||
	   ( !gridGeometryFace(grid, nodes[1]) && 
	     !gridGeometryFace(grid, nodes[3]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[3] ) ) ) ||
	   ( !gridGeometryFace(grid, nodes[2]) && 
	     !gridGeometryFace(grid, nodes[3]) &&
	     ( NULL != gridSwapEdge( grid, NULL, nodes[2], nodes[3] ) ) ) ) {
	it = adjFirst(gridCellAdj(grid),node);
	nswap++;
      }else{
	it = adjNext(it);
      } 
    }else{
      it = adjNext(it);
    }
    if (nswap>100) {
      printf("node %d swap out.",node);
      return NULL;
    }
  }

  return grid;
}

Grid *gridSwap(Grid *grid, double improvementLimit)
{
  int cellId, maxcell;
  int nodes[4];
  GridBool swap;

  if ( improvementLimit < 0.0 ) improvementLimit = 0.5;

  maxcell = gridMaxCell(grid);

  for (cellId=0;cellId<maxcell;cellId++){
    gridRemoveTwoFaceCell(grid, NULL, cellId );
    if ( grid == gridCell( grid, cellId, nodes) && 
	 gridAR(grid, nodes) < improvementLimit) {
      swap = TRUE;
      if (swap) swap = (grid != gridSwapFace(grid, NULL,nodes[1],nodes[2],nodes[3]) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = (grid != gridSwapFace(grid, NULL,nodes[0],nodes[2],nodes[3]) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = (grid != gridSwapFace(grid, NULL,nodes[0],nodes[1],nodes[3]) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = (grid != gridSwapFace(grid, NULL,nodes[0],nodes[1],nodes[2]) )
		  || ( grid == gridCell( grid, cellId, nodes) );

      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[0], nodes[1] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[0], nodes[2] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[0], nodes[3] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[1], nodes[2] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[1], nodes[3] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
      if (swap) swap = ( grid != gridSwapEdge( grid, NULL, nodes[2], nodes[3] ) )
		  || ( grid == gridCell( grid, cellId, nodes) );
    }
  }

  return grid;
}

Grid *gridSwapEdge3(Grid *grid, Queue *queue, int n0, int n1 )
{
  int i, nodes[2][4];
  double origcost, bestcost;

  gridGemAR(grid, &origcost);

  nodes[0][0]=n0;
  nodes[0][1]=gridEqu(grid,0);
  nodes[0][2]=gridEqu(grid,1);
  nodes[0][3]=gridEqu(grid,2);
  nodes[1][0]=n1;
  nodes[1][1]=gridEqu(grid,0);
  nodes[1][2]=gridEqu(grid,2);
  nodes[1][3]=gridEqu(grid,1);

  bestcost = MIN( gridAR( grid, nodes[0] ), gridAR( grid, nodes[1] ) );

  if ( bestcost > origcost && bestcost > COSTLIMIT ) {

    gridRemoveGemAndQueue(grid,queue);

    for ( i = 0 ; i < 2 ; i++ )
      gridAddCellAndQueue( grid, queue, 
			   nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );

    return grid;
  }
  return NULL;
}

Grid *gridSwapEdge4(Grid *grid, Queue *queue, int n0, int n1 )
{
  int i, nodes[4][4], bestindex;
  double cost, origcost, currentcost, bestcost;

  gridGemAR(grid, &origcost);

  nodes[0][0]=n0;
  nodes[0][1]=gridEqu(grid,0);
  nodes[0][2]=gridEqu(grid,1);
  nodes[0][3]=gridEqu(grid,2);
  nodes[1][0]=n0;
  nodes[1][1]=gridEqu(grid,2);
  nodes[1][2]=gridEqu(grid,3);
  nodes[1][3]=gridEqu(grid,4);
  nodes[2][0]=n1;
  nodes[2][1]=gridEqu(grid,0);
  nodes[2][2]=gridEqu(grid,2);
  nodes[2][3]=gridEqu(grid,1);
  nodes[3][0]=n1;
  nodes[3][1]=gridEqu(grid,2);
  nodes[3][2]=gridEqu(grid,0);
  nodes[3][3]=gridEqu(grid,3);

  currentcost = 2.0;

  for ( i = 0 ; i < 4 ; i++ ) {
    cost = gridAR( grid, nodes[i] );
    currentcost = MIN(currentcost,cost);
  }

  bestcost = currentcost;
  bestindex = 0;

  nodes[0][0]=n0;
  nodes[0][1]=gridEqu(grid,1);
  nodes[0][2]=gridEqu(grid,3);
  nodes[0][3]=gridEqu(grid,0);
  nodes[1][0]=n0;
  nodes[1][1]=gridEqu(grid,3);
  nodes[1][2]=gridEqu(grid,1);
  nodes[1][3]=gridEqu(grid,2);
  nodes[2][0]=n1;
  nodes[2][1]=gridEqu(grid,3);
  nodes[2][2]=gridEqu(grid,1);
  nodes[2][3]=gridEqu(grid,0);
  nodes[3][0]=n1;
  nodes[3][1]=gridEqu(grid,1);
  nodes[3][2]=gridEqu(grid,3);
  nodes[3][3]=gridEqu(grid,2);

  currentcost = 2.0;

  for ( i = 0 ; i < 4 ; i++ ) {
    cost = gridAR( grid, nodes[i] );
    currentcost = MIN(currentcost,cost);
  }

  if ( currentcost > bestcost ) {
    bestcost = currentcost;
    bestindex = 1;
  }

  if ( bestcost > origcost && bestcost > COSTLIMIT ) {

    if (bestindex == 0){
      nodes[0][0]=n0;
      nodes[0][1]=gridEqu(grid,0);
      nodes[0][2]=gridEqu(grid,1);
      nodes[0][3]=gridEqu(grid,2);
      nodes[1][0]=n0;
      nodes[1][1]=gridEqu(grid,2);
      nodes[1][2]=gridEqu(grid,3);
      nodes[1][3]=gridEqu(grid,4);
      nodes[2][0]=n1;
      nodes[2][1]=gridEqu(grid,0);
      nodes[2][2]=gridEqu(grid,2);
      nodes[2][3]=gridEqu(grid,1);
      nodes[3][0]=n1;
      nodes[3][1]=gridEqu(grid,2);
      nodes[3][2]=gridEqu(grid,0);
      nodes[3][3]=gridEqu(grid,3);
    }else{
      nodes[0][0]=n0;
      nodes[0][1]=gridEqu(grid,1);
      nodes[0][2]=gridEqu(grid,3);
      nodes[0][3]=gridEqu(grid,0);
      nodes[1][0]=n0;
      nodes[1][1]=gridEqu(grid,3);
      nodes[1][2]=gridEqu(grid,1);
      nodes[1][3]=gridEqu(grid,2);
      nodes[2][0]=n1;
      nodes[2][1]=gridEqu(grid,3);
      nodes[2][2]=gridEqu(grid,1);
      nodes[2][3]=gridEqu(grid,0);
      nodes[3][0]=n1;
      nodes[3][1]=gridEqu(grid,1);
      nodes[3][2]=gridEqu(grid,3);
      nodes[3][3]=gridEqu(grid,2);
    }

    gridRemoveGemAndQueue(grid,queue);
    
    for ( i = 0 ; i < 4 ; i++ )
      gridAddCellAndQueue( grid, queue,
			   nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );

    return grid;
  }
  return NULL;
}

Grid *gridSwapEdge5(Grid *grid, Queue *queue, int n0, int n1 )
{
  int i;
  int currentindex, bestindex, nodes[6][4];

  double cost, origcost, currentcost, bestcost;

  gridGemAR(grid, &origcost);

  bestcost  =  -2.0;
  bestindex = -1;

  for ( currentindex = 0 ; currentindex < 5 ; currentindex++ ) {
    nodes[0][0]=n0;
    nodes[0][1]=gridEqu(grid,0);
    nodes[0][2]=gridEqu(grid,1);
    nodes[0][3]=gridEqu(grid,2);
    nodes[1][0]=n0;
    nodes[1][1]=gridEqu(grid,0);
    nodes[1][2]=gridEqu(grid,2);
    nodes[1][3]=gridEqu(grid,3);
    nodes[2][0]=n0;
    nodes[2][1]=gridEqu(grid,0);
    nodes[2][2]=gridEqu(grid,3);
    nodes[2][3]=gridEqu(grid,4);
    nodes[3][0]=n1;
    nodes[3][1]=gridEqu(grid,0);
    nodes[3][2]=gridEqu(grid,2);
    nodes[3][3]=gridEqu(grid,1);
    nodes[4][0]=n1;
    nodes[4][1]=gridEqu(grid,0);
    nodes[4][2]=gridEqu(grid,3);
    nodes[4][3]=gridEqu(grid,2);
    nodes[5][0]=n1;
    nodes[5][1]=gridEqu(grid,0);
    nodes[5][2]=gridEqu(grid,4);
    nodes[5][3]=gridEqu(grid,3);

    currentcost = 2.0;

    for ( i = 0 ; i < 6 ; i++ ) {
      cost = gridAR( grid, nodes[i] );
      currentcost = MIN(currentcost,cost);
    } 

    if ( currentcost > bestcost ) {
      bestcost = currentcost;
      bestindex = currentindex;
    }

    gridCycleEquator( grid );
  }

  if ( bestcost > origcost && bestcost > COSTLIMIT ) {

    for ( i = 0 ; i < bestindex ; i++ ) 
      gridCycleEquator( grid );
    
    nodes[0][0]=n0;
    nodes[0][1]=gridEqu(grid,0);
    nodes[0][2]=gridEqu(grid,1);
    nodes[0][3]=gridEqu(grid,2);
    nodes[1][0]=n0;
    nodes[1][1]=gridEqu(grid,0);
    nodes[1][2]=gridEqu(grid,2);
    nodes[1][3]=gridEqu(grid,3);
    nodes[2][0]=n0;
    nodes[2][1]=gridEqu(grid,0);
    nodes[2][2]=gridEqu(grid,3);
    nodes[2][3]=gridEqu(grid,4);
    nodes[3][0]=n1;
    nodes[3][1]=gridEqu(grid,0);
    nodes[3][2]=gridEqu(grid,2);
    nodes[3][3]=gridEqu(grid,1);
    nodes[4][0]=n1;
    nodes[4][1]=gridEqu(grid,0);
    nodes[4][2]=gridEqu(grid,3);
    nodes[4][3]=gridEqu(grid,2);
    nodes[5][0]=n1;
    nodes[5][1]=gridEqu(grid,0);
    nodes[5][2]=gridEqu(grid,4);
    nodes[5][3]=gridEqu(grid,3);

    gridRemoveGemAndQueue(grid,queue);
    
    for ( i = 0 ; i < 6 ; i++ )
      gridAddCellAndQueue( grid, queue,
			   nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );
    return grid;
  }  
  return NULL;
}

Grid *gridGetCombo6( Grid *grid, int nodes[40][4], double costs[20], 
		     double *bestcost, int best[4] );
Grid *gridConstructTet6( Grid *grid, int n0, int n1, 
			 int nodes[40][4], double costs[20] );

Grid *gridSwapEdge6( Grid *grid, Queue *queue, int n0, int n1 )
{
  int i;
  int nodes[40][4], bestcombo[4];

  double costs[20], origcost, bestcost;
  
  gridGemAR(grid, &origcost);

  gridConstructTet6( grid, n0, n1, nodes, costs );
  
  gridGetCombo6( grid, nodes, costs, &bestcost, bestcombo );
  
  if ( bestcost > origcost && bestcost > COSTLIMIT ) {
      
    gridRemoveGemAndQueue(grid,queue);
    
    for ( i = 0 ; i < 4 ; i++ ){
      gridAddCellAndQueue( grid, queue,
			   nodes[bestcombo[i]][0], 
			   nodes[bestcombo[i]][1], 
			   nodes[bestcombo[i]][2], 
			   nodes[bestcombo[i]][3] );
      gridAddCellAndQueue( grid,  queue,
			   nodes[bestcombo[i]+20][0], 
			   nodes[bestcombo[i]+20][1], 
			   nodes[bestcombo[i]+20][2], 
			   nodes[bestcombo[i]+20][3] );
    }
    return grid;
  }
  return NULL;
}

Grid *gridGetCombo6( Grid *grid, int nodes[40][4], double costs[20], 
		     double *bestcost, int best[4] )
{  
  int i, j, tet[4];
  double cost;
  *bestcost = -2.0;

  for ( i = 0 ; i < 6 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0]+6;
    tet[2] = i+4; if ( tet[2] > 5 ) tet[2] -= 6;
    tet[3] = tet[2] + 12;
    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], costs[tet[3]] ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 4 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap6 %d currentcost %f\n",i,cost);
#endif
  }
  
  for ( i = 0 ; i < 3 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 12;
    tet[2] = i+3; if ( tet[2] > 5 ) tet[2] -= 6;
    tet[3] = tet[2] + 12;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], costs[tet[3]] ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 4 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap6 %d currentcost %f\n",i,cost);
#endif
  }
  
  for ( i = 0 ; i < 3 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 6;
    tet[2] = i+3; if ( tet[2] > 5 ) tet[2] -= 6;
    tet[3] = tet[2] + 6;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], costs[tet[3]] ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 4 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap6 %d currentcost %f\n",i,cost);
#endif
  }
  
  
  for ( i = 0 ; i < 2 ; i++ ) {
    tet[0] = i;
    tet[1] = i+2;
    tet[2] = i+4;
    tet[3] = tet[0] + 18;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], costs[tet[3]] ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 4 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap6 %d currentcost %f\n",i,cost);
#endif
  }  
  return grid;
}

Grid *gridConstructTet6( Grid *grid, int n0, int n1, 
			 int nodes[40][4], double costs[20] )
{
  int i;

  /* make the small triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i   ][0]=gridEqu(grid,i);
    nodes[i   ][1]=gridEqu(grid,i+5);
    nodes[i   ][2]=gridEqu(grid,i+1);
    nodes[i   ][3]=n0;
    nodes[i+20][0]=gridEqu(grid,i+1);
    nodes[i+20][1]=gridEqu(grid,i+5);
    nodes[i+20][2]=gridEqu(grid,i);
    nodes[i+20][3]=n1;
  }

  /* make the next triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i+06][0]=gridEqu(grid,i+5);
    nodes[i+06][1]=gridEqu(grid,i+2);
    nodes[i+06][2]=gridEqu(grid,i+1);
    nodes[i+06][3]=n0;
    nodes[i+26][0]=gridEqu(grid,i+1);
    nodes[i+26][1]=gridEqu(grid,i+2);
    nodes[i+26][2]=gridEqu(grid,i+5);
    nodes[i+26][3]=n1;
  }

  /* make the previous triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i+12][0]=gridEqu(grid,i+5);
    nodes[i+12][1]=gridEqu(grid,i+4);
    nodes[i+12][2]=gridEqu(grid,i+1);
    nodes[i+12][3]=n0;
    nodes[i+32][0]=gridEqu(grid,i+1);
    nodes[i+32][1]=gridEqu(grid,i+4);
    nodes[i+32][2]=gridEqu(grid,i+5);
    nodes[i+32][3]=n1;
  }

  /* make the big triangles */
  for ( i = 0; i < 2; i++ ){
    nodes[i+18][0]=gridEqu(grid,i+5);
    nodes[i+18][1]=gridEqu(grid,i+3);
    nodes[i+18][2]=gridEqu(grid,i+1);
    nodes[i+18][3]=n0;
    nodes[i+38][0]=gridEqu(grid,i+1);
    nodes[i+38][1]=gridEqu(grid,i+3);
    nodes[i+38][2]=gridEqu(grid,i+5);
    nodes[i+38][3]=n1;
  }

  for ( i = 0; i < 20; i++ )
    costs[i] = MIN( gridAR(grid, nodes[i]), gridAR(grid, nodes[i+20]) );
  
  return grid;
}

Grid *gridGetCombo7( Grid *grid, int nodes[70][4], double costs[35], 
		     double *bestcost, int best[4] );
Grid *gridConstructTet7( Grid *grid, int n0, int n1, 
			 int nodes[70][4], double costs[35] );

Grid *gridSwapEdge7( Grid *grid, Queue *queue, int n0, int n1 )
{
  int i;
  int nodes[70][4], bestcombo[5];

  double costs[35], origcost, bestcost;
  
  gridGemAR(grid, &origcost);

  gridConstructTet7( grid, n0, n1, nodes, costs );

  gridGetCombo7( grid, nodes, costs, &bestcost, bestcombo );
  
  if ( bestcost > origcost && bestcost > COSTLIMIT ) {

    gridRemoveGemAndQueue(grid,queue);
    
    for ( i = 0 ; i < 5 ; i++ ){
      gridAddCellAndQueue( grid, queue,
			   nodes[bestcombo[i]][0], 
			   nodes[bestcombo[i]][1], 
			   nodes[bestcombo[i]][2], 
			   nodes[bestcombo[i]][3] );
      gridAddCellAndQueue( grid,  queue,
			   nodes[bestcombo[i]+35][0], 
			   nodes[bestcombo[i]+35][1], 
			   nodes[bestcombo[i]+35][2], 
			   nodes[bestcombo[i]+35][3] );
    }
    return grid;
  }  
  return NULL;
}

Grid *gridGetCombo7( Grid *grid, int nodes[70][4], double costs[35], 
		     double *bestcost, int best[5] )
{  
  int i, j, tet[5];
  double cost;
  *bestcost = -2.0;
  //case1
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0]+7;
    tet[2] = i+6; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 21;
    tet[3] = i+5; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = tet[3] + 14;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
#ifdef DEBUGSWAP
      printf("case1 i %d tet %d %d %d %d %d \n",i,tet[0],tet[1],tet[2],tet[3],tet[4]);
#endif
    }
  }
  //case2
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0]+7;
    tet[2] = i+2; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 21;
    tet[3] = i+4; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = tet[3] + 14;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap7 %d currentcost %f\n",i,cost);
#endif
  }
  //case3
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 14;
    tet[2] = i+5; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 21;
    tet[3] = i+4; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = tet[3] + 14;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap7 %d currentcost %f\n",i,cost);
#endif
  }
  //case4  
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 7;
    tet[2] = i+6; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 21;
    tet[3] = i+4; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = tet[3] + 7;

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap7 %d currentcost %f\n",i,cost);
#endif
  }
  //case5
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 7;
    tet[2] = i+4; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 28;
    tet[3] = i+3; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = i+5; if ( tet[4] > 6 ) tet[4] -= 7; 

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap7 %d currentcost %f\n",i,cost);
#endif
  }
  //case6 
  for ( i = 0 ; i < 7 ; i++ ) {
    tet[0] = i;
    tet[1] = tet[0] + 14;
    tet[2] = i+3; if ( tet[2] > 6 ) tet[2] -= 7; tet[2] += 28;
    tet[3] = i+2; if ( tet[3] > 6 ) tet[3] -= 7; 
    tet[4] = i+4; if ( tet[4] > 6 ) tet[4] -= 7; 

    cost = MIN( MIN( costs[tet[0]], costs[tet[1]] ), 
		MIN( costs[tet[2]], MIN( costs[tet[3]], costs[tet[4]] ) ) );
    if ( cost > *bestcost ) {
      *bestcost = cost;
      for ( j = 0 ; j < 5 ; j++ ) best[j] = tet[j]; 
    }
#ifdef DEBUGSWAP
        printf("swap7 %d currentcost %f\n",i,cost);
#endif
  }

  return grid;
}

Grid *gridConstructTet7( Grid *grid, int n0, int n1,
			 int nodes[70][4], double costs[35] )
{
  int i;

  /* make the small triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i   ][0]=gridEqu(grid,i);
    nodes[i   ][1]=gridEqu(grid,i+6);
    nodes[i   ][2]=gridEqu(grid,i+1);
    nodes[i   ][3]=n0;
    nodes[i+35][0]=gridEqu(grid,i+1);
    nodes[i+35][1]=gridEqu(grid,i+6);
    nodes[i+35][2]=gridEqu(grid,i);
    nodes[i+35][3]=n1;
  }

  /* make the next triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+07][0]=gridEqu(grid,i+6);
    nodes[i+07][1]=gridEqu(grid,i+2);
    nodes[i+07][2]=gridEqu(grid,i+1);
    nodes[i+07][3]=n0;
    nodes[i+42][0]=gridEqu(grid,i+1);
    nodes[i+42][1]=gridEqu(grid,i+2);
    nodes[i+42][2]=gridEqu(grid,i+6);
    nodes[i+42][3]=n1;
  }

  /* make the previous triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+14][0]=gridEqu(grid,i+6);
    nodes[i+14][1]=gridEqu(grid,i+5);
    nodes[i+14][2]=gridEqu(grid,i+1);
    nodes[i+14][3]=n0;
    nodes[i+49][0]=gridEqu(grid,i+1);
    nodes[i+49][1]=gridEqu(grid,i+5);
    nodes[i+49][2]=gridEqu(grid,i+6);
    nodes[i+49][3]=n1;
  }

  /* make the spike */
  for ( i = 0; i < 7; i++ ){
    nodes[i+21][0]=gridEqu(grid,i);
    nodes[i+21][1]=gridEqu(grid,i+4);
    nodes[i+21][2]=gridEqu(grid,i+3);
    nodes[i+21][3]=n0;
    nodes[i+56][0]=gridEqu(grid,i+3);
    nodes[i+56][1]=gridEqu(grid,i+4);
    nodes[i+56][2]=gridEqu(grid,i);
    nodes[i+56][3]=n1;
  }

  /* make the big triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+28][0]=gridEqu(grid,i);
    nodes[i+28][1]=gridEqu(grid,i+5);
    nodes[i+28][2]=gridEqu(grid,i+2);
    nodes[i+28][3]=n0;
    nodes[i+63][0]=gridEqu(grid,i+2);
    nodes[i+63][1]=gridEqu(grid,i+5);
    nodes[i+63][2]=gridEqu(grid,i);
    nodes[i+63][3]=n1;
  }

  for ( i = 0; i < 35; i++ )
    costs[i] = MIN( gridAR( grid, nodes[i] ) , gridAR( grid, nodes[i+35]) ); 
  
  return grid;
}
