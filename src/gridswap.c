
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
#include "gridStruct.h"
#include "gridmetric.h"
#include "gridswap.h"

Grid *gridSwapEdge3(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge4(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge5(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge6(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge7(Grid *grid, int n0, int n1 );

Grid *gridSwapEdge(Grid *grid, int n0, int n1 )
{
  int gap0, gap1, face0, face1, faceId0, faceId1, newFaceId0, newFaceId1;
  double origMR, newMR;
  Grid *swapStatus;

  if ( gridNodeFrozen( grid, n0 ) && gridNodeFrozen( grid, n1 ) )return NULL;  
  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;
  
  //test face
  if ( grid->nequ != grid->ngem ){
    gap0 = grid->equ[0];
    gap1 = grid->equ[grid->ngem];
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    faceId0 = gridFaceId(grid, n0, n1, gap0 );
    faceId1 = gridFaceId(grid, n0, n1, gap1 );
    
    if ( faceId0 == EMPTY || faceId1 == EMPTY ) return NULL;
    if ( faceId0 != faceId1 ) return NULL;

    newFaceId0 = gridFaceId(grid, n0, gap0, gap1 );
    newFaceId1 = gridFaceId(grid, n1, gap0, gap1 );
    if ( newFaceId0 != EMPTY || newFaceId1 != EMPTY ) return NULL;

    // make sure that face MR can improve
    origMR = MIN( gridFaceMR(grid, n0, n1, gap0 ), 
		  gridFaceMR(grid, n0, n1, gap1 ) );
    newMR  = MIN( gridFaceMR(grid, n0, gap0, gap1 ),
		  gridFaceMR(grid, n1, gap0, gap1 ) );
    if ( origMR > newMR ) return NULL;
  }

  switch (grid->nequ) {
  case 3: 
    swapStatus = gridSwapEdge3(grid, n0, n1); break;
  case 4: 
    swapStatus = gridSwapEdge4(grid, n0, n1); break;
  case 5:
    swapStatus = gridSwapEdge5(grid, n0, n1); break;
  case 6:
    swapStatus = gridSwapEdge6(grid, n0, n1); break;
  case 7:
    swapStatus = gridSwapEdge7(grid, n0, n1); break;
  default:
    swapStatus = NULL; break;
  }

  if ( grid->nequ != grid->ngem && swapStatus != NULL ) {
    double n0uv[2], n1uv[2], gap0uv[2], gap1uv[2]; 

    gridNodeUV(grid, n0,   faceId0, n0uv);
    gridNodeUV(grid, n1,   faceId0, n1uv);
    gridNodeUV(grid, gap0, faceId0, gap0uv);
    gridNodeUV(grid, gap1, faceId0, gap1uv);

    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );

    gridAddFaceUV(grid, 
		  n0,   n0uv[0],   n0uv[1], 
		  gap1, gap1uv[0], gap1uv[1], 
		  gap0, gap0uv[0], gap0uv[1], 
		  faceId0 );
    gridAddFaceUV(grid, 
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
  it = adjFirst(grid->cellAdj,node);
  while ( adjValid(it) ){
    gridCell( grid, adjItem(it), nodes);
    if ( gridAR(grid, nodes) < 0.5 ) {
      if ( ( NULL != gridSwapEdge( grid, nodes[0], nodes[1] ) ) ||
	   ( NULL != gridSwapEdge( grid, nodes[0], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, nodes[0], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, nodes[1], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, nodes[1], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, nodes[2], nodes[3] ) ) ) {
	it = adjFirst(grid->cellAdj,node);
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
  it = adjFirst(grid->cellAdj,node);
  while ( adjValid(it) ){
    gridCell( grid, adjItem(it), nodes);
    if ( gridAR(grid, nodes) < 0.5 ) {
      if ( !gridGeometryFace(grid, nodes[0]) && 
	   !gridGeometryFace(grid, nodes[1]) &&
	   ( NULL != gridSwapEdge( grid, nodes[0], nodes[1] ) ) ||
	   !gridGeometryFace(grid, nodes[0]) && 
	   !gridGeometryFace(grid, nodes[2]) &&
	   ( NULL != gridSwapEdge( grid, nodes[0], nodes[2] ) ) ||
	   !gridGeometryFace(grid, nodes[0]) && 
	   !gridGeometryFace(grid, nodes[3]) &&
	   ( NULL != gridSwapEdge( grid, nodes[0], nodes[3] ) ) ||
	   !gridGeometryFace(grid, nodes[1]) && 
	   !gridGeometryFace(grid, nodes[2]) &&
	   ( NULL != gridSwapEdge( grid, nodes[1], nodes[2] ) ) ||
	   !gridGeometryFace(grid, nodes[1]) && 
	   !gridGeometryFace(grid, nodes[3]) &&
	   ( NULL != gridSwapEdge( grid, nodes[1], nodes[3] ) ) ||
	   !gridGeometryFace(grid, nodes[2]) && 
	   !gridGeometryFace(grid, nodes[3]) &&
	   ( NULL != gridSwapEdge( grid, nodes[2], nodes[3] ) ) ) {
	it = adjFirst(grid->cellAdj,node);
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

Grid *gridSwap(Grid *grid)
{
  int cellId, nodes[4];
  for (cellId=0;cellId<grid->maxcell;cellId++){
    if ( NULL != gridCell( grid, cellId, nodes) )
      if ( gridAR(grid, nodes) < 0.5 ) {
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[0], nodes[1] );
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[0], nodes[2] );
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[0], nodes[3] );
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[1], nodes[2] );
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[1], nodes[3] );
	if ( NULL != gridCell( grid, cellId, nodes) )
	  gridSwapEdge( grid, nodes[2], nodes[3] );
      }
  }
  return grid;
}

Grid *gridSwapEdge3(Grid *grid, int n0, int n1 )
{
  int i, nodes[2][4];
  double cost, origcost, bestcost;

  origcost = 2.0;

  for ( i = 0 ; i < grid->ngem ; i++ ){
    cost = gridAR( grid, &grid->c2n[4*grid->gem[i]] );
    origcost = MIN(origcost,cost);
  }

  nodes[0][0]=n0;
  nodes[0][1]=grid->equ[0];
  nodes[0][2]=grid->equ[1];
  nodes[0][3]=grid->equ[2];
  nodes[1][0]=n1;
  nodes[1][1]=grid->equ[0];
  nodes[1][2]=grid->equ[2];
  nodes[1][3]=grid->equ[1];

  bestcost = MIN( gridAR( grid, nodes[0] ), gridAR( grid, nodes[1] ) );

  if ( bestcost > origcost ) {

    for ( i = 0 ; i < grid->ngem ; i++ ) 
      gridRemoveCell( grid, grid->gem[i] );
    
    for ( i = 0 ; i < 2 ; i++ )
      gridAddCell( grid, nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );

    return grid;
  }
  return NULL;
}

Grid *gridSwapEdge4(Grid *grid, int n0, int n1 )
{
  int i, nodes[4][4], bestindex;
  double cost, origcost, currentcost, bestcost;

  origcost = 2.0;

  for ( i = 0 ; i < grid->ngem ; i++ ){
    cost = gridAR( grid, &grid->c2n[4*grid->gem[i]] );
    origcost = MIN(origcost,cost);
  }

  nodes[0][0]=n0;
  nodes[0][1]=grid->equ[0];
  nodes[0][2]=grid->equ[1];
  nodes[0][3]=grid->equ[2];
  nodes[1][0]=n0;
  nodes[1][1]=grid->equ[2];
  nodes[1][2]=grid->equ[3];
  nodes[1][3]=grid->equ[4];
  nodes[2][0]=n1;
  nodes[2][1]=grid->equ[0];
  nodes[2][2]=grid->equ[2];
  nodes[2][3]=grid->equ[1];
  nodes[3][0]=n1;
  nodes[3][1]=grid->equ[2];
  nodes[3][2]=grid->equ[0];
  nodes[3][3]=grid->equ[3];

  currentcost = 2.0;

  for ( i = 0 ; i < 4 ; i++ ) {
    cost = gridAR( grid, nodes[i] );
    currentcost = MIN(currentcost,cost);
  }

  bestcost = currentcost;
  bestindex = 0;

  nodes[0][0]=n0;
  nodes[0][1]=grid->equ[1];
  nodes[0][2]=grid->equ[3];
  nodes[0][3]=grid->equ[0];
  nodes[1][0]=n0;
  nodes[1][1]=grid->equ[3];
  nodes[1][2]=grid->equ[1];
  nodes[1][3]=grid->equ[2];
  nodes[2][0]=n1;
  nodes[2][1]=grid->equ[3];
  nodes[2][2]=grid->equ[1];
  nodes[2][3]=grid->equ[0];
  nodes[3][0]=n1;
  nodes[3][1]=grid->equ[1];
  nodes[3][2]=grid->equ[3];
  nodes[3][3]=grid->equ[2];

  currentcost = 2.0;

  for ( i = 0 ; i < 4 ; i++ ) {
    cost = gridAR( grid, nodes[i] );
    currentcost = MIN(currentcost,cost);
  }

  if ( currentcost > bestcost ) {
    bestcost = currentcost;
    bestindex = 1;
  }

  if ( bestcost > origcost ) {

    if (bestindex == 0){
      nodes[0][0]=n0;
      nodes[0][1]=grid->equ[0];
      nodes[0][2]=grid->equ[1];
      nodes[0][3]=grid->equ[2];
      nodes[1][0]=n0;
      nodes[1][1]=grid->equ[2];
      nodes[1][2]=grid->equ[3];
      nodes[1][3]=grid->equ[4];
      nodes[2][0]=n1;
      nodes[2][1]=grid->equ[0];
      nodes[2][2]=grid->equ[2];
      nodes[2][3]=grid->equ[1];
      nodes[3][0]=n1;
      nodes[3][1]=grid->equ[2];
      nodes[3][2]=grid->equ[0];
      nodes[3][3]=grid->equ[3];
    }else{
      nodes[0][0]=n0;
      nodes[0][1]=grid->equ[1];
      nodes[0][2]=grid->equ[3];
      nodes[0][3]=grid->equ[0];
      nodes[1][0]=n0;
      nodes[1][1]=grid->equ[3];
      nodes[1][2]=grid->equ[1];
      nodes[1][3]=grid->equ[2];
      nodes[2][0]=n1;
      nodes[2][1]=grid->equ[3];
      nodes[2][2]=grid->equ[1];
      nodes[2][3]=grid->equ[0];
      nodes[3][0]=n1;
      nodes[3][1]=grid->equ[1];
      nodes[3][2]=grid->equ[3];
      nodes[3][3]=grid->equ[2];
    }

    for ( i = 0 ; i < grid->ngem ; i++ ) 
      gridRemoveCell( grid, grid->gem[i] );
    
    for ( i = 0 ; i < 4 ; i++ )
      gridAddCell( grid, nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );

    return grid;
  }
  return NULL;
}

Grid *gridCycleEquator( Grid *grid )
{
  int i;

  for ( i = grid->nequ ; i > 0 ; i-- )
    grid->equ[i] = grid->equ[i-1];

  grid->equ[0] = grid->equ[grid->nequ];

  return grid;
}

Grid *gridSwapEdge5(Grid *grid, int n0, int n1 )
{
  int i;
  int currentindex, bestindex, nodes[6][4];

  double cost, origcost, currentcost, bestcost;

  origcost = 2.0;

  for ( i = 0 ; i < grid->ngem ; i++ ){
    cost = gridAR( grid, &grid->c2n[4*grid->gem[i]] );
    origcost = MIN(origcost,cost);
  }

  bestcost  =  -2.0;
  bestindex = -1;

  for ( currentindex = 0 ; currentindex < 5 ; currentindex++ ) {
    nodes[0][0]=n0;
    nodes[0][1]=grid->equ[0];
    nodes[0][2]=grid->equ[1];
    nodes[0][3]=grid->equ[2];
    nodes[1][0]=n0;
    nodes[1][1]=grid->equ[0];
    nodes[1][2]=grid->equ[2];
    nodes[1][3]=grid->equ[3];
    nodes[2][0]=n0;
    nodes[2][1]=grid->equ[0];
    nodes[2][2]=grid->equ[3];
    nodes[2][3]=grid->equ[4];
    nodes[3][0]=n1;
    nodes[3][1]=grid->equ[0];
    nodes[3][2]=grid->equ[2];
    nodes[3][3]=grid->equ[1];
    nodes[4][0]=n1;
    nodes[4][1]=grid->equ[0];
    nodes[4][2]=grid->equ[3];
    nodes[4][3]=grid->equ[2];
    nodes[5][0]=n1;
    nodes[5][1]=grid->equ[0];
    nodes[5][2]=grid->equ[4];
    nodes[5][3]=grid->equ[3];

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

  if (bestindex == -1 ) 
    printf("ERROR in bestindex, file %s line %d \n",__FILE__, __LINE__ ); 

  if ( bestcost > origcost ) {

    for ( i = 0 ; i < bestindex ; i++ ) 
      gridCycleEquator( grid );
    
    nodes[0][0]=n0;
    nodes[0][1]=grid->equ[0];
    nodes[0][2]=grid->equ[1];
    nodes[0][3]=grid->equ[2];
    nodes[1][0]=n0;
    nodes[1][1]=grid->equ[0];
    nodes[1][2]=grid->equ[2];
    nodes[1][3]=grid->equ[3];
    nodes[2][0]=n0;
    nodes[2][1]=grid->equ[0];
    nodes[2][2]=grid->equ[3];
    nodes[2][3]=grid->equ[4];
    nodes[3][0]=n1;
    nodes[3][1]=grid->equ[0];
    nodes[3][2]=grid->equ[2];
    nodes[3][3]=grid->equ[1];
    nodes[4][0]=n1;
    nodes[4][1]=grid->equ[0];
    nodes[4][2]=grid->equ[3];
    nodes[4][3]=grid->equ[2];
    nodes[5][0]=n1;
    nodes[5][1]=grid->equ[0];
    nodes[5][2]=grid->equ[4];
    nodes[5][3]=grid->equ[3];

    for ( i = 0 ; i < grid->ngem ; i++ ) 
      gridRemoveCell( grid, grid->gem[i] );
    
    for ( i = 0 ; i < 6 ; i++ )
      gridAddCell( grid, nodes[i][0], nodes[i][1], nodes[i][2], nodes[i][3] );
  return grid;
  }  
  return NULL;
}

Grid *gridGetCombo6( Grid *grid, int nodes[40][4], double costs[20], 
		     double *bestcost, int best[4] );
Grid *gridConstructTet6( Grid *grid, int n0, int n1, 
			 int nodes[40][4], double costs[20] );

Grid *gridSwapEdge6( Grid *grid, int n0, int n1 )
{
  int i;
  int nodes[40][4], bestcombo[4];

  double cost, costs[20], origcost, bestcost;
  
  origcost = 2.0;

  for ( i = 0 ; i < grid->ngem ; i++ ){
    cost = gridAR( grid, &grid->c2n[4*grid->gem[i]] );
    origcost = MIN(origcost,cost);
  }
  
  gridConstructTet6( grid, n0, n1, nodes, costs );
  
  gridGetCombo6( grid, nodes, costs, &bestcost, bestcombo );
  
  if ( bestcost > origcost ) {
      
    for ( i = 0 ; i < grid->ngem ; i++ ) 
      gridRemoveCell( grid, grid->gem[i] );
    
    for ( i = 0 ; i < 4 ; i++ ){
      gridAddCell( grid, 
		   nodes[bestcombo[i]][0], 
		   nodes[bestcombo[i]][1], 
		   nodes[bestcombo[i]][2], 
		   nodes[bestcombo[i]][3] );
      gridAddCell( grid, 
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

  for ( i = 0; i < 6; i++ )
    grid->equ[i+6]=grid->equ[i];

  /* make the small triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i   ][0]=grid->equ[i];
    nodes[i   ][1]=grid->equ[i+5];
    nodes[i   ][2]=grid->equ[i+1];
    nodes[i   ][3]=n0;
    nodes[i+20][0]=grid->equ[i+1];
    nodes[i+20][1]=grid->equ[i+5];
    nodes[i+20][2]=grid->equ[i];
    nodes[i+20][3]=n1;
  }

  /* make the next triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i+06][0]=grid->equ[i+5];
    nodes[i+06][1]=grid->equ[i+2];
    nodes[i+06][2]=grid->equ[i+1];
    nodes[i+06][3]=n0;
    nodes[i+26][0]=grid->equ[i+1];
    nodes[i+26][1]=grid->equ[i+2];
    nodes[i+26][2]=grid->equ[i+5];
    nodes[i+26][3]=n1;
  }

  /* make the previous triangles */
  for ( i = 0; i < 6; i++ ){
    nodes[i+12][0]=grid->equ[i+5];
    nodes[i+12][1]=grid->equ[i+4];
    nodes[i+12][2]=grid->equ[i+1];
    nodes[i+12][3]=n0;
    nodes[i+32][0]=grid->equ[i+1];
    nodes[i+32][1]=grid->equ[i+4];
    nodes[i+32][2]=grid->equ[i+5];
    nodes[i+32][3]=n1;
  }

  /* make the big triangles */
  for ( i = 0; i < 2; i++ ){
    nodes[i+18][0]=grid->equ[i+5];
    nodes[i+18][1]=grid->equ[i+3];
    nodes[i+18][2]=grid->equ[i+1];
    nodes[i+18][3]=n0;
    nodes[i+38][0]=grid->equ[i+1];
    nodes[i+38][1]=grid->equ[i+3];
    nodes[i+38][2]=grid->equ[i+5];
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

Grid *gridSwapEdge7( Grid *grid, int n0, int n1 )
{
  int i;
  int nodes[70][4], bestcombo[5];

  double cost, costs[35], origcost, bestcost;
  
  origcost = 2.0;

  for ( i = 0 ; i < grid->ngem ; i++ ){
    cost = gridAR( grid, &grid->c2n[4*grid->gem[i]] );
    origcost = MIN(origcost,cost);
  }

  gridConstructTet7( grid, n0, n1, nodes, costs );

  gridGetCombo7( grid, nodes, costs, &bestcost, bestcombo );
  
  if ( bestcost > origcost ) {

    for ( i = 0 ; i < grid->ngem ; i++ ) 
      gridRemoveCell( grid, grid->gem[i] );
      
    for ( i = 0 ; i < 5 ; i++ ){
      gridAddCell( grid, 
		   nodes[bestcombo[i]][0], 
		   nodes[bestcombo[i]][1], 
		   nodes[bestcombo[i]][2], 
		   nodes[bestcombo[i]][3] );
      gridAddCell( grid, 
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

  for ( i = 0; i < 7; i++ )
    grid->equ[i+7]=grid->equ[i];

  /* make the small triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i   ][0]=grid->equ[i];
    nodes[i   ][1]=grid->equ[i+6];
    nodes[i   ][2]=grid->equ[i+1];
    nodes[i   ][3]=n0;
    nodes[i+35][0]=grid->equ[i+1];
    nodes[i+35][1]=grid->equ[i+6];
    nodes[i+35][2]=grid->equ[i];
    nodes[i+35][3]=n1;
  }

  /* make the next triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+07][0]=grid->equ[i+6];
    nodes[i+07][1]=grid->equ[i+2];
    nodes[i+07][2]=grid->equ[i+1];
    nodes[i+07][3]=n0;
    nodes[i+42][0]=grid->equ[i+1];
    nodes[i+42][1]=grid->equ[i+2];
    nodes[i+42][2]=grid->equ[i+6];
    nodes[i+42][3]=n1;
  }

  /* make the previous triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+14][0]=grid->equ[i+6];
    nodes[i+14][1]=grid->equ[i+5];
    nodes[i+14][2]=grid->equ[i+1];
    nodes[i+14][3]=n0;
    nodes[i+49][0]=grid->equ[i+1];
    nodes[i+49][1]=grid->equ[i+5];
    nodes[i+49][2]=grid->equ[i+6];
    nodes[i+49][3]=n1;
  }

  /* make the spike */
  for ( i = 0; i < 7; i++ ){
    nodes[i+21][0]=grid->equ[i];
    nodes[i+21][1]=grid->equ[i+4];
    nodes[i+21][2]=grid->equ[i+3];
    nodes[i+21][3]=n0;
    nodes[i+56][0]=grid->equ[i+3];
    nodes[i+56][1]=grid->equ[i+4];
    nodes[i+56][2]=grid->equ[i];
    nodes[i+56][3]=n1;
  }

  /* make the big triangles */
  for ( i = 0; i < 7; i++ ){
    nodes[i+28][0]=grid->equ[i];
    nodes[i+28][1]=grid->equ[i+5];
    nodes[i+28][2]=grid->equ[i+2];
    nodes[i+28][3]=n0;
    nodes[i+63][0]=grid->equ[i+2];
    nodes[i+63][1]=grid->equ[i+5];
    nodes[i+63][2]=grid->equ[i];
    nodes[i+63][3]=n1;
  }

  for ( i = 0; i < 35; i++ )
    costs[i] = MIN( gridAR( grid, nodes[i] ) , gridAR( grid, nodes[i+35]) ); 
  
  return grid;
}

Grid *gridSwapCellFaceArea(Grid *grid, int cell)
{
  int nodes[4];
  int i, nface, emptyIndex, faceCheck[4], face[4], haveFace[4];
  double origArea, newArea;
  int faceNodes[4][3]= {{0,1,2},{0,3,1},{1,3,2},{2,3,0}};

  if ( NULL == gridCell( grid, cell, nodes) ) return NULL;

  for (i=0;i<4;i++){
    faceCheck[i] = 
      gridFindFace(grid,faceNodes[i][0],faceNodes[i][1],faceNodes[i][2]);
  }

  nface = 0;
  emptyIndex = 3;
  for (i=0;i<4;i++){
    if (EMPTY != faceCheck[i]) {
      face[nface] = faceCheck[i];
      haveFace[nface] = i;
      nface++;
    }else{
      haveFace[emptyIndex] = i;
      emptyIndex--;
    }
  }

  if ( 2 != nface ) return NULL;
  if ( grid->faceId[face[0]] != grid->faceId[face[1]] )return NULL;
  
  i = haveFace[0];
  origArea = 
    gridFaceArea(grid,faceNodes[i][0],faceNodes[i][1],faceNodes[i][2]);
  i = haveFace[1];
  origArea = origArea + 
    gridFaceArea(grid,faceNodes[i][0],faceNodes[i][1],faceNodes[i][2]);

  i = haveFace[2];
  newArea = 
    gridFaceArea(grid,faceNodes[i][0],faceNodes[i][1],faceNodes[i][2]);
  i = haveFace[3];
  newArea = newArea + 
    gridFaceArea(grid,faceNodes[i][0],faceNodes[i][1],faceNodes[i][2]);

  //printf("gridSwapCellFaceArea: %f\n",newArea/origArea);

  if (newArea>origArea) return NULL;

  return grid;
}

Grid *gridSwapFaceArea(Grid *grid)
{
  int cell;

  for (cell=0;cell<grid->maxcell;cell++) {
    if ( EMPTY != grid->c2n[4*cell] ) gridSwapCellFaceArea( grid, cell );
  }
  return grid;
}

