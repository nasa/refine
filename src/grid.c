
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
#include "grid.h"
#include "adj.h"

#define MAXDEG 200

struct Grid {
  int maxnode, nnode;
  double *xyz;

  int maxcell, ncell;
  int blankc2n;
  int *c2n;
  Adj *cellAdj;

  int maxface, nface;
  int blankf2n;
  int *f2n;
  int *faceId;
  Adj *faceAdj;

  int ngem;
  int gem[MAXDEG];

  int nequ;
  int equ[MAXDEG];
};

//#define EBUG

Grid* gridCreate(int maxnode, int maxcell, int maxface)
{
  int i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->maxnode = MAX(maxnode,1);
  grid->nnode   = 0;
  grid->maxcell = MAX(maxcell,1);
  grid->ncell   = 0;
  grid->maxface = MAX(maxface,1);
  grid->nface   = 0;

  grid->xyz = malloc(3 * grid->maxnode * sizeof(double));

  // cells
  grid->c2n = malloc(4 * grid->maxcell * sizeof(int));
  for (i=0;i < grid->maxcell; i++ ) {
    grid->c2n[0+4*i] = EMPTY; 
    grid->c2n[1+4*i] = i+1; 
  }
  grid->c2n[1+4*(grid->maxcell-1)] = EMPTY; 
  grid->blankc2n = 0;

  grid->cellAdj = adjCreate(grid->maxnode,grid->maxcell*4);

  // face
  grid->f2n    = malloc(3 * grid->maxface * sizeof(int));
  grid->faceId = malloc(1 * grid->maxface * sizeof(int));
  for (i=0;i < grid->maxface; i++ ) {
    grid->f2n[0+3*i] = EMPTY; 
    grid->f2n[1+3*i] = i+1; 
    grid->faceId[i] = EMPTY; 
  }
  grid->f2n[1+3*(grid->maxface-1)] = EMPTY; 
  grid->blankf2n = 0;

  grid->faceAdj = adjCreate(grid->maxnode,grid->maxface*3);

  grid->ngem = 0;

  return  grid;
}

Grid *gridImport(int nnode, int nface, int maxcell, int ncell,
		 double *xyz, int *f2n, int *faceId, int *c2n )
{
  int i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->maxnode = nnode;
  grid->nnode   = nnode;
  grid->maxcell = maxcell;
  grid->ncell   = ncell;
  grid->maxface = nface;
  grid->nface   = nface;

  printf("gridImport: %d nodes %d faces %d cells\n",
	 grid->nnode,grid->nface,grid->ncell);

  grid->xyz = xyz;

  // cells
  grid->c2n = c2n;
  for ( i=grid->ncell ; i < grid->maxcell ; i++ ) {
    grid->c2n[0+4*i] = EMPTY; 
    grid->c2n[1+4*i] = i+1; 
  }
  if (grid->maxcell != grid->ncell) grid->c2n[1+4*(grid->maxcell-1)] = EMPTY; 
  grid->blankc2n = grid->ncell;

  grid->cellAdj = adjCreate(grid->maxnode,grid->maxcell*4);

  for ( i=0 ; i < grid->ncell ; i++ ) {
    adjRegister(grid->cellAdj,grid->c2n[0+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[1+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[2+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[3+4*i],i);
  }

  grid->f2n    = f2n;
  grid->faceId = faceId;
  grid->blankf2n = EMPTY;

  grid->faceAdj = adjCreate(grid->maxnode,grid->maxface*3);

  //addface

  for ( i=0 ; i < grid->nface ; i++ ) {
    adjRegister(grid->faceAdj,grid->f2n[0+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[1+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[2+3*i],i);
  }

  grid->ngem = 0;

  return  grid;
}

void gridFree(Grid *grid)
{
  adjFree(grid->faceAdj);
  free(grid->faceId);
  free(grid->f2n);
  adjFree(grid->cellAdj);
  free(grid->c2n);
  free(grid->xyz);
  free(grid);
}

int gridMaxNode(Grid *grid)
{
  return grid->maxnode;
}

int gridNNode(Grid *grid)
{
  return grid->nnode;
}

int gridMaxCell(Grid *grid)
{
  return grid->maxcell;
}

int gridNCell(Grid *grid)
{
  return grid->ncell;
}

int gridMaxFace(Grid *grid)
{
  return grid->maxface;
}

int gridNFace(Grid *grid)
{
  return grid->nface;
}

int gridNGem(Grid *grid)
{
  return grid->ngem;
}

int gridGem(Grid *grid, int index)
{
  return grid->gem[index];
}

int gridNEqu(Grid *grid)
{
  return grid->nequ;
}

int gridEqu(Grid *grid, int index)
{
  return grid->equ[index];
}

int gridCellDegree(Grid *grid, int id)
{
  return adjDegree(grid->cellAdj, id);
}

Grid *gridAddCell(Grid *grid, int n0, int n1, int n2, int n3)
{
  int cellId,icell;
  if ( grid->blankc2n == EMPTY ) return NULL;
  cellId = grid->blankc2n;
  grid->blankc2n = grid->c2n[1+4*cellId];
  grid->ncell++;
  
  grid->c2n[0+4*cellId] = n0;
  grid->c2n[1+4*cellId] = n1;
  grid->c2n[2+4*cellId] = n2;
  grid->c2n[3+4*cellId] = n3;
  
  if ( NULL == adjRegister( grid->cellAdj, n0, cellId ) ) return NULL;
  if ( NULL == adjRegister( grid->cellAdj, n1, cellId ) ) return NULL;
  if ( NULL == adjRegister( grid->cellAdj, n2, cellId ) ) return NULL;
  if ( NULL == adjRegister( grid->cellAdj, n3, cellId ) ) return NULL;
  
  return grid;
}

Grid *gridRemoveCell(Grid *grid, int cellId )
{
  if (cellId >= grid->maxcell || cellId < 0 ) return NULL;
  if (EMPTY == grid->c2n[4*cellId]) return NULL;
  
  if ( grid->ncell <= 0) return NULL;
  grid->ncell--;

  if( ( NULL == adjRemove( grid->cellAdj, grid->c2n[0+4*cellId], cellId ) ) || 
      ( NULL == adjRemove( grid->cellAdj, grid->c2n[1+4*cellId], cellId ) ) ||
      ( NULL == adjRemove( grid->cellAdj, grid->c2n[2+4*cellId], cellId ) ) || 
      ( NULL == adjRemove( grid->cellAdj, grid->c2n[3+4*cellId], cellId ) ) ) 
    return NULL;  

  grid->c2n[0+4*cellId] = EMPTY;
  grid->c2n[1+4*cellId] = grid->blankc2n;
  grid->blankc2n = cellId;

  return grid;
}

Grid *gridCell(Grid *grid, int cellId, int *nodes )
{
  if ( cellId >= grid->maxcell ) return NULL;
  if ( grid->c2n[4*cellId] == EMPTY ) return NULL;

  nodes[0] = grid->c2n[0+4*cellId];
  nodes[1] = grid->c2n[1+4*cellId];
  nodes[2] = grid->c2n[2+4*cellId];
  nodes[3] = grid->c2n[3+4*cellId];

  return grid;
}

Grid *gridAddFace(Grid *grid, int n0, int n1, int n2, int faceId )
{
  int face;
  face = grid->nface;
  if (grid->nface >= grid->maxface) return NULL;
  grid->nface++;

  grid->f2n[0+3*face] = n0;
  grid->f2n[1+3*face] = n1;
  grid->f2n[2+3*face] = n2;
  grid->faceId[face]  = faceId;

  if ( NULL == adjRegister( grid->faceAdj, n0, face ) ) return NULL;
  if ( NULL == adjRegister( grid->faceAdj, n1, face ) ) return NULL;
  if ( NULL == adjRegister( grid->faceAdj, n2, face ) ) return NULL;

 return grid;
}

Grid *gridRemoveFace(Grid *grid, int face )
{
  if (face >= grid->maxface || face < 0) return NULL;
  if (EMPTY == grid->f2n[3*face]) return NULL;
  
  if ( grid->nface <= 0) return NULL;
  grid->nface--;

  if( ( NULL == adjRemove( grid->faceAdj, grid->f2n[0+3*face], face ) ) || 
      ( NULL == adjRemove( grid->faceAdj, grid->f2n[1+3*face], face ) ) || 
      ( NULL == adjRemove( grid->faceAdj, grid->f2n[2+3*face], face ) ) )
    return NULL;  

  grid->f2n[0+3*face] = EMPTY;
  grid->f2n[1+3*face] = grid->blankf2n;
  grid->blankf2n = face;

  return grid;
}

int gridFindFace(Grid *grid, int n0, int n1, int n2 )
{
  AdjIterator it0, it1, it2;
  Adj *adj=grid->faceAdj;

  for ( it0 = adjFirst(adj,n0); adjValid(it0); it0 = adjNext(it0) )
    for ( it1 = adjFirst(adj,n1); adjValid(it1); it1 = adjNext(it1) )
      if ( adjItem(it0) == adjItem(it1) )
	for ( it2 = adjFirst(adj,n2); adjValid(it2); it2 = adjNext(it2) )
	  if ( adjItem(it0)==adjItem(it2) ) return adjItem(it2);

  return EMPTY;
}

int gridFaceId(Grid *grid, int n0, int n1, int n2 )
{
  int face = gridFindFace(grid, n0, n1, n2 );
  if ( face == EMPTY ) return EMPTY;
  return grid->faceId[face];
}

Grid *gridMakeGem(Grid *grid, int n0, int n1 )
{
  AdjIterator it;
  int cellId;
  grid->ngem = 0;

  for ( it = adjFirst(grid->cellAdj,n0); 
	adjValid(it); 
	it = adjNext(it)) {

    cellId = adjItem(it);
    if ( n1 == grid->c2n[0+4*cellId] ||
	 n1 == grid->c2n[1+4*cellId] ||
	 n1 == grid->c2n[2+4*cellId] ||
	 n1 == grid->c2n[3+4*cellId] ) {
      if (grid->ngem >= MAXDEG) { 
	grid->ngem = 0; 
	return NULL; 
      }
      grid->gem[grid->ngem] = cellId;
      grid->ngem++;
    }
  }

  return grid;
}

Grid *gridOrient(Grid *grid, int *c, int *n )
{

  /* 0 leads */
  if ( n[0] == c[0] && n[1] == c[1] ){ 
    n[2] = c[2]; 
    n[3] = c[3]; 
    return grid; 
  }
  if ( n[0] == c[0] && n[1] == c[2] ){ 
    n[2] = c[3]; 
    n[3] = c[1]; 
    return grid; 
  }
  if ( n[0] == c[0] && n[1] == c[3] ) {
    n[2] = c[1];
    n[3] = c[2];
    return grid;
  }

  /* 1 leads */
  if ( n[0] == c[1] && n[1] == c[0] ) {
    n[2] = c[3];
    n[3] = c[2];
    return grid;
  }
  if ( n[0] == c[2] && n[1] == c[0] ) {
    n[2] = c[1];
    n[3] = c[3];
    return grid;
  }
  if ( n[0] == c[3] && n[1] == c[0] ) {
    n[2] = c[2];
    n[3] = c[1];
    return grid;
  }

  /* 2 leads */
  if ( n[0] == c[2] && n[1] == c[3] ) {
    n[2] = c[0];
    n[3] = c[1];
    return grid;
  }
  if ( n[0] == c[3] && n[1] == c[1] ) {
    n[2] = c[0];
    n[3] = c[2];
    return grid;
  }
  if ( n[0] == c[1] && n[1] == c[2] ) {
    n[2] = c[0];
    n[3] = c[3];
    return grid;
  }

  /* 3 leads */
  if ( n[0] == c[3] && n[1] == c[2] ) {
    n[2] = c[1];
    n[3] = c[0];
    return grid;
  }
  if ( n[0] == c[1] && n[1] == c[3] ) {
    n[2] = c[2];
    n[3] = c[0];
    return grid;
  }
  if ( n[0] == c[2] && n[1] == c[1] ) {
    n[2] = c[3];
    n[3] = c[0];
    return grid;
  }

  return NULL;
}

Grid *gridEquator(Grid *grid, int n0, int n1 )
{
  int igem, iequ, cell[4], nodes[4];
  bool gap, found;
  grid->nequ = 0;

  if ( NULL == gridMakeGem( grid, n0, n1) ) return NULL;

  if ( grid->ngem == (MAXDEG-1) ) return NULL;
  if ( grid->ngem == 0 ) return grid;

  nodes[0] = n0;
  nodes[1] = n1;

  gap = FALSE;

  gridOrient( grid, &grid->c2n[4*grid->gem[0]], nodes );

  grid->equ[0]= nodes[3];

  // put and not found in loops

  for ( iequ=1 ; iequ <= grid->ngem && !gap ; iequ++ ){
    found = FALSE;
    for( igem=0 ; igem<grid->ngem ; igem++ ){
      gridOrient( grid, &grid->c2n[4*grid->gem[igem]], nodes );
      if ( grid->equ[iequ-1] == nodes[2] ) {
	grid->equ[iequ] = nodes[3];
	found = TRUE;
      }
    }
    if (!found) {
      gap = TRUE;
      grid->equ[grid->ngem] = grid->equ[iequ-1];
    }
  }
  
  if ( gap ) {

    for ( iequ=grid->ngem-1 ; iequ >= 0 ; iequ-- ){
      found = FALSE;
      for( igem=0 ; igem<grid->ngem ; igem++ ){
	gridOrient( grid, &grid->c2n[4*grid->gem[igem]], nodes );
	if ( grid->equ[iequ+1] == nodes[3] ) {
	  grid->equ[iequ] = nodes[2];
	  found = TRUE;
	}
      }
      if (!found) return NULL;
    }
    grid->nequ = grid->ngem+1;
    grid->equ[grid->ngem+1] = grid->equ[0];

  }else{
    grid->nequ = grid->ngem;
  }

  return grid;
}

Grid *gridSwapEdge4(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge5(Grid *grid, int n0, int n1 );
Grid *gridSwapEdge6(Grid *grid, int n0, int n1 );

Grid *gridSwapEdge(Grid *grid, int n0, int n1 )
{
  int gap0, gap1, face0, face1, faceId0, faceId1, newFaceId0, newFaceId1;
  Grid *swapStatus;

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
  }

  if (grid->nequ==4) swapStatus = gridSwapEdge4(grid, n0, n1);
  if (grid->nequ==5) swapStatus = gridSwapEdge5(grid, n0, n1);
  if (grid->nequ==6) swapStatus = gridSwapEdge6(grid, n0, n1);

  if ( grid->nequ != grid->ngem && swapStatus != NULL ) {
    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );
    gridAddFace(grid, n0, gap0, gap1, faceId0 );
    gridAddFace(grid, n1, gap0, gap1, faceId0 );
  }

  return swapStatus;
}

Grid *gridSwap(Grid *grid)
{
  int cellId, nodes[4];
  for (cellId=0;cellId<grid->maxcell;cellId++){
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
  return grid;
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
  }
  
  return grid;
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
  }
  
  return( grid );
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
  }
  
  return grid;
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

int gridAddNode(Grid *grid, double x, double y, double z )
{
  int nodeId;
  if (grid->nnode >= grid->maxnode) return EMPTY;
  nodeId = grid->nnode;
  grid->nnode++;

  grid->xyz[0+3*nodeId] = x;
  grid->xyz[1+3*nodeId] = y;
  grid->xyz[2+3*nodeId] = z;

  return nodeId;
}

double gridVolume(Grid *grid, int *nodes )
{
  int ixyz;
  double edge1[3], edge2[3], edge3[3], norm[3], volume; 
  
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){

    edge1[ixyz] = grid->xyz[ixyz+3*nodes[1]]
                - grid->xyz[ixyz+3*nodes[0]];
    edge2[ixyz] = grid->xyz[ixyz+3*nodes[2]]
                - grid->xyz[ixyz+3*nodes[0]];
    edge3[ixyz] = grid->xyz[ixyz+3*nodes[3]]
                - grid->xyz[ixyz+3*nodes[0]];
      
  }

  norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; 
  norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; 
  norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 

  return  (norm[0]*edge3[0]+norm[1]*edge3[1]+norm[2]*edge3[2])/6.0;

}

double gridAR(Grid *grid, int *nodes )
{
  double x1, x2, x3, x4; 
  double y1, y2, y3, y4; 
  double z1, z2, z3, z4; 
  double s1, s2, s3, s4, det;
  double xr, yr, zr;
  double circ;
  double nx1, ny1, nz1, rmag1;
  double nx2, ny2, nz2, rmag2;
  double nx3, ny3, nz3, rmag3;
  double nx4, ny4, nz4, rmag4;
  double xins;
  double aspect, cost;

  x1 = grid->xyz[0+3*nodes[0]];
  y1 = grid->xyz[1+3*nodes[0]];
  z1 = grid->xyz[2+3*nodes[0]];

  x2 = grid->xyz[0+3*nodes[1]];
  y2 = grid->xyz[1+3*nodes[1]];
  z2 = grid->xyz[2+3*nodes[1]];

  x3 = grid->xyz[0+3*nodes[2]];
  y3 = grid->xyz[1+3*nodes[2]];
  z3 = grid->xyz[2+3*nodes[2]];

  x4 = grid->xyz[0+3*nodes[3]];
  y4 = grid->xyz[1+3*nodes[3]];
  z4 = grid->xyz[2+3*nodes[3]];

  /* Compute the aspect ratios */

  det = (x4-x1)*((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
      + (y4-y1)*((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
      + (z4-z1)*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
  s1 = ((x2*x2 + y2*y2 + z2*z2) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  s2 = ((x3*x3 + y3*y3 + z3*z3) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  s3 = ((x4*x4 + y4*y4 + z4*z4) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  xr  =(s3     *((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
       	+ (y4-y1)*(s2     *(z2-z1) - s1     *(z3-z1))
	+ (z4-z1)*(s1     *(y3-y1) - s2     *(y2-y1)))/det;
  yr  =((x4-x1)*(s1     *(z3-z1) - s2     *(z2-z1))
	+ s3     *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
	+ (z4-z1)*((x2-x1)*s2      - (x3-x1)*s1     ))/det;
  zr  =((x4-x1)*((y2-y1)*s2      - (y3-y1)*s1     )
	+ (y4-y1)*((x3-x1)*s1      - (x2-x1)*s2     )
	+ s3     *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)))/det;
  circ = sqrt((x1-xr)*(x1-xr) + (y1-yr)*(y1-yr) + (z1-zr)*(z1-zr));

  /* Get the in-circle */

  nx1 = (y3 - y1)*(z2 - z1) - (y2 - y1)*(z3 - z1);
  ny1 =-(x3 - x1)*(z2 - z1) + (x2 - x1)*(z3 - z1);
  nz1 = (x3 - x1)*(y2 - y1) - (x2 - x1)*(y3 - y1);
  rmag1 = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);

  nx2 = (y2 - y1)*(z4 - z1) - (y4 - y1)*(z2 - z1);
  ny2 =-(x2 - x1)*(z4 - z1) + (x4 - x1)*(z2 - z1);
  nz2 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);
  rmag2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);

  nx3 = (y4 - y1)*(z3 - z1) - (y3 - y1)*(z4 - z1);
  ny3 =-(x4 - x1)*(z3 - z1) + (x3 - x1)*(z4 - z1);
  nz3 = (x4 - x1)*(y3 - y1) - (x3 - x1)*(y4 - y1);
  rmag3 = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3);

  nx4 = (y3 - y2)*(z4 - z2) - (y4 - y2)*(z3 - z2);
  ny4 =-(x3 - x2)*(z4 - z2) + (x4 - x2)*(z3 - z2);
  nz4 = (x3 - x2)*(y4 - y2) - (x4 - x2)*(y3 - y2);
  rmag4 = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4);
  nx1 = nx1/rmag1;
  ny1 = ny1/rmag1;
  nz1 = nz1/rmag1;
  nx2 = nx2/rmag2;
  ny2 = ny2/rmag2;
  nz2 = nz2/rmag2;
  nx3 = nx3/rmag3;
  ny3 = ny3/rmag3;
  nz3 = nz3/rmag3;
  nx4 = nx4/rmag4;
  ny4 = ny4/rmag4;
  nz4 = nz4/rmag4;
  det= -(nx3*ny2*nz1) + nx4*ny2*nz1 + nx2*ny3*nz1 - nx4*ny3*nz1
    -nx2*ny4*nz1 + nx3*ny4*nz1 + nx3*ny1*nz2 - nx4*ny1*nz2
    -nx1*ny3*nz2 + nx4*ny3*nz2 + nx1*ny4*nz2 - nx3*ny4*nz2
    -nx2*ny1*nz3 + nx4*ny1*nz3 + nx1*ny2*nz3 - nx4*ny2*nz3
    -nx1*ny4*nz3 + nx2*ny4*nz3 + nx2*ny1*nz4 - nx3*ny1*nz4
    -nx1*ny2*nz4 + nx3*ny2*nz4 + nx1*ny3*nz4 - nx2*ny3*nz4;
  s1 = nx1*x1 + ny1*y1 + nz1*z1;
  s2 = nx2*x1 + ny2*y1 + nz2*z1;
  s3 = nx3*x1 + ny3*y1 + nz3*z1;
  s4 = nx4*x4 + ny4*y4 + nz4*z4;
  xins = (nx4*ny3*nz2*s1 - nx3*ny4*nz2*s1 - nx4*ny2*nz3*s1 +
	  nx2*ny4*nz3*s1 +
	  nx3*ny2*nz4*s1 - nx2*ny3*nz4*s1 - nx4*ny3*nz1*s2 +
	  nx3*ny4*nz1*s2 +
	  nx4*ny1*nz3*s2 - nx1*ny4*nz3*s2 - nx3*ny1*nz4*s2 +
	  nx1*ny3*nz4*s2 +
	  nx4*ny2*nz1*s3 - nx2*ny4*nz1*s3 - nx4*ny1*nz2*s3 +
	  nx1*ny4*nz2*s3 +
	  nx2*ny1*nz4*s3 - nx1*ny2*nz4*s3 - nx3*ny2*nz1*s4 +
	  nx2*ny3*nz1*s4 +
	  nx3*ny1*nz2*s4 - nx1*ny3*nz2*s4 - nx2*ny1*nz3*s4 +
	  nx1*ny2*nz3*s4)/det;

  aspect = xins/circ*3.0;

  if ( gridVolume( grid, nodes ) <= 1.0e-14) aspect = -1.0;
  return aspect;
}

double gridMinVolume( Grid *grid )
{
  int cellId, nodes[4];
  double minVol;
  minVol = 999.0;
  for (cellId=0;cellId<grid->maxcell;cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      minVol = MIN(minVol,gridVolume(grid, nodes) );
  return minVol;
}

double gridMinAR( Grid *grid )
{
  int cellId, nodes[4];
  double minAR;
  minAR = 999.0;
  for (cellId=0;cellId<grid->maxcell;cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      minAR = MIN(minAR, gridAR(grid, nodes) );
  return minAR;
}
