
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
#include "grid.h"
#include "adj.h"
#include "gridStruct.h"

//#define EBUG

Grid* gridCreate(int maxnode, int maxcell, int maxface, int maxedge)
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
  grid->maxedge = MAX(maxedge,1);
  grid->nedge   = 0;
  grid->nGeomNode   = 0;

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
  grid->faceU  = malloc(3 * grid->maxface * sizeof(double));
  grid->faceV  = malloc(3 * grid->maxface * sizeof(double));
  grid->faceId = malloc(1 * grid->maxface * sizeof(int));

  for (i=0;i < grid->maxface; i++ ) {
    grid->f2n[0+3*i] = EMPTY; 
    grid->f2n[1+3*i] = i+1; 
    grid->faceId[i] = EMPTY; 
    grid->faceU[0+3*i] = DBL_MAX;
    grid->faceU[1+3*i] = DBL_MAX;
    grid->faceU[2+3*i] = DBL_MAX;
    grid->faceV[0+3*i] = DBL_MAX;
    grid->faceV[1+3*i] = DBL_MAX;
    grid->faceV[2+3*i] = DBL_MAX;
  }
  grid->f2n[1+3*(grid->maxface-1)] = EMPTY; 
  grid->blankf2n = 0;

  grid->faceAdj = adjCreate(grid->maxnode,grid->maxface*3);

  // edge
  grid->e2n    = malloc(2 * grid->maxedge * sizeof(int));
  grid->edgeId = malloc(1 * grid->maxedge * sizeof(int));
  grid->edgeT  = malloc(2 * grid->maxedge * sizeof(double));
  for (i=0;i < grid->maxedge; i++ ) {
    grid->e2n[0+2*i] = EMPTY; 
    grid->e2n[1+2*i] = i+1; 
    grid->edgeId[i] = EMPTY; 
    grid->edgeT[0+2*i] = DBL_MAX; 
    grid->edgeT[1+2*i] = DBL_MAX; 
  }
  grid->e2n[1+2*(grid->maxedge-1)] = EMPTY; 
  grid->blanke2n = 0;

  grid->edgeAdj = adjCreate(grid->maxnode,grid->maxedge*2);

  grid->ngem = 0;

  return grid;
}

Grid *gridImport(int maxnode, int nnode, 
		 int maxface, int nface, 
		 int maxcell, int ncell,
		 int maxedge, 
		 double *xyz, int *f2n, int *faceId, int *c2n )
{
  int i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->maxnode = MAX(maxnode,1);
  grid->nnode   = nnode;
  grid->maxcell = MAX(maxcell,1);
  grid->ncell   = ncell;
  grid->maxface = MAX(maxface,1);
  grid->nface   = nface;
  grid->maxedge = MAX(maxedge,1);
  grid->nedge   = 0;
  grid-> nGeomNode = 0;

  grid->xyz = xyz;

  // cells
  grid->c2n = c2n;
  for ( i=grid->ncell ; i < grid->maxcell ; i++ ) {
    grid->c2n[0+4*i] = EMPTY; 
    grid->c2n[1+4*i] = i+1; 
  }
  if (grid->maxcell == grid->ncell) {
    grid->blankc2n = EMPTY;
  }else{
    grid->c2n[1+4*(grid->maxcell-1)] = EMPTY; 
    grid->blankc2n = grid->ncell;
  }
  grid->cellAdj = adjCreate(grid->maxnode,grid->maxcell*4);

  for ( i=0 ; i < grid->ncell ; i++ ) {
    adjRegister(grid->cellAdj,grid->c2n[0+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[1+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[2+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[3+4*i],i);
  }

  grid->f2n    = f2n;
  for ( i=grid->nface ; i < grid->maxface ; i++ ) {
    grid->f2n[0+3*i] = EMPTY; 
    grid->f2n[1+3*i] = i+1; 
  }
  if (grid->maxface == grid->nface) {
    grid->blankf2n = EMPTY;
  }else{
    grid->f2n[1+3*(grid->maxface-1)] = EMPTY; 
    grid->blankf2n = grid->nface;
  }
  grid->faceU  = malloc(3 * grid->maxface * sizeof(double));
  grid->faceV  = malloc(3 * grid->maxface * sizeof(double));
  for (i=0;i < grid->maxface; i++ ) {
    grid->faceU[0+3*i] = DBL_MAX;
    grid->faceU[1+3*i] = DBL_MAX;
    grid->faceU[2+3*i] = DBL_MAX;
    grid->faceV[0+3*i] = DBL_MAX;
    grid->faceV[1+3*i] = DBL_MAX;
    grid->faceV[2+3*i] = DBL_MAX;
  }

  grid->faceId = faceId;

  grid->faceAdj = adjCreate(grid->maxnode,grid->maxface*3);

  //addface

  for ( i=0 ; i < grid->nface ; i++ ) {
    adjRegister(grid->faceAdj,grid->f2n[0+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[1+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[2+3*i],i);
  }

  // edge
  grid->e2n    = malloc(2 * grid->maxedge * sizeof(int));
  grid->edgeId = malloc(1 * grid->maxedge * sizeof(int));
  grid->edgeT  = malloc(2 * grid->maxedge * sizeof(double));
  for (i=0;i < grid->maxedge; i++ ) {
    grid->e2n[0+2*i] = EMPTY; 
    grid->e2n[1+2*i] = i+1; 
    grid->edgeId[i] = EMPTY; 
    grid->edgeT[0+2*i] = DBL_MAX; 
    grid->edgeT[1+2*i] = DBL_MAX; 
  }
  grid->e2n[1+2*(grid->maxedge-1)] = EMPTY; 
  grid->blanke2n = 0;

  grid->edgeAdj = adjCreate(grid->maxnode,grid->maxedge*2);

  grid->ngem = 0;

  return  grid;
}

Grid *gridImportFAST( char *filename )
{
  FILE *file;
  int i, nnode, nface, maxcell, ncell;
  double *xyz;
  int *f2n, *faceId;
  int *c2n;

  file = fopen(filename,"r");
  fscanf(file,"%d %d %d",&nnode,&nface,&ncell);
  printf("fast size: %d nodes %d faces %d cells.\n",nnode,nface,ncell);

  printf("reading xyz...\n");
  
  xyz = malloc(3*nnode*sizeof(double));

  for( i=0; i<nnode ; i++ ) fscanf(file,"%lf",&xyz[0+3*i]);
  for( i=0; i<nnode ; i++ ) fscanf(file,"%lf",&xyz[1+3*i]);
  for( i=0; i<nnode ; i++ ) fscanf(file,"%lf",&xyz[2+3*i]);

  printf("reading faces...\n");
  
  f2n = malloc(3*nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    fscanf(file,"%d",&f2n[0+3*i]);
    fscanf(file,"%d",&f2n[1+3*i]);
    fscanf(file,"%d",&f2n[2+3*i]);
    f2n[0+3*i]--;
    f2n[1+3*i]--;
    f2n[2+3*i]--;
  }

  printf("reading face ID tags...\n");
  
  faceId = malloc(nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    fscanf(file,"%d",&faceId[i]);
  }

  printf("reading cells...\n");
  
  maxcell = ncell*2;

  c2n = malloc(4*maxcell*sizeof(int));

  for( i=0; i<ncell ; i++ ) {
    fscanf(file,"%d",&c2n[0+4*i]);
    fscanf(file,"%d",&c2n[1+4*i]);
    fscanf(file,"%d",&c2n[2+4*i]);
    fscanf(file,"%d",&c2n[3+4*i]);
    c2n[0+4*i]--;
    c2n[1+4*i]--;
    c2n[2+4*i]--;
    c2n[3+4*i]--;
  }

  fclose(file);

  return gridImport( nnode, nnode, nface, nface, maxcell, ncell, 0,
		     xyz, f2n, faceId, c2n );
}

Grid *gridExport(Grid *grid, int *nnode, int *nface, int *ncell,
		 double **xyz, int **f2n, int **faceId, int **c2n )
{
  int i, origcell, packcell, origface, packface;
  int iface, n0, n1, n2, id;
  bool emptyFace;

  printf("gridExport: %d nodes %d faces %d cells\n",
	 grid->nnode,grid->nface,grid->ncell);

  *nnode = grid->nnode;
  *ncell = grid->ncell;
  *nface = grid->nface;

  *xyz = grid->xyz;

  packcell =0;
  for ( origcell=0 ; origcell < grid->maxcell ; origcell++ )
    if ( grid->c2n[0+4*origcell] != EMPTY) {
      grid->c2n[0+4*packcell] = grid->c2n[0+4*origcell];
      grid->c2n[1+4*packcell] = grid->c2n[1+4*origcell];
      grid->c2n[2+4*packcell] = grid->c2n[2+4*origcell];
      grid->c2n[3+4*packcell] = grid->c2n[3+4*origcell];
      packcell++;
    } 
  
  if (grid->ncell != packcell) {
    printf("ERROR: grid->ncell != packcell, file %s line %d \n",
	   __FILE__, __LINE__ );
    return NULL;
  }
  *c2n = grid->c2n;

  packface=0;

  emptyFace = FALSE;
  iface=0;
  while (!emptyFace) {
    emptyFace = TRUE;
    iface++;
    for ( origface=0 ; origface < grid->maxface ; origface++ ){ 
      if ( grid->faceId[origface]==iface && grid->f2n[0+3*origface] != EMPTY) {
	emptyFace = FALSE;
	n0 = grid->f2n[0+3*origface];
	n1 = grid->f2n[1+3*origface];
	n2 = grid->f2n[2+3*origface];
	id = grid->faceId[origface];

	grid->f2n[0+3*origface]	= grid->f2n[0+3*packface];
	grid->f2n[1+3*origface]	= grid->f2n[1+3*packface];
	grid->f2n[2+3*origface]	= grid->f2n[2+3*packface];
	grid->faceId[origface]	= grid->faceId[packface];

	grid->f2n[0+3*packface] = n0;
	grid->f2n[1+3*packface] = n1;
	grid->f2n[2+3*packface] = n2;
	grid->faceId[packface]  = id;
	packface++;
      } 
    }
  }
  printf("gridExport: %d geometry faces detected.\n",iface);

  if (grid->nface != packface) {
    printf("ERROR: grid->nface %d != packface %d, file %s line %d \n",
	   grid->nface, packface, __FILE__, __LINE__ );
    return NULL;
  }

  *f2n    = grid->f2n;
  *faceId = grid->faceId;

  return  grid;
}

void gridFree(Grid *grid)
{
  adjFree(grid->edgeAdj);
  free(grid->edgeId);
  free(grid->edgeT);
  free(grid->e2n);
  adjFree(grid->faceAdj);
  free(grid->faceId);
  free(grid->faceV);
  free(grid->faceU);
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

int gridMaxEdge(Grid *grid)
{
  return grid->maxedge;
}

int gridNEdge(Grid *grid)
{
  return grid->nedge;
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
  int cellId;
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
 return gridAddFaceUV(grid, 
		      n0, DBL_MAX, DBL_MAX,
		      n1, DBL_MAX, DBL_MAX,
		      n2, DBL_MAX, DBL_MAX,
		      faceId );
}

Grid *gridAddFaceUV(Grid *grid, 
		    int n0, double u0, double v0,
		    int n1, double u1, double v1,
		    int n2, double u2, double v2, int faceId )
{
  int face;
  if ( grid->blankf2n == EMPTY ) return NULL;
  face = grid->blankf2n;
  grid->blankf2n = grid->f2n[1+3*face];
  grid->nface++;

  grid->f2n[0+3*face] = n0;
  grid->f2n[1+3*face] = n1;
  grid->f2n[2+3*face] = n2;
  grid->faceU[0+3*face] = u0;
  grid->faceU[1+3*face] = u1;
  grid->faceU[2+3*face] = u2;
  grid->faceV[0+3*face] = v0;
  grid->faceV[1+3*face] = v1;
  grid->faceV[2+3*face] = v2;
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

Grid *gridNodeUV(Grid *grid, int  node, int faceId, double *uv )
{
  AdjIterator it;
  int face, i;

  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    if ( grid->faceId[face] == faceId ) {
      for ( i=0 ; i<3 ; i++ ) {
	if (grid->f2n[i+3*face] == node){
	  uv[0] = grid->faceU[i+3*face];
	  uv[1] = grid->faceV[i+3*face];
	  return grid;
	}
      }
    }
  }

  uv[0] = DBL_MAX;
  uv[1] = DBL_MAX;
  return NULL;
}

Grid *gridSetNodeUV(Grid *grid, int  node, int faceId, double u, double v )
{
  AdjIterator it;
  int face, i;
  bool found;
  found = FALSE;

  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    if ( grid->faceId[face] == faceId ) {
      for ( i=0 ; i<3 ; i++ ) {
	if (grid->f2n[i+3*face] == node){
	  grid->faceU[i+3*face] = u;
	  grid->faceV[i+3*face] = v;
	  found = TRUE;
	}
      }
    }
  }

  if (found) return grid;
  return NULL;
}

double gridNodeU(Grid *grid, int node, int faceId)
{
  double uv[2];
  if ( grid != gridNodeUV(grid, node, faceId, uv ) ) return DBL_MAX;
  return uv[0];
}

double gridNodeV(Grid *grid, int node, int faceId)
{
  double uv[2];
  if ( grid != gridNodeUV(grid, node, faceId, uv ) ) return DBL_MAX;
  return uv[1];
}

Grid *gridNodeT(Grid *grid, int  node, int edgeId, double *t )
{
  AdjIterator it;
  int edge;

  for ( it = adjFirst(grid->edgeAdj,node); adjValid(it); it = adjNext(it) ){
    edge = adjItem(it);
    if ( grid->edgeId[edge] == edgeId ) {
      if (grid->e2n[0+2*edge] == node){
	*t = grid->edgeT[0+2*edge];
	return grid;
      }
      if (grid->e2n[1+2*edge] == node){
	*t = grid->edgeT[1+2*edge];
	return grid;
      }
    }
  }

  *t = DBL_MAX;
  return NULL;
}

Grid *gridSetNodeT(Grid *grid, int  node, int edgeId, double t )
{
  AdjIterator it;
  int edge;
  bool found;
  found = FALSE;

  for ( it = adjFirst(grid->edgeAdj,node); adjValid(it); it = adjNext(it) ){
    edge = adjItem(it);
    if ( grid->edgeId[edge] == edgeId ) {
      if (grid->e2n[0+2*edge] == node){
	grid->edgeT[0+2*edge] = t;
	found = TRUE;
      }
      if (grid->e2n[1+2*edge] == node){
	grid->edgeT[1+2*edge] = t;
	found = TRUE;
      }
    }
  }

  if (found) return grid;
  return NULL;
}

Grid *gridAddEdge(Grid *grid, int n0, int n1, 
		  int edgeId, double t0, double t1 )
{
  int edge;
  if ( grid->blanke2n == EMPTY ) return NULL;
  edge = grid->blanke2n;
  grid->blanke2n = grid->e2n[1+2*edge];
  grid->nedge++;

  grid->e2n[0+2*edge] = n0;
  grid->e2n[1+2*edge] = n1;
  grid->edgeId[edge]  = edgeId;
  grid->edgeT[0+2*edge] = t0;
  grid->edgeT[1+2*edge] = t1;

  if ( NULL == adjRegister( grid->edgeAdj, n0, edge ) ) return NULL;
  if ( NULL == adjRegister( grid->edgeAdj, n1, edge ) ) return NULL;

  return grid;
}

Grid *gridRemoveEdge(Grid *grid, int edge )
{
  if (edge >= grid->maxedge || edge < 0) return NULL;
  if (EMPTY == grid->e2n[2*edge]) return NULL;
  
  if ( grid->nedge <= 0) return NULL;
  grid->nedge--;

  if( ( NULL == adjRemove( grid->edgeAdj, grid->e2n[0+2*edge], edge ) ) || 
      ( NULL == adjRemove( grid->edgeAdj, grid->e2n[1+2*edge], edge ) ) )
    return NULL;  

  grid->e2n[0+2*edge] = EMPTY;
  grid->e2n[1+2*edge] = grid->blanke2n;
  grid->blanke2n = edge;

  return grid;
}

int gridFindEdge(Grid *grid, int n0, int n1 )
{
  AdjIterator it0, it1;
  Adj *adj=grid->edgeAdj;

  for ( it0 = adjFirst(adj,n0); adjValid(it0); it0 = adjNext(it0) )
    for ( it1 = adjFirst(adj,n1); adjValid(it1); it1 = adjNext(it1) )
      if ( adjItem(it0) == adjItem(it1) ) return adjItem(it0);

  return EMPTY;
}

int gridEdgeId(Grid *grid, int n0, int n1 )
{
  int edge = gridFindEdge(grid, n0, n1 );
  if ( edge == EMPTY ) return EMPTY;
  return grid->edgeId[edge];
}

int gridGeomCurveSize( Grid *grid, int edgeId, int startNode )
{
  AdjIterator it;
  int node, lastnode, edge, n1, nedgenode;
  bool found;

  nedgenode=0;
  node = startNode;
  lastnode = EMPTY;
  found = TRUE;
  while (found) {
    found = FALSE;
    for ( it = adjFirst(grid->edgeAdj,node); 
	  adjValid(it) && !found; 
	  it = adjNext(it)) {
      edge = adjItem(it);
      if (grid->edgeId[edge] == edgeId) {
	if ( node == grid->e2n[0+2*edge] ) {
	  n1 = grid->e2n[1+2*edge];
	}else{
	  n1 = grid->e2n[0+2*edge];	  
	}
	if ( n1 != lastnode ) { 
	  found = TRUE;
	  lastnode = node;
	  node = n1;
	  nedgenode++;
	}
      }
    }
  }
  if (nedgenode>0) nedgenode++;
  return nedgenode;
}

Grid *gridGeomCurve( Grid *grid, int edgeId, int startNode, int *curve )
{
  AdjIterator it;
  int node, lastnode, edge, n1, nedgenode;
  bool found;

  node = startNode;
  nedgenode=0;
  curve[nedgenode]=node;
  nedgenode++;
  lastnode = EMPTY;
  found = TRUE;
  while (found) {
    found = FALSE;
    for ( it = adjFirst(grid->edgeAdj,node); 
	  adjValid(it) && !found; 
	  it = adjNext(it)) {
      edge = adjItem(it);
      if (grid->edgeId[edge] == edgeId) {
	if ( node == grid->e2n[0+2*edge] ) {
	  n1 = grid->e2n[1+2*edge];
	}else{
	  n1 = grid->e2n[0+2*edge];	  
	}
	if ( n1 != lastnode ) { 
	  found = TRUE;
	  lastnode = node;
	  node = n1;
	  curve[nedgenode]=node;
	  nedgenode++;
	}
      }
    }
  }
  return grid;
}

Grid *gridGeomCurveT( Grid *grid, int edgeId, int startNode, double *curve )
{
  AdjIterator it;
  int node, lastnode, edge, n1, nedgenode;
  double t0, t1;
  bool found;

  node = startNode;
  nedgenode=1;
  lastnode = EMPTY;
  found = TRUE;
  while (found) {
    found = FALSE;
    for ( it = adjFirst(grid->edgeAdj,node); 
	  adjValid(it) && !found; 
	  it = adjNext(it)) {
      edge = adjItem(it);
      if (grid->edgeId[edge] == edgeId) {
	if ( node == grid->e2n[0+2*edge] ) {
	  n1 = grid->e2n[1+2*edge];
	  t0 = grid->edgeT[0+2*edge];
	  t1 = grid->edgeT[1+2*edge];
	}else{
	  n1 = grid->e2n[0+2*edge];	  
	  t0 = grid->edgeT[1+2*edge];
	  t1 = grid->edgeT[0+2*edge];
	}
	if ( n1 != lastnode ) { 
	  found = TRUE;
	  lastnode = node;
	  node = n1;
	  curve[nedgenode-1]=t0;
	  curve[nedgenode]=t1;
	  nedgenode++;
	}
      }
    }
  }
  return grid;
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

Grid *gridNodeXYZ(Grid *grid, int  node, double *xyz )
{
  if (node >=grid->nnode) return NULL;
  xyz[0] = grid->xyz[0+3*node];
  xyz[1] = grid->xyz[1+3*node];
  xyz[2] = grid->xyz[2+3*node];
  return grid;
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

Grid *gridARDerivative(Grid *grid, int node, double *ar, double *dARdx )
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

  double s1_dx, s2_dx, s3_dx, s4_dx, det_dx;
  double s1_dy, s2_dy, s3_dy, s4_dy, det_dy;
  double s1_dz, s2_dz, s3_dz, s4_dz, det_dz;
  double xr_dx, yr_dx, zr_dx;
  double xr_dy, yr_dy, zr_dy;
  double xr_dz, yr_dz, zr_dz;
  double circ_dx, circ_dy, circ_dz; 

  double nx1_dx, ny1_dx, nz1_dx, rmag1_dx;
  double nx1_dy, ny1_dy, nz1_dy, rmag1_dy;
  double nx1_dz, ny1_dz, nz1_dz, rmag1_dz;

  double nx2_dx, ny2_dx, nz2_dx, rmag2_dx;
  double nx2_dy, ny2_dy, nz2_dy, rmag2_dy;
  double nx2_dz, ny2_dz, nz2_dz, rmag2_dz;

  double nx3_dx, ny3_dx, nz3_dx, rmag3_dx;
  double nx3_dy, ny3_dy, nz3_dy, rmag3_dy;
  double nx3_dz, ny3_dz, nz3_dz, rmag3_dz;

  double nx4_dx, ny4_dx, nz4_dx;
  double nx4_dy, ny4_dy, nz4_dy;
  double nx4_dz, ny4_dz, nz4_dz;

  double xins_dx, xins_dy, xins_dz;

  x1 = grid->xyz[0+3*grid->c2n[0+4*node]];
  y1 = grid->xyz[1+3*grid->c2n[0+4*node]];
  z1 = grid->xyz[2+3*grid->c2n[0+4*node]];

  x2 = grid->xyz[0+3*grid->c2n[1+4*node]];
  y2 = grid->xyz[1+3*grid->c2n[1+4*node]];
  z2 = grid->xyz[2+3*grid->c2n[1+4*node]];

  x3 = grid->xyz[0+3*grid->c2n[2+4*node]];
  y3 = grid->xyz[1+3*grid->c2n[2+4*node]];
  z3 = grid->xyz[2+3*grid->c2n[2+4*node]];

  x4 = grid->xyz[0+3*grid->c2n[3+4*node]];
  y4 = grid->xyz[1+3*grid->c2n[3+4*node]];
  z4 = grid->xyz[2+3*grid->c2n[3+4*node]];

  /* Compute the aspect ratios */

        det = (x4-x1)*((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
            + (y4-y1)*((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
            + (z4-z1)*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

	det_dx = -((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
   	       + (y4-y1)*(-(z2-z1) + (z3-z1))
	       + (z4-z1)*(-(y3-y1) + (y2-y1));

        det_dy = (x4-x1)*(-(z3-z1) +(z2-z1))
            - ((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
            + (z4-z1)*(-(x2-x1) + (x3-x1));

        det_dz = (x4-x1)*(-(y2-y1) + (y3-y1))
            + (y4-y1)*(-(x3-x1) + (x2-x1))
            -((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

         s1 = ((x2*x2 + y2*y2 + z2*z2) - (x1*x1 + y1*y1 + z1*z1))/2.0;
         s2 = ((x3*x3 + y3*y3 + z3*z3) - (x1*x1 + y1*y1 + z1*z1))/2.0;
         s3 = ((x4*x4 + y4*y4 + z4*z4) - (x1*x1 + y1*y1 + z1*z1))/2.0;

	 s1_dx = -x1;
	 s1_dy = -y1;
	 s1_dz = -z1;
	 s2_dx = -x1;
	 s2_dy = -y1;
	 s2_dz = -z1;
	 s3_dx = -x1;
	 s3_dy = -y1;
	 s3_dz = -z1;


         xr  = (   s3    * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
               + (y4-y1) * (    s2  *(z2-z1) -   s1   *(z3-z1))
               + (z4-z1) * (    s1  *(y3-y1) -   s2   *(y2-y1)) ) / det;


	 xr_dx = (  s3_dx  * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
                 + (y4-y1) * (  s2_dx *(z2-z1) -  s1_dx *(z3-z1))
                 + (z4-z1) * (  s1_dx *(y3-y1) -  s2_dx *(y2-y1)) ) ;

	 xr_dx = det * xr_dx - (xr*det) * det_dx ;
	 xr_dx = xr_dx / det / det;

	 xr_dy = (   s3    * (        -(z3-z1) +         (z2-z1))
                 +  s3_dy  * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
                 + (y4-y1) * (  s2_dy *(z2-z1) -  s1_dy *(z3-z1))
                 -           (    s2  *(z2-z1) -   s1   *(z3-z1))
                 + (z4-z1) * ( s1_dy  *(y3-y1) -  s2_dy *(y2-y1)) 
                 + (z4-z1) * (        -s1      +        s2     ) ) ;

	 xr_dy = det * xr_dy - (xr*det) * det_dy ;
	 xr_dy = xr_dy / det / det;


         xr_dz = (   s3    * (           -(y2-y1) +            (y3-y1))
	         +  s3_dz  * (    (y2-y1)*(z3-z1) -    (y3-y1)*(z2-z1))
                 + (y4-y1) * (    -s2             +    s1             )
                 + (y4-y1) * (    s2_dz  *(z2-z1) -   s1_dz   *(z3-z1))
                 + (z4-z1) * (    s1_dz  *(y3-y1) -   s2_dz   *(y2-y1)) 
                 -           (    s1     *(y3-y1) -   s2      *(y2-y1)) );
	 xr_dz = det * xr_dz - (xr*det) * det_dz ;
	 xr_dz = xr_dz / det / det;

         yr  =((x4-x1)*(s1     *(z3-z1) - s2     *(z2-z1))
             + s3     *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
             + (z4-z1)*((x2-x1)*s2      - (x3-x1)*s1     ))/det;

         yr_dx = ( (x4-x1)*(s1_dx  *(z3-z1) - s2_dx  *(z2-z1))
		 -         (s1     *(z3-z1) - s2     *(z2-z1))
                 + s3     *(       -(z2-z1)          +(z3-z1))
                 + s3_dx  *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                 + (z4-z1)*(       -s2               +s1     )
                 + (z4-z1)*((x2-x1)*s2_dx   - (x3-x1)*s1_dx  ));

	 yr_dx = det * yr_dx - (yr*det) * det_dx ;
	 yr_dx = yr_dx / det / det;

         yr_dy  = ( (x4-x1)*(s1_dy  *(z3-z1) - s2_dy  *(z2-z1))
                  +  s3_dy *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                  + (z4-z1)*((x2-x1)*s2_dy   - (x3-x1)*s1_dy  ) );

	 yr_dy = det * yr_dy - (yr*det) * det_dy ;
	 yr_dy = yr_dy / det / det;


         yr_dz = ( (x4-x1)*(  -s1           + s2             )
		 + (x4-x1)*( s1_dz *(z3-z1) - s2_dz  *(z2-z1))
                 + s3     *(-(x3-x1)        + (x2-x1)        )
                 + s3_dz  *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                 + (z4-z1)*((x2-x1)*s2_dz   - (x3-x1)*s1_dz  )
                 -         ((x2-x1)*s2      - (x3-x1)*s1     ) );

	 yr_dz = det * yr_dz - (yr*det) * det_dz ;
	 yr_dz = yr_dz / det / det;

         zr  =((x4-x1)*((y2-y1)*s2      - (y3-y1)*s1     )
             + (y4-y1)*((x3-x1)*s1      - (x2-x1)*s2     )
             + s3     *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)))/det;

         zr_dx = ( (x4-x1)*((y2-y1)*s2_dx   - (y3-y1)*s1_dx  )
		 -         ((y2-y1)*s2      - (y3-y1)*s1     )
                 + (y4-y1)*(       -s1      +         s2     )
                 + (y4-y1)*((x3-x1)*s1_dx   - (x2-x1)*s2_dx  )
                 + s3     *(       -(y3-y1) +         (y2-y1))
                 + s3_dx  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dx = det * zr_dx - (zr*det) * det_dx ;
	 zr_dx = zr_dx / det / det;

         zr_dy = ( (x4-x1)*((y2-y1)*s2_dy   - (y3-y1)*s1_dy  )
		 + (x4-x1)*(       -s2      +         s1     )
                 + (y4-y1)*((x3-x1)*s1_dy   - (x2-x1)*s2_dy  )
                 -         ((x3-x1)*s1      - (x2-x1)*s2     )
                 + s3     *(-(x2-x1)        + (x3-x1)        )
                 + s3_dy  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dy = det * zr_dy - (zr*det) * det_dy ;
	 zr_dy = zr_dy / det / det;

         zr_dz = ( (x4-x1)*((y2-y1)*s2_dz   - (y3-y1)*s1_dz  )
                 + (y4-y1)*((x3-x1)*s1_dz   - (x2-x1)*s2_dz  )
                 + s3_dz  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dz = det * zr_dz - (zr*det) * det_dz ;
	 zr_dz = zr_dz / det / det;

         circ = sqrt((x1-xr)*(x1-xr) + (y1-yr)*(y1-yr) + (z1-zr)*(z1-zr));

	 circ_dx = 2.0*(x1-xr)*(1.0-xr_dx) 
	         + 2.0*(y1-yr)*(   -yr_dx) 
	         + 2.0*(z1-zr)*(   -zr_dx);

	 circ_dx = 0.5 / circ * circ_dx;

	 circ_dy = 2.0*(x1-xr)*(   -xr_dy) 
	         + 2.0*(y1-yr)*(1.0-yr_dy) 
	         + 2.0*(z1-zr)*(   -zr_dy);

	 circ_dy = 0.5 / circ * circ_dy;

	 circ_dz = 2.0*(x1-xr)*(   -xr_dz) 
	         + 2.0*(y1-yr)*(   -yr_dz) 
	         + 2.0*(z1-zr)*(1.0-zr_dz);

	 circ_dz = 0.5 / circ * circ_dz;
 
  /* Get the in-circle */

         nx1 = (y3 - y1)*(z2 - z1) - (y2 - y1)*(z3 - z1);
         ny1 =-(x3 - x1)*(z2 - z1) + (x2 - x1)*(z3 - z1);
         nz1 = (x3 - x1)*(y2 - y1) - (x2 - x1)*(y3 - y1);

	 nx1_dx = 0.0;
	 ny1_dx = (z2 - z1) - (z3 - z1);
         nz1_dx =-(y2 - y1) + (y3 - y1);

	 nx1_dy = -(z2 - z1) + (z3 - z1);
	 ny1_dy = 0.0;
         nz1_dy = -(x3 - x1) + (x2 - x1);
 
	 nx1_dz = -(y3 - y1) + (y2 - y1);
	 ny1_dz =  (x3 - x1) - (x2 - x1);
         nz1_dz = 0.0;

         nx2 = (y2 - y1)*(z4 - z1) - (y4 - y1)*(z2 - z1);
         ny2 =-(x2 - x1)*(z4 - z1) + (x4 - x1)*(z2 - z1);
         nz2 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);

         nx2_dx =  0.0;
         ny2_dx =  (z4 - z1) - (z2 - z1);
         nz2_dx = -(y4 - y1) + (y2 - y1);

         nx2_dy = -(z4 - z1) + (z2 - z1);
         ny2_dy =  0.0;
         nz2_dy = -(x2 - x1) + (x4 - x1);

         nx2_dz = -(y2 - y1) + (y4 - y1);
         ny2_dz =  (x2 - x1) - (x4 - x1);
         nz2_dz =  0.0;

         nx3 = (y4 - y1)*(z3 - z1) - (y3 - y1)*(z4 - z1);
         ny3 =-(x4 - x1)*(z3 - z1) + (x3 - x1)*(z4 - z1);
         nz3 = (x4 - x1)*(y3 - y1) - (x3 - x1)*(y4 - y1);

         nx3_dx =  0.0;
         ny3_dx =  (z3 - z1) - (z4 - z1);
         nz3_dx = -(y3 - y1) + (y4 - y1);

         nx3_dy = -(z3 - z1) + (z4 - z1);
         ny3_dy = 0.0;
         nz3_dy = -(x4 - x1) + (x3 - x1);

         nx3_dz =-(y4 - y1) + (y3 - y1);
         ny3_dz = (x4 - x1) - (x3 - x1);
         nz3_dz = 0.0;

         nx4 = (y3 - y2)*(z4 - z2) - (y4 - y2)*(z3 - z2);
         ny4 =-(x3 - x2)*(z4 - z2) + (x4 - x2)*(z3 - z2);
         nz4 = (x3 - x2)*(y4 - y2) - (x4 - x2)*(y3 - y2);

         rmag1 = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);

	 rmag1_dx = (nx1*nx1_dx + ny1*ny1_dx + nz1*nz1_dx) / rmag1;
	 rmag1_dy = (nx1*nx1_dy + ny1*ny1_dy + nz1*nz1_dy) / rmag1;
	 rmag1_dz = (nx1*nx1_dz + ny1*ny1_dz + nz1*nz1_dz) / rmag1;

         rmag2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);

	 rmag2_dx = (nx2*nx2_dx + ny2*ny2_dx + nz2*nz2_dx) / rmag2;
	 rmag2_dy = (nx2*nx2_dy + ny2*ny2_dy + nz2*nz2_dy) / rmag2;
	 rmag2_dz = (nx2*nx2_dz + ny2*ny2_dz + nz2*nz2_dz) / rmag2;

         rmag3 = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3);

	 rmag3_dx = (nx3*nx3_dx + ny3*ny3_dx + nz3*nz3_dx) / rmag3;
	 rmag3_dy = (nx3*nx3_dy + ny3*ny3_dy + nz3*nz3_dy) / rmag3;
	 rmag3_dz = (nx3*nx3_dz + ny3*ny3_dz + nz3*nz3_dz) / rmag3;

         rmag4 = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4);
	
	 /* note that the derivatives procede the normalization */
	 /* because I need the un-nornmalized value for deriavtives */
 
         nx1_dx = ( rmag1 * nx1_dx - nx1 * rmag1_dx ) / rmag1 / rmag1;
         ny1_dx = ( rmag1 * ny1_dx - ny1 * rmag1_dx ) / rmag1 / rmag1;
         nz1_dx = ( rmag1 * nz1_dx - nz1 * rmag1_dx ) / rmag1 / rmag1;

         nx1_dy = ( rmag1 * nx1_dy - nx1 * rmag1_dy ) / rmag1 / rmag1;
         ny1_dy = ( rmag1 * ny1_dy - ny1 * rmag1_dy ) / rmag1 / rmag1;
         nz1_dy = ( rmag1 * nz1_dy - nz1 * rmag1_dy ) / rmag1 / rmag1;

         nx1_dz = ( rmag1 * nx1_dz - nx1 * rmag1_dz ) / rmag1 / rmag1;
         ny1_dz = ( rmag1 * ny1_dz - ny1 * rmag1_dz ) / rmag1 / rmag1;
         nz1_dz = ( rmag1 * nz1_dz - nz1 * rmag1_dz ) / rmag1 / rmag1;

         nx1 = nx1/rmag1;
         ny1 = ny1/rmag1;
         nz1 = nz1/rmag1;

         nx2_dx = ( rmag2 * nx2_dx - nx2 * rmag2_dx ) / rmag2 / rmag2;
         ny2_dx = ( rmag2 * ny2_dx - ny2 * rmag2_dx ) / rmag2 / rmag2;
         nz2_dx = ( rmag2 * nz2_dx - nz2 * rmag2_dx ) / rmag2 / rmag2;

         nx2_dy = ( rmag2 * nx2_dy - nx2 * rmag2_dy ) / rmag2 / rmag2;
         ny2_dy = ( rmag2 * ny2_dy - ny2 * rmag2_dy ) / rmag2 / rmag2;
         nz2_dy = ( rmag2 * nz2_dy - nz2 * rmag2_dy ) / rmag2 / rmag2;

         nx2_dz = ( rmag2 * nx2_dz - nx2 * rmag2_dz ) / rmag2 / rmag2;
         ny2_dz = ( rmag2 * ny2_dz - ny2 * rmag2_dz ) / rmag2 / rmag2;
         nz2_dz = ( rmag2 * nz2_dz - nz2 * rmag2_dz ) / rmag2 / rmag2;

         nx2 = nx2/rmag2;
         ny2 = ny2/rmag2;
         nz2 = nz2/rmag2;

         nx3_dx = ( rmag3 * nx3_dx - nx3 * rmag3_dx ) / rmag3 / rmag3;
         ny3_dx = ( rmag3 * ny3_dx - ny3 * rmag3_dx ) / rmag3 / rmag3;
         nz3_dx = ( rmag3 * nz3_dx - nz3 * rmag3_dx ) / rmag3 / rmag3;

         nx3_dy = ( rmag3 * nx3_dy - nx3 * rmag3_dy ) / rmag3 / rmag3;
         ny3_dy = ( rmag3 * ny3_dy - ny3 * rmag3_dy ) / rmag3 / rmag3;
         nz3_dy = ( rmag3 * nz3_dy - nz3 * rmag3_dy ) / rmag3 / rmag3;

         nx3_dz = ( rmag3 * nx3_dz - nx3 * rmag3_dz ) / rmag3 / rmag3;
         ny3_dz = ( rmag3 * ny3_dz - ny3 * rmag3_dz ) / rmag3 / rmag3;
         nz3_dz = ( rmag3 * nz3_dz - nz3 * rmag3_dz ) / rmag3 / rmag3;

         nx3 = nx3/rmag3;
         ny3 = ny3/rmag3;
         nz3 = nz3/rmag3;

         nx4_dx = 0.0;
         ny4_dx = 0.0;
         nz4_dx = 0.0;

         nx4_dy = 0.0;
         ny4_dy = 0.0;
         nz4_dy = 0.0;

         nx4_dz = 0.0;
         ny4_dz = 0.0;
         nz4_dz = 0.0;

         nx4 = nx4/rmag4;
         ny4 = ny4/rmag4;
         nz4 = nz4/rmag4;

         det= -nx3*ny2*nz1 + nx4*ny2*nz1 + nx2*ny3*nz1 - nx4*ny3*nz1
              -nx2*ny4*nz1 + nx3*ny4*nz1 + nx3*ny1*nz2 - nx4*ny1*nz2
              -nx1*ny3*nz2 + nx4*ny3*nz2 + nx1*ny4*nz2 - nx3*ny4*nz2
              -nx2*ny1*nz3 + nx4*ny1*nz3 + nx1*ny2*nz3 - nx4*ny2*nz3
              -nx1*ny4*nz3 + nx2*ny4*nz3 + nx2*ny1*nz4 - nx3*ny1*nz4
              -nx1*ny2*nz4 + nx3*ny2*nz4 + nx1*ny3*nz4 - nx2*ny3*nz4;

         det_dx = 
	   - nx3_dx*ny2*nz1 - nx3*ny2_dx*nz1 - nx3*ny2*nz1_dx 
	   + nx4_dx*ny2*nz1 + nx4*ny2_dx*nz1 + nx4*ny2*nz1_dx 
	   + nx2_dx*ny3*nz1 + nx2*ny3_dx*nz1 + nx2*ny3*nz1_dx 
	   - nx4_dx*ny3*nz1 - nx4*ny3_dx*nz1 - nx4*ny3*nz1_dx
	   - nx2_dx*ny4*nz1 - nx2*ny4_dx*nz1 - nx2*ny4*nz1_dx 
	   + nx3_dx*ny4*nz1 + nx3*ny4_dx*nz1 + nx3*ny4*nz1_dx
	   + nx3_dx*ny1*nz2 + nx3*ny1_dx*nz2 + nx3*ny1*nz2_dx 
	   - nx4_dx*ny1*nz2 - nx4*ny1_dx*nz2 - nx4*ny1*nz2_dx
	   - nx1_dx*ny3*nz2 - nx1*ny3_dx*nz2 - nx1*ny3*nz2_dx 
	   + nx4_dx*ny3*nz2 + nx4*ny3_dx*nz2 + nx4*ny3*nz2_dx 
	   + nx1_dx*ny4*nz2 + nx1*ny4_dx*nz2 + nx1*ny4*nz2_dx 
	   - nx3_dx*ny4*nz2 - nx3*ny4_dx*nz2 - nx3*ny4*nz2_dx
	   - nx2_dx*ny1*nz3 - nx2*ny1_dx*nz3 - nx2*ny1*nz3_dx 
	   + nx4_dx*ny1*nz3 + nx4*ny1_dx*nz3 + nx4*ny1*nz3_dx 
	   + nx1_dx*ny2*nz3 + nx1*ny2_dx*nz3 + nx1*ny2*nz3_dx 
	   - nx4_dx*ny2*nz3 - nx4*ny2_dx*nz3 - nx4*ny2*nz3_dx
	   - nx1_dx*ny4*nz3 - nx1*ny4_dx*nz3 - nx1*ny4*nz3_dx 
	   + nx2_dx*ny4*nz3 + nx2*ny4_dx*nz3 + nx2*ny4*nz3_dx 
	   + nx2_dx*ny1*nz4 + nx2*ny1_dx*nz4 + nx2*ny1*nz4_dx 
	   - nx3_dx*ny1*nz4 - nx3*ny1_dx*nz4 - nx3*ny1*nz4_dx
	   - nx1_dx*ny2*nz4 - nx1*ny2_dx*nz4 - nx1*ny2*nz4_dx 
	   + nx3_dx*ny2*nz4 + nx3*ny2_dx*nz4 + nx3*ny2*nz4_dx 
	   + nx1_dx*ny3*nz4 + nx1*ny3_dx*nz4 + nx1*ny3*nz4_dx 
	   - nx2_dx*ny3*nz4 - nx2*ny3_dx*nz4 - nx2*ny3*nz4_dx;

         det_dy = 
	   - nx3_dy*ny2*nz1 - nx3*ny2_dy*nz1 - nx3*ny2*nz1_dy 
	   + nx4_dy*ny2*nz1 + nx4*ny2_dy*nz1 + nx4*ny2*nz1_dy 
	   + nx2_dy*ny3*nz1 + nx2*ny3_dy*nz1 + nx2*ny3*nz1_dy 
	   - nx4_dy*ny3*nz1 - nx4*ny3_dy*nz1 - nx4*ny3*nz1_dy
	   - nx2_dy*ny4*nz1 - nx2*ny4_dy*nz1 - nx2*ny4*nz1_dy 
	   + nx3_dy*ny4*nz1 + nx3*ny4_dy*nz1 + nx3*ny4*nz1_dy
	   + nx3_dy*ny1*nz2 + nx3*ny1_dy*nz2 + nx3*ny1*nz2_dy 
	   - nx4_dy*ny1*nz2 - nx4*ny1_dy*nz2 - nx4*ny1*nz2_dy
	   - nx1_dy*ny3*nz2 - nx1*ny3_dy*nz2 - nx1*ny3*nz2_dy 
	   + nx4_dy*ny3*nz2 + nx4*ny3_dy*nz2 + nx4*ny3*nz2_dy 
	   + nx1_dy*ny4*nz2 + nx1*ny4_dy*nz2 + nx1*ny4*nz2_dy 
	   - nx3_dy*ny4*nz2 - nx3*ny4_dy*nz2 - nx3*ny4*nz2_dy
	   - nx2_dy*ny1*nz3 - nx2*ny1_dy*nz3 - nx2*ny1*nz3_dy 
	   + nx4_dy*ny1*nz3 + nx4*ny1_dy*nz3 + nx4*ny1*nz3_dy 
	   + nx1_dy*ny2*nz3 + nx1*ny2_dy*nz3 + nx1*ny2*nz3_dy 
	   - nx4_dy*ny2*nz3 - nx4*ny2_dy*nz3 - nx4*ny2*nz3_dy
	   - nx1_dy*ny4*nz3 - nx1*ny4_dy*nz3 - nx1*ny4*nz3_dy 
	   + nx2_dy*ny4*nz3 + nx2*ny4_dy*nz3 + nx2*ny4*nz3_dy 
	   + nx2_dy*ny1*nz4 + nx2*ny1_dy*nz4 + nx2*ny1*nz4_dy 
	   - nx3_dy*ny1*nz4 - nx3*ny1_dy*nz4 - nx3*ny1*nz4_dy
	   - nx1_dy*ny2*nz4 - nx1*ny2_dy*nz4 - nx1*ny2*nz4_dy 
	   + nx3_dy*ny2*nz4 + nx3*ny2_dy*nz4 + nx3*ny2*nz4_dy 
	   + nx1_dy*ny3*nz4 + nx1*ny3_dy*nz4 + nx1*ny3*nz4_dy 
	   - nx2_dy*ny3*nz4 - nx2*ny3_dy*nz4 - nx2*ny3*nz4_dy;

         det_dz = 
	   - nx3_dz*ny2*nz1 - nx3*ny2_dz*nz1 - nx3*ny2*nz1_dz 
	   + nx4_dz*ny2*nz1 + nx4*ny2_dz*nz1 + nx4*ny2*nz1_dz 
	   + nx2_dz*ny3*nz1 + nx2*ny3_dz*nz1 + nx2*ny3*nz1_dz 
	   - nx4_dz*ny3*nz1 - nx4*ny3_dz*nz1 - nx4*ny3*nz1_dz
	   - nx2_dz*ny4*nz1 - nx2*ny4_dz*nz1 - nx2*ny4*nz1_dz 
	   + nx3_dz*ny4*nz1 + nx3*ny4_dz*nz1 + nx3*ny4*nz1_dz
	   + nx3_dz*ny1*nz2 + nx3*ny1_dz*nz2 + nx3*ny1*nz2_dz 
	   - nx4_dz*ny1*nz2 - nx4*ny1_dz*nz2 - nx4*ny1*nz2_dz
	   - nx1_dz*ny3*nz2 - nx1*ny3_dz*nz2 - nx1*ny3*nz2_dz 
	   + nx4_dz*ny3*nz2 + nx4*ny3_dz*nz2 + nx4*ny3*nz2_dz 
	   + nx1_dz*ny4*nz2 + nx1*ny4_dz*nz2 + nx1*ny4*nz2_dz 
	   - nx3_dz*ny4*nz2 - nx3*ny4_dz*nz2 - nx3*ny4*nz2_dz
	   - nx2_dz*ny1*nz3 - nx2*ny1_dz*nz3 - nx2*ny1*nz3_dz 
	   + nx4_dz*ny1*nz3 + nx4*ny1_dz*nz3 + nx4*ny1*nz3_dz 
	   + nx1_dz*ny2*nz3 + nx1*ny2_dz*nz3 + nx1*ny2*nz3_dz 
	   - nx4_dz*ny2*nz3 - nx4*ny2_dz*nz3 - nx4*ny2*nz3_dz
	   - nx1_dz*ny4*nz3 - nx1*ny4_dz*nz3 - nx1*ny4*nz3_dz 
	   + nx2_dz*ny4*nz3 + nx2*ny4_dz*nz3 + nx2*ny4*nz3_dz 
	   + nx2_dz*ny1*nz4 + nx2*ny1_dz*nz4 + nx2*ny1*nz4_dz 
	   - nx3_dz*ny1*nz4 - nx3*ny1_dz*nz4 - nx3*ny1*nz4_dz
	   - nx1_dz*ny2*nz4 - nx1*ny2_dz*nz4 - nx1*ny2*nz4_dz 
	   + nx3_dz*ny2*nz4 + nx3*ny2_dz*nz4 + nx3*ny2*nz4_dz 
	   + nx1_dz*ny3*nz4 + nx1*ny3_dz*nz4 + nx1*ny3*nz4_dz 
	   - nx2_dz*ny3*nz4 - nx2*ny3_dz*nz4 - nx2*ny3*nz4_dz;

         s1 = nx1*x1 + ny1*y1 + nz1*z1;

         s1_dx = nx1_dx*x1 + nx1 + ny1_dx*y1       + nz1_dx*z1;
         s1_dy = nx1_dy*x1       + ny1_dy*y1 + ny1 + nz1_dy*z1;
         s1_dz = nx1_dz*x1       + ny1_dz*y1       + nz1_dz*z1 + nz1;

         s2 = nx2*x1 + ny2*y1 + nz2*z1;

         s2_dx = nx2_dx*x1 + nx2 + ny2_dx*y1       + nz2_dx*z1;
         s2_dy = nx2_dy*x1       + ny2_dy*y1 + ny2 + nz2_dy*z1;
         s2_dz = nx2_dz*x1       + ny2_dz*y1       + nz2_dz*z1 + nz2;

         s3 = nx3*x1 + ny3*y1 + nz3*z1;

         s3_dx = nx3_dx*x1 + nx3 + ny3_dx*y1       + nz3_dx*z1;
         s3_dy = nx3_dy*x1       + ny3_dy*y1 + ny3 + nz3_dy*z1;
         s3_dz = nx3_dz*x1       + ny3_dz*y1       + nz3_dz*z1 + nz3;

         s4 = nx4*x4 + ny4*y4 + nz4*z4;

         s4_dx = 0.0;
         s4_dy = 0.0;
         s4_dz = 0.0;

         xins = nx4*ny3*nz2*s1 - nx3*ny4*nz2*s1 - nx4*ny2*nz3*s1 +
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
                nx1*ny2*nz3*s4;

         xins_dx = 
  nx4_dx*ny3*nz2*s1 + nx4*ny3_dx*nz2*s1 + nx4*ny3*nz2_dx*s1 + nx4*ny3*nz2*s1_dx 
- nx3_dx*ny4*nz2*s1 - nx3*ny4_dx*nz2*s1 - nx3*ny4*nz2_dx*s1 - nx3*ny4*nz2*s1_dx 
- nx4_dx*ny2*nz3*s1 - nx4*ny2_dx*nz3*s1 - nx4*ny2*nz3_dx*s1 - nx4*ny2*nz3*s1_dx 
+ nx2_dx*ny4*nz3*s1 + nx2*ny4_dx*nz3*s1 + nx2*ny4*nz3_dx*s1 + nx2*ny4*nz3*s1_dx 
+ nx3_dx*ny2*nz4*s1 + nx3*ny2_dx*nz4*s1 + nx3*ny2*nz4_dx*s1 + nx3*ny2*nz4*s1_dx 
- nx2_dx*ny3*nz4*s1 - nx2*ny3_dx*nz4*s1 - nx2*ny3*nz4_dx*s1 - nx2*ny3*nz4*s1_dx 
- nx4_dx*ny3*nz1*s2 - nx4*ny3_dx*nz1*s2 - nx4*ny3*nz1_dx*s2 - nx4*ny3*nz1*s2_dx 
+ nx3_dx*ny4*nz1*s2 + nx3*ny4_dx*nz1*s2 + nx3*ny4*nz1_dx*s2 + nx3*ny4*nz1*s2_dx 
+ nx4_dx*ny1*nz3*s2 + nx4*ny1_dx*nz3*s2 + nx4*ny1*nz3_dx*s2 + nx4*ny1*nz3*s2_dx 
- nx1_dx*ny4*nz3*s2 - nx1*ny4_dx*nz3*s2 - nx1*ny4*nz3_dx*s2 - nx1*ny4*nz3*s2_dx 
- nx3_dx*ny1*nz4*s2 - nx3*ny1_dx*nz4*s2 - nx3*ny1*nz4_dx*s2 - nx3*ny1*nz4*s2_dx 
+ nx1_dx*ny3*nz4*s2 + nx1*ny3_dx*nz4*s2 + nx1*ny3*nz4_dx*s2 + nx1*ny3*nz4*s2_dx 
+ nx4_dx*ny2*nz1*s3 + nx4*ny2_dx*nz1*s3 + nx4*ny2*nz1_dx*s3 + nx4*ny2*nz1*s3_dx 
- nx2_dx*ny4*nz1*s3 - nx2*ny4_dx*nz1*s3 - nx2*ny4*nz1_dx*s3 - nx2*ny4*nz1*s3_dx 
- nx4_dx*ny1*nz2*s3 - nx4*ny1_dx*nz2*s3 - nx4*ny1*nz2_dx*s3 - nx4*ny1*nz2*s3_dx 
+ nx1_dx*ny4*nz2*s3 + nx1*ny4_dx*nz2*s3 + nx1*ny4*nz2_dx*s3 + nx1*ny4*nz2*s3_dx 
+ nx2_dx*ny1*nz4*s3 + nx2*ny1_dx*nz4*s3 + nx2*ny1*nz4_dx*s3 + nx2*ny1*nz4*s3_dx 
- nx1_dx*ny2*nz4*s3 - nx1*ny2_dx*nz4*s3 - nx1*ny2*nz4_dx*s3 - nx1*ny2*nz4*s3_dx 
- nx3_dx*ny2*nz1*s4 - nx3*ny2_dx*nz1*s4 - nx3*ny2*nz1_dx*s4 - nx3*ny2*nz1*s4_dx 
+ nx2_dx*ny3*nz1*s4 + nx2*ny3_dx*nz1*s4 + nx2*ny3*nz1_dx*s4 + nx2*ny3*nz1*s4_dx 
+ nx3_dx*ny1*nz2*s4 + nx3*ny1_dx*nz2*s4 + nx3*ny1*nz2_dx*s4 + nx3*ny1*nz2*s4_dx 
- nx1_dx*ny3*nz2*s4 - nx1*ny3_dx*nz2*s4 - nx1*ny3*nz2_dx*s4 - nx1*ny3*nz2*s4_dx 
- nx2_dx*ny1*nz3*s4 - nx2*ny1_dx*nz3*s4 - nx2*ny1*nz3_dx*s4 - nx2*ny1*nz3*s4_dx 
+ nx1_dx*ny2*nz3*s4 + nx1*ny2_dx*nz3*s4 + nx1*ny2*nz3_dx*s4 + nx1*ny2*nz3*s4_dx;

         xins_dy = 
  nx4_dy*ny3*nz2*s1 + nx4*ny3_dy*nz2*s1 + nx4*ny3*nz2_dy*s1 + nx4*ny3*nz2*s1_dy 
- nx3_dy*ny4*nz2*s1 - nx3*ny4_dy*nz2*s1 - nx3*ny4*nz2_dy*s1 - nx3*ny4*nz2*s1_dy 
- nx4_dy*ny2*nz3*s1 - nx4*ny2_dy*nz3*s1 - nx4*ny2*nz3_dy*s1 - nx4*ny2*nz3*s1_dy 
+ nx2_dy*ny4*nz3*s1 + nx2*ny4_dy*nz3*s1 + nx2*ny4*nz3_dy*s1 + nx2*ny4*nz3*s1_dy 
+ nx3_dy*ny2*nz4*s1 + nx3*ny2_dy*nz4*s1 + nx3*ny2*nz4_dy*s1 + nx3*ny2*nz4*s1_dy 
- nx2_dy*ny3*nz4*s1 - nx2*ny3_dy*nz4*s1 - nx2*ny3*nz4_dy*s1 - nx2*ny3*nz4*s1_dy 
- nx4_dy*ny3*nz1*s2 - nx4*ny3_dy*nz1*s2 - nx4*ny3*nz1_dy*s2 - nx4*ny3*nz1*s2_dy 
+ nx3_dy*ny4*nz1*s2 + nx3*ny4_dy*nz1*s2 + nx3*ny4*nz1_dy*s2 + nx3*ny4*nz1*s2_dy 
+ nx4_dy*ny1*nz3*s2 + nx4*ny1_dy*nz3*s2 + nx4*ny1*nz3_dy*s2 + nx4*ny1*nz3*s2_dy 
- nx1_dy*ny4*nz3*s2 - nx1*ny4_dy*nz3*s2 - nx1*ny4*nz3_dy*s2 - nx1*ny4*nz3*s2_dy 
- nx3_dy*ny1*nz4*s2 - nx3*ny1_dy*nz4*s2 - nx3*ny1*nz4_dy*s2 - nx3*ny1*nz4*s2_dy 
+ nx1_dy*ny3*nz4*s2 + nx1*ny3_dy*nz4*s2 + nx1*ny3*nz4_dy*s2 + nx1*ny3*nz4*s2_dy 
+ nx4_dy*ny2*nz1*s3 + nx4*ny2_dy*nz1*s3 + nx4*ny2*nz1_dy*s3 + nx4*ny2*nz1*s3_dy 
- nx2_dy*ny4*nz1*s3 - nx2*ny4_dy*nz1*s3 - nx2*ny4*nz1_dy*s3 - nx2*ny4*nz1*s3_dy 
- nx4_dy*ny1*nz2*s3 - nx4*ny1_dy*nz2*s3 - nx4*ny1*nz2_dy*s3 - nx4*ny1*nz2*s3_dy 
+ nx1_dy*ny4*nz2*s3 + nx1*ny4_dy*nz2*s3 + nx1*ny4*nz2_dy*s3 + nx1*ny4*nz2*s3_dy 
+ nx2_dy*ny1*nz4*s3 + nx2*ny1_dy*nz4*s3 + nx2*ny1*nz4_dy*s3 + nx2*ny1*nz4*s3_dy 
- nx1_dy*ny2*nz4*s3 - nx1*ny2_dy*nz4*s3 - nx1*ny2*nz4_dy*s3 - nx1*ny2*nz4*s3_dy 
- nx3_dy*ny2*nz1*s4 - nx3*ny2_dy*nz1*s4 - nx3*ny2*nz1_dy*s4 - nx3*ny2*nz1*s4_dy 
+ nx2_dy*ny3*nz1*s4 + nx2*ny3_dy*nz1*s4 + nx2*ny3*nz1_dy*s4 + nx2*ny3*nz1*s4_dy 
+ nx3_dy*ny1*nz2*s4 + nx3*ny1_dy*nz2*s4 + nx3*ny1*nz2_dy*s4 + nx3*ny1*nz2*s4_dy 
- nx1_dy*ny3*nz2*s4 - nx1*ny3_dy*nz2*s4 - nx1*ny3*nz2_dy*s4 - nx1*ny3*nz2*s4_dy 
- nx2_dy*ny1*nz3*s4 - nx2*ny1_dy*nz3*s4 - nx2*ny1*nz3_dy*s4 - nx2*ny1*nz3*s4_dy 
+ nx1_dy*ny2*nz3*s4 + nx1*ny2_dy*nz3*s4 + nx1*ny2*nz3_dy*s4 + nx1*ny2*nz3*s4_dy;

         xins_dz = 
  nx4_dz*ny3*nz2*s1 + nx4*ny3_dz*nz2*s1 + nx4*ny3*nz2_dz*s1 + nx4*ny3*nz2*s1_dz 
- nx3_dz*ny4*nz2*s1 - nx3*ny4_dz*nz2*s1 - nx3*ny4*nz2_dz*s1 - nx3*ny4*nz2*s1_dz 
- nx4_dz*ny2*nz3*s1 - nx4*ny2_dz*nz3*s1 - nx4*ny2*nz3_dz*s1 - nx4*ny2*nz3*s1_dz 
+ nx2_dz*ny4*nz3*s1 + nx2*ny4_dz*nz3*s1 + nx2*ny4*nz3_dz*s1 + nx2*ny4*nz3*s1_dz 
+ nx3_dz*ny2*nz4*s1 + nx3*ny2_dz*nz4*s1 + nx3*ny2*nz4_dz*s1 + nx3*ny2*nz4*s1_dz 
- nx2_dz*ny3*nz4*s1 - nx2*ny3_dz*nz4*s1 - nx2*ny3*nz4_dz*s1 - nx2*ny3*nz4*s1_dz 
- nx4_dz*ny3*nz1*s2 - nx4*ny3_dz*nz1*s2 - nx4*ny3*nz1_dz*s2 - nx4*ny3*nz1*s2_dz 
+ nx3_dz*ny4*nz1*s2 + nx3*ny4_dz*nz1*s2 + nx3*ny4*nz1_dz*s2 + nx3*ny4*nz1*s2_dz 
+ nx4_dz*ny1*nz3*s2 + nx4*ny1_dz*nz3*s2 + nx4*ny1*nz3_dz*s2 + nx4*ny1*nz3*s2_dz 
- nx1_dz*ny4*nz3*s2 - nx1*ny4_dz*nz3*s2 - nx1*ny4*nz3_dz*s2 - nx1*ny4*nz3*s2_dz 
- nx3_dz*ny1*nz4*s2 - nx3*ny1_dz*nz4*s2 - nx3*ny1*nz4_dz*s2 - nx3*ny1*nz4*s2_dz 
+ nx1_dz*ny3*nz4*s2 + nx1*ny3_dz*nz4*s2 + nx1*ny3*nz4_dz*s2 + nx1*ny3*nz4*s2_dz 
+ nx4_dz*ny2*nz1*s3 + nx4*ny2_dz*nz1*s3 + nx4*ny2*nz1_dz*s3 + nx4*ny2*nz1*s3_dz 
- nx2_dz*ny4*nz1*s3 - nx2*ny4_dz*nz1*s3 - nx2*ny4*nz1_dz*s3 - nx2*ny4*nz1*s3_dz 
- nx4_dz*ny1*nz2*s3 - nx4*ny1_dz*nz2*s3 - nx4*ny1*nz2_dz*s3 - nx4*ny1*nz2*s3_dz 
+ nx1_dz*ny4*nz2*s3 + nx1*ny4_dz*nz2*s3 + nx1*ny4*nz2_dz*s3 + nx1*ny4*nz2*s3_dz 
+ nx2_dz*ny1*nz4*s3 + nx2*ny1_dz*nz4*s3 + nx2*ny1*nz4_dz*s3 + nx2*ny1*nz4*s3_dz 
- nx1_dz*ny2*nz4*s3 - nx1*ny2_dz*nz4*s3 - nx1*ny2*nz4_dz*s3 - nx1*ny2*nz4*s3_dz 
- nx3_dz*ny2*nz1*s4 - nx3*ny2_dz*nz1*s4 - nx3*ny2*nz1_dz*s4 - nx3*ny2*nz1*s4_dz 
+ nx2_dz*ny3*nz1*s4 + nx2*ny3_dz*nz1*s4 + nx2*ny3*nz1_dz*s4 + nx2*ny3*nz1*s4_dz 
+ nx3_dz*ny1*nz2*s4 + nx3*ny1_dz*nz2*s4 + nx3*ny1*nz2_dz*s4 + nx3*ny1*nz2*s4_dz 
- nx1_dz*ny3*nz2*s4 - nx1*ny3_dz*nz2*s4 - nx1*ny3*nz2_dz*s4 - nx1*ny3*nz2*s4_dz 
- nx2_dz*ny1*nz3*s4 - nx2*ny1_dz*nz3*s4 - nx2*ny1*nz3_dz*s4 - nx2*ny1*nz3*s4_dz 
+ nx1_dz*ny2*nz3*s4 + nx1*ny2_dz*nz3*s4 + nx1*ny2*nz3_dz*s4 + nx1*ny2*nz3*s4_dz;

	 xins_dx = ( det * xins_dx - xins * det_dx ) / det / det ;

	 xins_dy = ( det * xins_dy - xins * det_dy ) / det / det ;

	 xins_dz = ( det * xins_dz - xins * det_dz ) / det / det ;

	 xins = xins / det;

	 dARdx[0] = ( circ * xins_dx - xins * circ_dx ) / circ / circ * 3.0;
	 dARdx[1] = ( circ * xins_dy - xins * circ_dy ) / circ / circ * 3.0;
	 dARdx[2] = ( circ * xins_dz - xins * circ_dz ) / circ / circ * 3.0;

	 *ar = xins/circ*3.0;

  return grid;
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

bool gridNegCellAroundNode( Grid *grid, int node )
{
  int cellId, nodes[4];
  AdjIterator it;

  for ( it = adjFirst(grid->cellAdj,node); adjValid(it); it = adjNext(it) ) {
    cellId = adjItem(it);
    gridCell( grid, cellId, nodes );
    if (gridVolume(grid, nodes) <= 0.0) return TRUE;
  }

  return FALSE;
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

int gridFindCellWithFace(Grid *grid, int face ){
  int n0, n1, n2;
  AdjIterator it0, it1, it2;
  Adj *adj=grid->cellAdj;

  if (face >= grid->maxface || face < 0 ) return EMPTY;
  if (grid->f2n[3*face] == EMPTY )return EMPTY;

  n0 = grid->f2n[0+3*face];
  n1 = grid->f2n[1+3*face];
  n2 = grid->f2n[2+3*face];

  for ( it0 = adjFirst(adj,n0); adjValid(it0); it0 = adjNext(it0) )
    for ( it1 = adjFirst(adj,n1); adjValid(it1); it1 = adjNext(it1) )
      if ( adjItem(it0) == adjItem(it1) )
	for ( it2 = adjFirst(adj,n2); adjValid(it2); it2 = adjNext(it2) )
	  if ( adjItem(it0)==adjItem(it2) ) return adjItem(it2);

  return EMPTY;
}

bool gridRightHandedFace(Grid *grid, int face ){
  int cell;
  int nodes[4];
  cell = gridFindCellWithFace(grid, face );
  if (cell == EMPTY) return FALSE;

  nodes[0] = grid->f2n[0+3*face];
  nodes[1] = grid->f2n[1+3*face];
  nodes[2] = grid->f2n[2+3*face];
  nodes[3] 
    = grid->c2n[0+4*cell] 
    + grid->c2n[1+4*cell] 
    + grid->c2n[2+4*cell] 
    + grid->c2n[3+4*cell]
    - nodes[0]
    - nodes[1]
    - nodes[2];

  return (gridVolume(grid, nodes) > 0.0);
}

bool gridRightHandedBoundary( Grid *grid )
{
  int face;

  for (face=0;face<grid->maxface;face++)
    if ( grid->f2n[3*face] != EMPTY )
      if ( !gridRightHandedFace(grid, face) ) return FALSE;

  return TRUE;
}

int gridNGeomNode(Grid *grid)
{
  return grid->nGeomNode;
}

Grid *gridSetNGeomNode(Grid *grid, int nGeomNode)
{
  grid->nGeomNode = nGeomNode;
  return grid;
}

bool gridGeometryNode(Grid *grid, int node)
{
  return (node<grid->nGeomNode);
}

bool gridGeometryEdge(Grid *grid, int node)
{
  AdjIterator it;

  for ( it = adjFirst(grid->edgeAdj,node); adjValid(it); it = adjNext(it) )
    return TRUE;
  
  return FALSE;
}

bool gridGeometryFace(Grid *grid, int node)
{
  AdjIterator it;

  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) )
    return TRUE;
  
  return FALSE;
}

