
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
  grid->nGeomNode = 0;
  grid->nGeomEdge = 0;
  grid->geomEdge = NULL;

  grid->xyz = malloc(3 * grid->maxnode * sizeof(double));
  for (i=0;i < grid->maxnode; i++ ) {
    grid->xyz[0+3*i] = DBL_MAX;
    grid->xyz[1+3*i] = (double)(i+1);
  }
  grid->xyz[1+3*(grid->maxnode-1)] = (double)(EMPTY);
  grid->blanknode = 0;

  grid->map = malloc(grid->maxnode * 6 * sizeof(double));
  for (i=0;i < grid->maxnode; i++ ) {
    grid->map[0+6*i] = 1.0;
    grid->map[1+6*i] = 0.0;
    grid->map[2+6*i] = 0.0;
    grid->map[3+6*i] = 1.0;
    grid->map[4+6*i] = 0.0;
    grid->map[5+6*i] = 1.0;
  }

  grid->frozen = malloc(grid->maxnode * sizeof(bool));
  for (i=0;i < grid->maxnode; i++ ) grid->frozen[i] = FALSE;

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
  grid->degAR = 0;

  grid->tecplotFileOpen = FALSE;

  grid->renumberFunc = NULL;
  grid->renumberData = NULL;

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
  grid->nGeomNode = 0;
  grid->nGeomEdge = 0;
  grid->geomEdge = NULL;

  grid->xyz = xyz;
  for (i=grid->nnode ; i < grid->maxnode; i++ ) {
    grid->xyz[0+3*i] = DBL_MAX;
    grid->xyz[1+3*i] = (double)(i+1);
  }
  if (grid->maxnode == grid->nnode) {
    grid->blanknode = EMPTY;
  }else{
    grid->xyz[1+3*(grid->maxnode-1)] = (double)(EMPTY);
    grid->blanknode = grid->nnode;
  }

  grid->map = malloc(grid->maxnode * 6 * sizeof(double));
  for (i=0;i < grid->maxnode; i++ ) {
    grid->map[0+6*i] = 1.0;
    grid->map[1+6*i] = 0.0;
    grid->map[2+6*i] = 0.0;
    grid->map[3+6*i] = 1.0;
    grid->map[4+6*i] = 0.0;
    grid->map[5+6*i] = 1.0;
  }

  grid->frozen = malloc(grid->maxnode * sizeof(bool));
  for (i=0;i < grid->maxnode; i++ ) grid->frozen[i] = FALSE;

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
  grid->faceId = faceId;

  for ( i=grid->nface ; i < grid->maxface ; i++ ) {
    grid->f2n[0+3*i] = EMPTY; 
    grid->f2n[1+3*i] = i+1; 
    grid->faceId[i] = EMPTY; 
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
  grid->degAR = 0;

  grid->tecplotFileOpen = FALSE;

  grid->renumberFunc = NULL;
  grid->renumberData = NULL;

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

  if (NULL == gridPack(grid)) {
    printf("gridExport, gridPack failed.\n");
    return NULL;
  }

  *nnode = grid->nnode;
  *ncell = grid->ncell;
  *nface = grid->nface;

  *xyz = grid->xyz;

  *c2n = grid->c2n;

  *f2n    = grid->f2n;
  *faceId = grid->faceId;

  return  grid;
}

Grid *gridExportFAST( Grid *grid, char *filename )
{
  FILE *file;
  int i;

  if (NULL == gridPack(grid)) {
    printf("gridExportFAST: gridPack failed.\n");
    return NULL;
  }

  printf("gridExportFAST: open file: %s\n",filename);
  file = fopen(filename,"w");

  fprintf(file,"%10d %10d %10d\n",grid->nnode,grid->nface,grid->ncell);

  printf("gridExportFAST: writing xyz...\n");
  
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[0+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[1+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[2+3*i]);

  printf("gridExportFAST: writing faces...\n");

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%10d %10d %10d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%4d",grid->faceId[i]);
  }

  printf("gridExportFAST: writing cells...\n");
  
  for( i=0; i<grid->ncell ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);

  }

  fclose(file);

  return grid;
}

Grid *gridImportAdapt( Grid *grid, char *filename )
{
  int i;
  FILE *file;

  file = fopen(filename,"r");

  for( i=0; i<grid->nnode ; i++ ) 
    fscanf(file,"%lf %lf %lf %lf %lf %lf\n",
       &grid->map[0+6*i], &grid->map[1+6*i], &grid->map[2+6*i],
                          &grid->map[3+6*i], &grid->map[4+6*i],
                                             &grid->map[5+6*i]);
  fclose(file);
  return grid;
}

Grid *gridAttachNodeSorter(Grid *grid, 
			   void (*renumberFunc)(void *renumberData, int *o2n), 
			   void *renumberData )
{
  grid->renumberFunc = renumberFunc;
  grid->renumberData = renumberData;
  return grid;
}

Grid *gridDetachNodeSorter(Grid *grid )
{
  grid->renumberFunc = NULL;
  grid->renumberData = NULL;
  return grid;
}

void gridFree(Grid *grid)
{
  if ( grid->tecplotFileOpen ) fclose(grid->tecplotFile);
  if ( NULL != grid->geomEdge) free(grid->geomEdge);
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
  free(grid->frozen);
  free(grid->map);
  free(grid->xyz);
  free(grid);
}

Grid *gridPack(Grid *grid)
{
  int i;
  int orignode, packnode, origcell, packcell;
  int origface, packface, origedge, packedge;
  int iface, n0, n1, id, n[3];
  double t0, t1, u[3], v[3];
  bool emptyFace;
  int *o2n;

  o2n = malloc(grid->maxnode*sizeof(int));
  if (o2n == NULL) {
    printf("ERROR: Pack: %s: %d: could not allocate o2n\n",__FILE__,__LINE__);
    return NULL;
  }
  for (i=0;i<grid->maxnode;i++) o2n[i] = EMPTY;

  packnode = 0;
  for ( orignode=0 ; orignode < grid->maxnode ; orignode++ )
    if ( grid->xyz[0+3*orignode] != DBL_MAX) {
      o2n[orignode] = packnode;
      grid->xyz[0+3*packnode] = grid->xyz[0+3*orignode];
      grid->xyz[1+3*packnode] = grid->xyz[1+3*orignode];
      grid->xyz[2+3*packnode] = grid->xyz[2+3*orignode];
      grid->map[0+6*packnode] = grid->map[0+6*orignode];
      grid->map[1+6*packnode] = grid->map[1+6*orignode];
      grid->map[2+6*packnode] = grid->map[2+6*orignode];
      grid->map[3+6*packnode] = grid->map[3+6*orignode];
      grid->map[4+6*packnode] = grid->map[4+6*orignode];
      grid->map[5+6*packnode] = grid->map[5+6*orignode];
      grid->frozen[packnode]  = grid->frozen[orignode];
      packnode++;
    } 
  
  if (grid->nnode != packnode) {
    printf("ERROR: Pack: %s: %d: grid->nnode != packnode\n",
	   __FILE__, __LINE__ );
    return NULL;
  }

  for (i=grid->nnode ; i < grid->maxnode; i++ ) {
    grid->xyz[0+3*i] = DBL_MAX;
    grid->xyz[1+3*i] = (double)(i+1);
  }
  if (grid->maxnode == grid->nnode) {
    grid->blanknode = EMPTY;
  }else{
    grid->xyz[1+3*(grid->maxnode-1)] = (double)(EMPTY);
    grid->blanknode = grid->nnode;
  }

  packcell = 0;
  for ( origcell=0 ; origcell < grid->maxcell ; origcell++ )
    if ( grid->c2n[0+4*origcell] != EMPTY) {
      adjRemove( grid->cellAdj, grid->c2n[0+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[1+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[2+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[3+4*origcell], origcell );
      grid->c2n[0+4*packcell] = o2n[grid->c2n[0+4*origcell]];
      grid->c2n[1+4*packcell] = o2n[grid->c2n[1+4*origcell]];
      grid->c2n[2+4*packcell] = o2n[grid->c2n[2+4*origcell]];
      grid->c2n[3+4*packcell] = o2n[grid->c2n[3+4*origcell]];
      adjRegister( grid->cellAdj, grid->c2n[0+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[1+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[2+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[3+4*packcell], packcell );
      packcell++;
    } 
  
  if (grid->ncell != packcell) {
    printf("ERROR: Pack: %s: %d: grid->ncell %d != packcell %d \n",
	   __FILE__, __LINE__, grid->ncell, packcell);
    return NULL;
  }

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

  packface=0;

  emptyFace = FALSE;
  iface=0;
  while (!emptyFace) {
    emptyFace = TRUE;
    iface++;
    for ( origface=0 ; origface < grid->maxface ; origface++ ){ 
      if ( grid->faceId[origface]==iface && grid->f2n[0+3*origface] != EMPTY) {
	emptyFace = FALSE;
	if (origface == packface) {
	  for (i=0;i<3;i++){
	    adjRemove( grid->faceAdj, grid->f2n[i+3*packface], packface );
	    grid->f2n[i+3*packface] = o2n[grid->f2n[i+3*packface]];
	    adjRegister( grid->faceAdj, grid->f2n[i+3*packface], packface );
	  }
	}else{
	  for (i=0;i<3;i++) {
	    n[i] = grid->f2n[i+3*origface];
	    adjRemove( grid->faceAdj, grid->f2n[i+3*origface], origface );
	    u[i] = grid->faceU[i+3*origface];
	    v[i] = grid->faceV[i+3*origface];
	  }
	  id = grid->faceId[origface];
	  
	  for (i=0;i<3;i++) {
	    grid->f2n[i+3*origface] = grid->f2n[i+3*packface];
	    if ( grid->f2n[3*packface] != EMPTY ) {
	      adjRemove( grid->faceAdj, grid->f2n[i+3*packface], packface );
	      adjRegister( grid->faceAdj, grid->f2n[i+3*origface], origface );
	    }
	    grid->faceU[i+3*origface] = grid->faceU[i+3*packface];
	    grid->faceV[i+3*origface] = grid->faceV[i+3*packface];
	  }

	  grid->faceId[origface]	= grid->faceId[packface];
	  
	  for (i=0;i<3;i++) {
	    grid->f2n[i+3*packface] = o2n[n[i]];
	    adjRegister( grid->faceAdj, grid->f2n[i+3*packface], packface );
	    grid->faceU[i+3*packface] = u[i];
	    grid->faceV[i+3*packface] = v[i];
	  }
	  grid->faceId[packface]  = id;
	}
	packface++;
      } 
    }
  }
  //iface--; printf("gridPack: %d geometry faces detected.\n",iface);

  if (grid->nface != packface) {
    printf("ERROR: Pack: %s: %d: grid->nface %d != packface %d \n",
	   __FILE__, __LINE__, grid->nface, packface );
    return NULL;
  }

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

  packedge = 0;
  for (origedge=0;origedge<grid->maxedge;origedge++){
    if (grid->e2n[0+2*origedge] != EMPTY){

      n0 = o2n[grid->e2n[0+2*origedge]];
      n1 = o2n[grid->e2n[1+2*origedge]];
      if (n0==EMPTY || n1==EMPTY) {
	printf("ERROR: Pack: %s: %d: n0 %d n1 %d: edge orig %d pack %d\n",
	       __FILE__, __LINE__, n0, n1, origedge, packedge);
	return NULL;
      }

      id = grid->edgeId[origedge];
      t0 = grid->edgeT[0+2*origedge];
      t1 = grid->edgeT[1+2*origedge];
      adjRemove( grid->edgeAdj, grid->e2n[0+2*origedge], origedge );
      adjRemove( grid->edgeAdj, grid->e2n[1+2*origedge], origedge );

      grid->e2n[0+2*packedge] = n0;
      grid->e2n[1+2*packedge] = n1;
      grid->edgeId[packedge] = id;
      grid->edgeT[0+2*packedge] = t0;
      grid->edgeT[1+2*packedge] = t1;
      adjRegister( grid->edgeAdj, n0, packedge );
      adjRegister( grid->edgeAdj, n1, packedge );
     
      packedge++;
    }
  }

  if (grid->nedge != packedge) {
    printf("ERROR: Pack: %s: %d: grid->nedge %d != packedge %d\n",
	   __FILE__, __LINE__, grid->nface, packface);
    return NULL;
  }

  for (i=grid->nedge;i < grid->maxedge; i++ ) {
    grid->e2n[0+2*i] = EMPTY; 
    grid->e2n[1+2*i] = i+1; 
  }
  if (grid->maxedge == grid->nedge) {
    grid->blanke2n = EMPTY;
  }else{
    grid->e2n[1+2*(grid->maxedge-1)] = EMPTY; 
    grid->blanke2n = grid->nedge;
  }

  if ( NULL != grid->renumberFunc ) 
    (*grid->renumberFunc)( grid->renumberData, o2n );

  free(o2n);

  return  grid;
}

Grid *gridSortNodeGridEx(Grid *grid)
{
  int i, newnode, edge, nCurveNode;
  int ixyz, node, face, cell, inode;
  int *o2n, *curve;
  double *temp_xyz;
  bool *temp_frozen;

  if (NULL == gridPack(grid)) {
    printf("gridSortNodeGridEx: gridPack failed.\n");
    return NULL;
  }

  o2n = malloc( grid->nnode * sizeof(int) );
  if (o2n == NULL) {
    printf("ERROR: gridSortNodeGridEx: %s: %d: could not allocate o2n\n",
	   __FILE__,__LINE__);
    return NULL;
  }
  for (i=0;i<grid->nnode;i++) o2n[i] = EMPTY;

  // geom nodes
  for (i=0;i<grid->nGeomNode;i++) o2n[i] = i;
  newnode = grid->nGeomNode;

  // edge stuff
  for (edge=1; edge<=grid->nGeomEdge; edge++){

    nCurveNode = gridGeomEdgeSize( grid, edge );
    curve = malloc( nCurveNode * sizeof(int) );
    gridGeomEdge( grid, edge, curve );

    for ( i=1; i<(nCurveNode-1); i++){ // skip end points
      if (o2n[curve[i]] != EMPTY) 
	printf("gridSortNodeGridEx: %s: %d: newnode error %d\n",
	       __FILE__, __LINE__, o2n[curve[i]] );
      o2n[curve[i]] = newnode;
      newnode++;
    }

    free(curve);
  }

  // face stuff - assuming that the bc faces are sorted.
  for ( face=0; face<grid->nface; face++ ){
    for ( i=0; i<3; i++ ){
      node = grid->f2n[i+3*face];
      if ( o2n[node] == EMPTY ) {
	o2n[node] = newnode;
	newnode++;
      }
    }
  }

  // interior nodes
  for ( node=0; node<grid->nnode; node++ ){
    if ( o2n[node] == EMPTY ) {
      o2n[node] = newnode;
      newnode++;
    }
  }

  if (newnode != grid->nnode) 
    printf("ERROR: gridSortNodeGridEx, newnode %d nnode %d, line %d of %s\n.",
	   newnode,grid->nnode,__LINE__, __FILE__);

  temp_xyz = malloc( grid->nnode * sizeof(double) );

  if (temp_xyz == NULL) {
    printf("ERROR: gridSortNodeGridEx: %s: %d: could not allocate temp_xyz\n",
	   __FILE__,__LINE__);
    return NULL;
  }

  for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
    for ( node = 0 ; node < grid->nnode ; node++ ){
      temp_xyz[o2n[node]] = grid->xyz[ixyz+3*node];
    }
    for ( node = 0 ; node < grid->nnode ; node++ ){
      grid->xyz[ixyz+3*node] = temp_xyz[node];
    }
  }
  for ( ixyz = 0; ixyz < 6 ; ixyz++ ){
    for ( node = 0 ; node < grid->nnode ; node++ ){
      temp_xyz[o2n[node]] = grid->map[ixyz+6*node];
    }
    for ( node = 0 ; node < grid->nnode ; node++ ){
      grid->map[ixyz+6*node] = temp_xyz[node];
    }
  }

  free(temp_xyz);

  temp_frozen = malloc( grid->nnode * sizeof(bool) );

  for ( node = 0 ; node < grid->nnode ; node++ ){
    temp_frozen[o2n[node]] = grid->frozen[node];
  }
  for ( node = 0 ; node < grid->nnode ; node++ ){
    grid->frozen[node] = temp_frozen[node];
  }

  free(temp_frozen);

  for ( cell = 0; cell < grid->ncell ; cell++ ){
    for ( inode = 0 ; inode < 4 ; inode++ ){
      adjRemove( grid->cellAdj, grid->c2n[inode+4*cell], cell );
      grid->c2n[inode+4*cell] = o2n[grid->c2n[inode+4*cell]];
      adjRegister( grid->cellAdj, grid->c2n[inode+4*cell], cell );
    }
  }

  for ( face = 0; face < grid->nface ; face++ ){
    for ( inode = 0 ; inode < 3 ; inode++ ){
      adjRemove( grid->faceAdj, grid->f2n[inode+3*face], face );
      grid->f2n[inode+3*face] = o2n[grid->f2n[inode+3*face]];
      adjRegister( grid->faceAdj, grid->f2n[inode+3*face], face );
    }
  }

  for ( edge = 0; edge < grid->nedge ; edge++ ){
    for ( inode = 0 ; inode < 2 ; inode++ ){
      adjRemove( grid->edgeAdj, grid->e2n[inode+2*edge], edge );
      grid->e2n[inode+2*edge] = o2n[grid->e2n[inode+2*edge]];
      adjRegister( grid->edgeAdj, grid->e2n[inode+2*edge], edge );
    }
  }

  if ( NULL != grid->renumberFunc ) 
    (*grid->renumberFunc)( grid->renumberData, o2n );

  free(o2n);

  return grid;
}

Grid *gridWriteTecplotSurfaceZone(Grid *grid)
{
  int i, nfacenode;
  if ( grid !=  gridSortNodeGridEx(grid) ) {
    printf("gridWriteTecplotSurfaceZone: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (!grid->tecplotFileOpen) {
    grid->tecplotFile = fopen("grid.t","w");
    grid->tecplotFileOpen = TRUE;
    fprintf(grid->tecplotFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotFile, "variables=\"X\",\"Y\",\"Z\"\n");
  }

  nfacenode=0;
  for(i=0;i<3*grid->nface;i++){
    nfacenode = MAX(nfacenode, grid->f2n[i]);
  }
  nfacenode++;

  fprintf(grid->tecplotFile, "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  nfacenode, grid->nface);

  for ( i=0; i<nfacenode ; i++ ){
    fprintf(grid->tecplotFile, "%23.15e%23.15e%23.15e\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i]);
  }

  fprintf(grid->tecplotFile, "\n");

  for ( i=0; i<grid->nface ; i++ ){
    fprintf(grid->tecplotFile, " %9d %9d %9d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  fflush(grid->tecplotFile);

  return grid;
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

Grid *gridReconnectCell(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int cell, i, node;

  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->cellAdj,oldNode);
  while (adjValid(it)){
    cell = adjItem(it);
    for (i=0;i<4;i++){
      node = grid->c2n[i+4*cell];
      if (oldNode == node) {
	grid->c2n[i+4*cell]=newNode;
	adjRemove( grid->cellAdj, oldNode, cell);
	adjRegister( grid->cellAdj, newNode, cell);
      }
    }
    it = adjFirst(grid->cellAdj,oldNode);
  }

  return grid;
}

Grid *gridReconnectCellUnlessFrozen(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int cell, i, node;
  bool frozen;

  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;
  
  it = adjFirst(grid->cellAdj,oldNode);
  while (adjValid(it)){
    cell = adjItem(it);
    frozen = ( gridNodeFrozen(grid, grid->c2n[0+4*cell]) &&
	       gridNodeFrozen(grid, grid->c2n[1+4*cell]) &&
	       gridNodeFrozen(grid, grid->c2n[2+4*cell]) &&
	       gridNodeFrozen(grid, grid->c2n[3+4*cell]) );
    if (!frozen ) {
      for (i=0;i<4;i++){
	node = grid->c2n[i+4*cell];
	if (oldNode == node) {
	  grid->c2n[i+4*cell]=newNode;
	  adjRemove( grid->cellAdj, oldNode, cell);
	  adjRegister( grid->cellAdj, newNode, cell);
	}
      }
      it = adjFirst(grid->cellAdj,oldNode);
    }else{
      it = adjNext(it);
    }      
  }
  
  return grid;
}

Grid *gridCell(Grid *grid, int cellId, int *nodes )
{
  if ( cellId < 0 || cellId >= grid->maxcell ) return NULL;
  if ( grid->c2n[4*cellId] == EMPTY ) return NULL;

  nodes[0] = grid->c2n[0+4*cellId];
  nodes[1] = grid->c2n[1+4*cellId];
  nodes[2] = grid->c2n[2+4*cellId];
  nodes[3] = grid->c2n[3+4*cellId];

  return grid;
}

bool gridCellEdge(Grid *grid, int n0, int n1 )
{
  AdjIterator it;
  int i, nodes[4];

  if ( n0 < 0 || n0 >= grid->maxnode ) return FALSE;
  if ( n1 < 0 || n1 >= grid->maxnode ) return FALSE;

  for ( it = adjFirst(grid->cellAdj,n0); adjValid(it); it = adjNext(it) ) {
    gridCell( grid, adjItem(it), nodes );
    for(i=0;i<4;i++) if (n1 == nodes[i]) return TRUE;
  }

  return FALSE;
}

bool gridCellFace(Grid *grid, int n0, int n1, int n2 )
{
  AdjIterator it;
  int i, j, nodes[4];

  if ( n0 < 0 || n0 >= grid->maxnode ) return FALSE;
  if ( n1 < 0 || n1 >= grid->maxnode ) return FALSE;
  if ( n2 < 0 || n2 >= grid->maxnode ) return FALSE;

  for ( it = adjFirst(grid->cellAdj,n0); adjValid(it); it = adjNext(it) ) {
    gridCell( grid, adjItem(it), nodes );
    for(i=0;i<3;i++) if (n0 == nodes[i]) for(j=i;j<3;j++) nodes[i]=nodes[i+1];
    if ( ( n1 == nodes[0] || n1 == nodes[1] || n1 == nodes[2] ) &&
	 ( n2 == nodes[0] || n2 == nodes[1] || n2 == nodes[2] ) ) return TRUE;
  }

  return FALSE;
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

Grid *gridReconnectFace(Grid *grid, int faceId, int oldNode, int newNode )
{
  AdjIterator it;
  int face, i, node;
  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->faceAdj,oldNode);
  while (adjValid(it)){
    face = adjItem(it);
    if (faceId == grid->faceId[face]) {
      for (i=0;i<3;i++){
	node = grid->f2n[i+3*face];
	if (oldNode == node) {
	  grid->f2n[i+3*face]=newNode;
	  adjRemove( grid->faceAdj, oldNode, face);
	  adjRegister( grid->faceAdj, newNode, face);
	}
      }
      it = adjFirst(grid->faceAdj,oldNode);
    }else{
      it = adjNext(it);
    }
  }

  return grid;
}

Grid *gridReconnectFaceUnlessFrozen(Grid *grid, int faceId, 
				    int oldNode, int newNode )
{
  AdjIterator it;
  int face, i, node;
  bool frozen;

  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->faceAdj,oldNode);
  while (adjValid(it)){
    face = adjItem(it);
    frozen = ( gridNodeFrozen(grid, grid->f2n[0+3*face]) &&
	       gridNodeFrozen(grid, grid->f2n[1+3*face]) &&
	       gridNodeFrozen(grid, grid->f2n[2+3*face]) );
    if (faceId == grid->faceId[face] && !frozen ) {
      for (i=0;i<3;i++){
	node = grid->f2n[i+3*face];
	if (oldNode == node) {
	  grid->f2n[i+3*face]=newNode;
	  adjRemove( grid->faceAdj, oldNode, face);
	  adjRegister( grid->faceAdj, newNode, face);
	}
      }
      it = adjFirst(grid->faceAdj,oldNode);
    }else{
      it = adjNext(it);
    }
  }

  return grid;
}

Grid *gridFace(Grid *grid, int face, int *nodes, int *id )
{
  if (face >= grid->maxface || face < 0) return NULL;
  if (EMPTY == grid->f2n[3*face]) return NULL;

  nodes[0] = grid->f2n[0+3*face];
  nodes[1] = grid->f2n[1+3*face];
  nodes[2] = grid->f2n[2+3*face];
  *id = grid->faceId[face];
  return grid;
}

Grid *gridNodeUV(Grid *grid, int node, int faceId, double *uv )
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

Grid *gridReconnectEdge(Grid *grid, int edgeId, int oldNode, int newNode )
{
  AdjIterator it;
  int edge, i, node;
  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->edgeAdj,oldNode);
  while (adjValid(it)){
    edge = adjItem(it);
    if (edgeId == grid->edgeId[edge]) {
      for (i=0;i<2;i++){
	node = grid->e2n[i+2*edge];
	if (oldNode == node) {
	  grid->e2n[i+2*edge]=newNode;
	  adjRemove( grid->edgeAdj, oldNode, edge);
	  adjRegister( grid->edgeAdj, newNode, edge);
	}
      }
      it = adjFirst(grid->edgeAdj,oldNode);
    }else{
      it = adjNext(it);
    }
  }

  return grid;
}

Grid *gridReconnectEdgeUnlessFrozen(Grid *grid, int edgeId, 
				    int oldNode, int newNode )
{
  AdjIterator it;
  int edge, i, node;
  bool frozen;

  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->edgeAdj,oldNode);
  while (adjValid(it)){
    edge = adjItem(it);
    frozen = ( gridNodeFrozen(grid, grid->e2n[0+2*edge]) &&
	       gridNodeFrozen(grid, grid->e2n[1+2*edge]) );
    if (edgeId == grid->edgeId[edge] && !frozen ) {
      for (i=0;i<2;i++){
	node = grid->e2n[i+2*edge];
	if (oldNode == node) {
	  grid->e2n[i+2*edge]=newNode;
	  adjRemove( grid->edgeAdj, oldNode, edge);
	  adjRegister( grid->edgeAdj, newNode, edge);
	}
      }
      it = adjFirst(grid->edgeAdj,oldNode);
    }else{
      it = adjNext(it);
    }
  }

  return grid;
}

Grid *gridEdge(Grid *grid, int edge, int *nodes, int *id )
{
  if (edge >= grid->maxedge || edge < 0) return NULL;
  if (EMPTY == grid->e2n[2*edge]) return NULL;

  nodes[0] = grid->e2n[0+2*edge];
  nodes[1] = grid->e2n[1+2*edge];
  *id = grid->edgeId[edge];
  return grid;
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

int gridNFrozen( Grid *grid )
{
  int node, nfrozen;
  nfrozen =0;
  for ( node=0 ; node<grid->maxnode ; node++ )
    if ( gridValidNode( grid, node) && gridNodeFrozen( grid, node ) ) 
      nfrozen++; 
  return nfrozen;
}

bool gridNodeFrozen( Grid *grid, int node )
{
  if ( !gridValidNode( grid, node) ) return TRUE;
  return grid->frozen[node];
}

Grid *gridFreezeNode( Grid *grid, int node )
{
  if ( !gridValidNode( grid, node) ) return NULL;
  grid->frozen[node] = TRUE;
  return grid;
}

Grid *gridThawNode( Grid *grid, int node )
{
  if ( !gridValidNode( grid, node) ) return NULL;
  grid->frozen[node] = FALSE;
  return grid;
}

Grid *gridFreezeAll( Grid *grid )
{
  int node;
  for ( node = 0; node<grid->maxnode; node++ ) grid->frozen[node] = TRUE;
  return grid;
}

Grid *gridThawAll( Grid *grid )
{
  int node;
  for ( node = 0; node<grid->maxnode; node++ ) grid->frozen[node] = FALSE;
  return grid;
}

Grid *gridThawSphere( Grid *grid, double x, double y, double z, double r )
{
  int node;
  double dx, dy, dz, distanceSquared, radiusSquared;
  radiusSquared = r*r;
  
  for ( node=0; node<grid->nnode; node++ ) {
    dx = grid->xyz[0+3*node] - x;
    dy = grid->xyz[1+3*node] - y;
    dz = grid->xyz[2+3*node] - z;
    distanceSquared = dx*dx + dy*dy + dz*dz;
    if (radiusSquared >= distanceSquared) gridThawNode(grid, node );
  }

  return grid;
}

Grid *gridFreezeBCFace( Grid *grid, int faceId )
{
  int face, node, i;
  for ( face=0; face<grid->nface; face++ ){
    if (grid->faceId[face] == faceId ){
      for ( i=0; i<3; i++ ){
	node = grid->f2n[i+3*face];
	gridFreezeNode(grid,node);
      }
    }
  }
  return grid;
}

Grid *gridThawNearBC( Grid *grid, double r, int faceId )
{
  int face, node, i;
  for ( face=0; face<grid->nface; face++ ){
    if (grid->faceId[face] == faceId ){
      for ( i=0; i<3; i++ ){
	node = grid->f2n[i+3*face];
	gridThawSphere(grid, grid->xyz[0+3*node], 
		       grid->xyz[1+3*node], grid->xyz[2+3*node], r);
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

int gridAddNode(Grid *grid, double x, double y, double z )
{
  int node;
  if (grid->nnode >= grid->maxnode) return EMPTY;
  node = grid->blanknode;
  if (EMPTY == node) return EMPTY;
  grid->blanknode = (int)grid->xyz[1+3*node];
  grid->nnode++;

  grid->xyz[0+3*node] = x;
  grid->xyz[1+3*node] = y;
  grid->xyz[2+3*node] = z;
  grid->frozen[node] = FALSE;

  return node;
}

Grid *gridRemoveNode(Grid *grid, int node )
{
  if (node>grid->maxnode) return NULL;
  if (DBL_MAX == grid->xyz[0+3*node]) return NULL;
  grid->nnode--;
  grid->xyz[0+3*node] = DBL_MAX;
  grid->xyz[1+3*node] = (double)grid->blanknode;
  grid->blanknode = node;
  return grid;
}

bool gridValidNode(Grid *grid, int node )
{
  if (node < 0 || node >=grid->maxnode) return FALSE;
  return (DBL_MAX != grid->xyz[0+3*node]);
}

Grid *gridNodeXYZ(Grid *grid, int node, double *xyz )
{
  if (!gridValidNode(grid,node)) return NULL;
  xyz[0] = grid->xyz[0+3*node];
  xyz[1] = grid->xyz[1+3*node];
  xyz[2] = grid->xyz[2+3*node];
  return grid;
}

Grid *gridSetNodeXYZ(Grid *grid, int node, double *xyz )
{
  if (!gridValidNode(grid,node)) return NULL;
  grid->xyz[0+3*node] = xyz[0];
  grid->xyz[1+3*node] = xyz[1];
  grid->xyz[2+3*node] = xyz[2];
  return grid;
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

int gridNGeomNode(Grid *grid)
{
  return grid->nGeomNode;
}

Grid *gridSetNGeomNode(Grid *grid, int nGeomNode)
{
  grid->nGeomNode = nGeomNode;
  return grid;
}

int gridNGeomEdge(Grid *grid)
{
  return grid->nGeomEdge;
}

Grid *gridSetNGeomEdge(Grid *grid, int nGeomEdge)
{
  grid->nGeomEdge = nGeomEdge;
  if ( NULL != grid->geomEdge) free(grid->geomEdge);
  grid->geomEdge = malloc(2*nGeomEdge*sizeof(int));
  return grid;
}

Grid *gridAddGeomEdge(Grid *grid, int edge, int n0, int n1 )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return NULL;
  grid->geomEdge[0+(edge-1)*2] = n0;
  grid->geomEdge[1+(edge-1)*2] = n1;
  return grid;
}

int gridGeomEdgeSize( Grid *grid, int edge )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return EMPTY;
  return gridGeomCurveSize( grid, edge, grid->geomEdge[0+2*(edge-1)]);
}

Grid *gridGeomEdge( Grid *grid, int edge, int *curve )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return NULL;
  return gridGeomCurve( grid, edge, grid->geomEdge[0+2*(edge-1)], curve );
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

