
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
#include <limits.h>
#include <values.h>
#include "grid.h"

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
  grid->nprism=0;
  grid->maxprism=0;
  grid->prism=NULL;

  grid->prismDeg=NULL;

  grid->npyramid=0;
  grid->maxpyramid=0;
  grid->pyramid=NULL;
  grid->nquad=0;
  grid->maxquad=0;
  grid->quad=NULL;
  grid->nGeomNode = 0;
  grid->nGeomEdge = 0;
  grid->nGeomFace = 0;
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

  grid->lines = linesCreate();

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
  grid->nprism=0;
  grid->maxprism=0;
  grid->prism=NULL;
  grid->prismDeg=NULL;
  grid->npyramid=0;
  grid->maxpyramid=0;
  grid->pyramid=NULL;
  grid->nquad=0;
  grid->maxquad=0;
  grid->quad=NULL;
  grid->nGeomNode = 0;
  grid->nGeomEdge = 0;
  grid->nGeomFace = 0;
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

  grid->lines = linesCreate();

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
    fprintf(file,"%4d\n",grid->faceId[i]);
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

Grid *gridExportAFLR3( Grid *grid, char *filename )
{
  FILE *file;
  int i;

  if (NULL == gridPack(grid)) {
    printf("gridExportAFLR3: gridPack failed.\n");
    return NULL;
  }

  printf("gridExportAFLR3: open file: %s\n",filename);
  file = fopen(filename,"w");

  fprintf(file,"%10d %10d %10d %10d %10d %10d %10d\n",
	  grid->nnode,grid->nface,grid->nquad,
	  grid->ncell,grid->npyramid,grid->nprism,0);

  printf("gridExportAFLR3: writing xyz...\n");
  
  for( i=0; i<grid->nnode ; i++ ) 
    fprintf(file,"%25.15e %25.15e %25.15e\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i]);

  printf("gridExportAFLR3: writing faces...\n");

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%10d %10d %10d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  for( i=0; i<grid->nquad ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->quad[i].nodes[0]+1,grid->quad[i].nodes[1]+1,
	    grid->quad[i].nodes[2]+1,grid->quad[i].nodes[3]+1);
  }

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%d\n",grid->faceId[i]);
  }

  for( i=0; i<grid->nquad ; i++ ) {
    fprintf(file,"%d\n",grid->quad[i].faceId);
  }

  printf("gridExportAFLR3: writing cells...\n");
  
  for( i=0; i<grid->ncell ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);

  }

  for( i=0; i<grid->npyramid ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d %10d\n",
	    grid->pyramid[i].nodes[0]+1,
	    grid->pyramid[i].nodes[1]+1,
	    grid->pyramid[i].nodes[2]+1,
	    grid->pyramid[i].nodes[3]+1,
	    grid->pyramid[i].nodes[4]+1);
  }

  for( i=0; i<grid->nprism ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d %10d %10d\n",
	    grid->prism[i].nodes[0]+1,
	    grid->prism[i].nodes[1]+1,
	    grid->prism[i].nodes[2]+1,
	    grid->prism[i].nodes[3]+1,
	    grid->prism[i].nodes[4]+1,
	    grid->prism[i].nodes[5]+1);
  }

  /* ain't got no hexes */

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
  if ( NULL != grid->lines ) linesFree(grid->lines);

  if ( grid->tecplotFileOpen ) fclose(grid->tecplotFile);
  if ( NULL != grid->geomEdge) free(grid->geomEdge);

  if (grid->prismDeg!=NULL) free(grid->prismDeg);

  if (grid->quad!=NULL) free(grid->quad);
  if (grid->pyramid!=NULL) free(grid->pyramid);
  if (grid->prism!=NULL) free(grid->prism);
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
  int nFaceId;
  int iface, n0, n1, id, n[3];
  double t0, t1, u[3], v[3];
  bool emptyFace;
  int *o2n;
  int prismIndex, pyramidIndex, quadIndex;

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

  nFaceId = 0;
  for ( origface=0 ; origface < grid->maxface ; origface++ ) 
    if (grid->f2n[0+3*origface] != EMPTY)
      nFaceId = MAX(nFaceId, grid->faceId[origface]);

  emptyFace = FALSE;
  iface=0;
  for (iface=1;iface <= nFaceId;iface++){
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

  for (prismIndex=0;prismIndex<gridNPrism(grid);prismIndex++){
    for (i=0;i<6;i++) grid->prism[prismIndex].nodes[i] =
			o2n[grid->prism[prismIndex].nodes[i]];
  }

  for (pyramidIndex=0;pyramidIndex<gridNPyramid(grid);pyramidIndex++){
    for (i=0;i<5;i++) grid->pyramid[pyramidIndex].nodes[i] =
			o2n[grid->pyramid[pyramidIndex].nodes[i]];
  }

  // note, these should be counted as boundaries
  for (quadIndex=0;quadIndex<gridNQuad(grid);quadIndex++){
    for (i=0;i<4;i++) grid->quad[quadIndex].nodes[i] =
			o2n[grid->quad[quadIndex].nodes[i]];
  }

  if ( NULL != grid->renumberFunc ) 
    (*grid->renumberFunc)( grid->renumberData, o2n );

  if ( NULL != gridLines(grid) ) linesRenumber(gridLines(grid),o2n);

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
  int prismIndex, pyramidIndex, quadIndex;

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
  newnode = MAX(0,grid->nGeomNode);

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

  for (prismIndex=0;prismIndex<gridNPrism(grid);prismIndex++){
    for (i=0;i<6;i++) grid->prism[prismIndex].nodes[i] =
			o2n[grid->prism[prismIndex].nodes[i]];
  }

  for (pyramidIndex=0;pyramidIndex<gridNPyramid(grid);pyramidIndex++){
    for (i=0;i<5;i++) grid->pyramid[pyramidIndex].nodes[i] =
			o2n[grid->pyramid[pyramidIndex].nodes[i]];
  }

  // note, these should be counted as boundaries
  for (quadIndex=0;quadIndex<gridNQuad(grid);quadIndex++){
    for (i=0;i<4;i++) grid->quad[quadIndex].nodes[i] =
			o2n[grid->quad[quadIndex].nodes[i]];
  }

  if ( NULL != grid->renumberFunc ) 
    (*grid->renumberFunc)( grid->renumberData, o2n );

  if ( NULL != gridLines(grid) ) linesRenumber(gridLines(grid),o2n);

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

int gridNPrism(Grid *grid)
{
  return grid->nprism;
}

int gridNPyramid(Grid *grid)
{
  return grid->npyramid;
}

int gridNQuad(Grid *grid)
{
  return grid->nquad;
}

int gridNGem(Grid *grid)
{
  return grid->ngem;
}

int gridGem(Grid *grid, int index)
{
  return grid->gem[index];
}

Grid *gridRemoveGem(Grid *grid) {
  int i;
  for ( i = 0 ; i < grid->ngem ; i++ ) 
    gridRemoveCell( grid, grid->gem[i] );
  return grid;
}

int gridCellDegree(Grid *grid, int id)
{
  return adjDegree(grid->cellAdj, id);
}

int gridAddCell(Grid *grid, int n0, int n1, int n2, int n3)
{
  int cellId;
  if ( grid->blankc2n == EMPTY ) {
    int cell, currentSize, chunkSize;
    chunkSize = 5000;
    currentSize = grid->maxcell;
    grid->maxcell += chunkSize;
    grid->c2n = (int *)realloc( grid->c2n, 4 * grid->maxcell * sizeof(int) );
    for (cell=currentSize;cell < grid->maxcell; cell++ ) {
      grid->c2n[0+4*cell] = EMPTY; 
      grid->c2n[1+4*cell] = cell+1; 
    }
    grid->c2n[1+4*(grid->maxcell-1)] = EMPTY; 
    grid->blankc2n = currentSize;
  }
  cellId = grid->blankc2n;
  grid->blankc2n = grid->c2n[1+4*cellId];
  grid->ncell++;
  
  grid->c2n[0+4*cellId] = n0;
  grid->c2n[1+4*cellId] = n1;
  grid->c2n[2+4*cellId] = n2;
  grid->c2n[3+4*cellId] = n3;
  
  if ( NULL == adjRegister( grid->cellAdj, n0, cellId ) ) return EMPTY;
  if ( NULL == adjRegister( grid->cellAdj, n1, cellId ) ) return EMPTY;
  if ( NULL == adjRegister( grid->cellAdj, n2, cellId ) ) return EMPTY;
  if ( NULL == adjRegister( grid->cellAdj, n3, cellId ) ) return EMPTY;
  
  return cellId;
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

Grid *gridReconnectAllCell(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int cell, i, node;

  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;
  if (newNode == oldNode) return grid;

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
  if (newNode == oldNode) return grid;
  
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

bool gridCellValid(Grid *grid, int cellId )
{
  if ( cellId < 0 || cellId >= grid->maxcell ) return FALSE;
  return (EMPTY != grid->c2n[4*cellId]);
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
  return (EMPTY != gridFindOtherCellWith3Nodes(grid, n0, n1, n2, EMPTY ) );
}

int gridFindOtherCellWith3Nodes(Grid *grid, int n0, int n1, int n2,
				int currentCell )
{
  AdjIterator it;
  int cell, i, j, nodes[4];

  for ( it = adjFirst(grid->cellAdj,n0); adjValid(it); it = adjNext(it) ) {
    cell = adjItem(it);
    if (cell != currentCell) {
      gridCell( grid, cell, nodes );
      if ( ( n1==nodes[0] || n1==nodes[1] || n1==nodes[2] || n1==nodes[3] ) &&
	   ( n2==nodes[0] || n2==nodes[1] || n2==nodes[2] || n2==nodes[3] ) )
	return cell;
    }
  }

  return EMPTY;
}

int gridFindCellWithFace(Grid *grid, int face ){
  int nodes[3], faceId;

  if (grid != gridFace(grid,face,nodes,&faceId)) return EMPTY;

  return gridFindOtherCellWith3Nodes(grid, nodes[0],nodes[1],nodes[2], EMPTY);
}

Grid *gridCheckCellConnections(Grid *grid){
  int cell, nodes[4];
  
  for (cell = 0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if (grid == gridCell(grid,cell,nodes)) {
      
    }
  }

  return grid;
}


Grid *gridDeleteThawedCells(Grid *grid){
  int cell, maxcell, nodes[4];

  maxcell = gridMaxCell(grid);
  for( cell=0 ; cell < maxcell ; cell++ ) {
    if (grid == gridCell(grid,cell,nodes)) {
      if ( !gridNodeFrozen(grid,nodes[0]) || 
	   !gridNodeFrozen(grid,nodes[1]) || 
	   !gridNodeFrozen(grid,nodes[2]) ||
	   !gridNodeFrozen(grid,nodes[3])  ) {
	gridRemoveCell(grid,cell);
      }
    }
  }

  return grid;
}

int gridAddFace(Grid *grid, int n0, int n1, int n2, int faceId )
{
 return gridAddFaceUV(grid, 
		      n0, DBL_MAX, DBL_MAX,
		      n1, DBL_MAX, DBL_MAX,
		      n2, DBL_MAX, DBL_MAX,
		      faceId );
}

int gridAddFaceUV(Grid *grid, 
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId )
{
  int face, i, origSize, chunkSize;
  if ( grid->blankf2n == EMPTY ) {
    chunkSize = 5000;
    origSize = grid->maxface;
    grid->maxface += chunkSize;
    grid->f2n    = realloc(grid->f2n,    3 * grid->maxface * sizeof(int));
    grid->faceU  = realloc(grid->faceU,  3 * grid->maxface * sizeof(double));
    grid->faceV  = realloc(grid->faceV,  3 * grid->maxface * sizeof(double));
    grid->faceId = realloc(grid->faceId, 1 * grid->maxface * sizeof(int));
    for (i=origSize;i < grid->maxface; i++ ) {
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
    grid->blankf2n = origSize;
  }
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

  if ( NULL == adjRegister( grid->faceAdj, n0, face ) ) return EMPTY;
  if ( NULL == adjRegister( grid->faceAdj, n1, face ) ) return EMPTY;
  if ( NULL == adjRegister( grid->faceAdj, n2, face ) ) return EMPTY;

 return face;
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

Grid *gridReconnectAllFace(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int face, i, node;
  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;
  if (newNode == oldNode) return grid;

  it = adjFirst(grid->faceAdj,oldNode);
  while (adjValid(it)){
    face = adjItem(it);
    for (i=0;i<3;i++){
      node = grid->f2n[i+3*face];
      if (oldNode == node) {
	grid->f2n[i+3*face]=newNode;
	adjRemove( grid->faceAdj, oldNode, face);
	adjRegister( grid->faceAdj, newNode, face);
      }
    }
    it = adjFirst(grid->faceAdj,oldNode);
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

Grid *gridDeleteThawedFaces(Grid *grid, int faceId ){
  int face, maxface, nodes[3], id;

  if (faceId < 1) return NULL;

  maxface = gridMaxFace(grid);
  for( face=0 ; face < maxface ; face++ ) {
    if (grid == gridFace(grid,face,nodes,&id)) {
      if ( id == faceId && 
	   ( !gridNodeFrozen(grid,nodes[0]) || 
	     !gridNodeFrozen(grid,nodes[1]) || 
	     !gridNodeFrozen(grid,nodes[2]) ) ){
	gridRemoveFace(grid,face);
      }
    }
  }

  return grid;
}

int gridNThawedFaces(Grid *grid, int faceId ){
  int face, maxface, nface, nodes[3], id;

  if (faceId < 1) return EMPTY;

  nface = 0;
  maxface = gridMaxFace(grid);
  for( face=0 ; face < maxface ; face++ ) {
    if (grid == gridFace(grid,face,nodes,&id)) {
      if ( id == faceId && 
	   ( !gridNodeFrozen(grid,nodes[0]) || 
	     !gridNodeFrozen(grid,nodes[1]) || 
	     !gridNodeFrozen(grid,nodes[2]) ) ){
	nface++;
      }
    }
  }

  return nface;
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

int gridAddEdge(Grid *grid, int n0, int n1, 
		int edgeId, double t0, double t1 )
{
  int edge;
  if ( grid->blanke2n == EMPTY ) return EMPTY;
  edge = grid->blanke2n;
  grid->blanke2n = grid->e2n[1+2*edge];
  grid->nedge++;

  grid->e2n[0+2*edge] = n0;
  grid->e2n[1+2*edge] = n1;
  grid->edgeId[edge]  = edgeId;
  grid->edgeT[0+2*edge] = t0;
  grid->edgeT[1+2*edge] = t1;

  if ( NULL == adjRegister( grid->edgeAdj, n0, edge ) ) return EMPTY;
  if ( NULL == adjRegister( grid->edgeAdj, n1, edge ) ) return EMPTY;

  return edge;
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

Grid *gridReconnectAllEdge(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int edge, i, node;
  if (oldNode < 0 || oldNode >= grid->maxnode ) return NULL;
  if (newNode < 0 || newNode >= grid->maxnode ) return NULL;

  it = adjFirst(grid->edgeAdj,oldNode);
  while (adjValid(it)){
    edge = adjItem(it);
    for (i=0;i<2;i++){
      node = grid->e2n[i+2*edge];
      if (oldNode == node) {
	grid->e2n[i+2*edge]=newNode;
	adjRemove( grid->edgeAdj, oldNode, edge);
	adjRegister( grid->edgeAdj, newNode, edge);
      }
    }
    it = adjFirst(grid->edgeAdj,oldNode);
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

Grid *gridDeleteThawedEdgeSegments(Grid *grid, int edgeId ){
  int edge, maxedge, nodes[2], id;

  if (edgeId < 1) return NULL;

  maxedge = gridMaxEdge(grid);
  for( edge=0 ; edge < maxedge ; edge++ ) {
    if (grid == gridEdge(grid,edge,nodes,&id)) {
      if ( id == edgeId && 
	   ( !gridNodeFrozen(grid,nodes[0]) || 
	     !gridNodeFrozen(grid,nodes[1]) ) ){
	gridRemoveEdge(grid,edge);
      }
    }
  }

  return grid;
}

int gridNThawedEdgeSegments(Grid *grid, int edgeId ){
  int edge, maxedge, nodes[2], id;
  int nSegments;

  if (edgeId < 1) return EMPTY;

  nSegments = 0;
  maxedge = gridMaxEdge(grid);
  for( edge=0 ; edge < maxedge ; edge++ ) {
    if (grid == gridEdge(grid,edge,nodes,&id)) {
      if ( id == edgeId && 
	   ( !gridNodeFrozen(grid,nodes[0]) || 
	     !gridNodeFrozen(grid,nodes[1]) ) ){
	nSegments++;
      }
    }
  }

  return nSegments;
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
	  if ( n1 == startNode ) found = FALSE; /* periodic edge check */
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
	  if ( n1 == startNode ) found = FALSE; /* periodic edge check */
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
	  if ( n1 == startNode ) found = FALSE; /* periodic edge check */
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

  for ( iequ = 0; iequ < grid->nequ; iequ++ )
    grid->equ[iequ+grid->nequ]=grid->equ[iequ];

  return grid;
}

int gridNEqu(Grid *grid)
{
  return grid->nequ;
}

int gridEqu(Grid *grid, int index)
{
  return grid->equ[index];
}

bool gridContinuousEquator(Grid *grid)
{
  return (gridNEqu(grid) == gridNGem(grid));
}

Grid *gridCycleEquator( Grid *grid )
{
  int i;

  for ( i = grid->nequ ; i > 0 ; i-- )
    grid->equ[i] = grid->equ[i-1];

  grid->equ[0] = grid->equ[grid->nequ];

  return grid;
}

int gridAddNode(Grid *grid, double x, double y, double z )
{
  int node, i, origSize, chunkSize;

  if (EMPTY == grid->blanknode) {
    chunkSize = 5000;
    origSize = grid->maxnode;
    grid->maxnode += chunkSize;
    
    grid->xyz = realloc(grid->xyz, 3 * grid->maxnode * sizeof(double));
    for (i=origSize;i < grid->maxnode; i++ ) {
      grid->xyz[0+3*i] = DBL_MAX;
      grid->xyz[1+3*i] = (double)(i+1);
    }
    grid->xyz[1+3*(grid->maxnode-1)] = (double)(EMPTY);
    grid->blanknode = origSize;

    grid->map = realloc(grid->map, grid->maxnode * 6 * sizeof(double));
    grid->frozen = realloc(grid->frozen,grid->maxnode * sizeof(bool));

    adjRealloc(grid->cellAdj,grid->maxnode);
    adjRealloc(grid->faceAdj,grid->maxnode);
    adjRealloc(grid->edgeAdj,grid->maxnode);

    if (grid->prismDeg != NULL ) {
      grid->prismDeg = realloc(grid->prismDeg,grid->maxnode * sizeof(int));
      for (i=origSize;i < grid->maxnode; i++ ) grid->prismDeg[i] = 0;
    }
  }

  node = grid->blanknode;
  grid->blanknode = (int)grid->xyz[1+3*node];
  grid->nnode++;

  grid->xyz[0+3*node] = x;
  grid->xyz[1+3*node] = y;
  grid->xyz[2+3*node] = z;
  grid->frozen[node] = FALSE;
  grid->map[0+6*node] = 1.0;
  grid->map[1+6*node] = 0.0;
  grid->map[2+6*node] = 0.0;
  grid->map[3+6*node] = 1.0;
  grid->map[4+6*node] = 0.0;
  grid->map[5+6*node] = 1.0;

  return node;
}

Grid *gridRemoveNode(Grid *grid, int node )
{
  if (node < 0 || node>grid->maxnode) return NULL;
  if (DBL_MAX == grid->xyz[0+3*node]) return NULL;
  grid->nnode--;
  grid->xyz[0+3*node] = DBL_MAX;
  grid->xyz[1+3*node] = (double)grid->blanknode;
  grid->blanknode = node;
  return grid;
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

Grid *gridDeleteNodesNotUsed(Grid *grid){
  int node, maxnode;
  bool prism;

  maxnode = gridMaxNode(grid);
  prism = FALSE;
  for( node=0 ; node < maxnode ; node++ ) {
    if ( grid->prismDeg != NULL ) prism = (grid->prismDeg[node] > 0);
    if ( gridValidNode(grid,node) && 
	 0 == gridCellDegree(grid, node) &&
	 !gridGeometryFace(grid, node) &&
	 !gridGeometryEdge(grid, node) &&
	 !prism )
      gridRemoveNode(grid,node);
  }

  return grid;
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
  int i;
  grid->nGeomEdge = nGeomEdge;
  if ( NULL != grid->geomEdge) free(grid->geomEdge);
  grid->geomEdge = malloc(2*nGeomEdge*sizeof(int));
  for (i=0;i<2*nGeomEdge;i++) grid->geomEdge[i]=EMPTY;
  return grid;
}

int gridNGeomFace(Grid *grid)
{
  return grid->nGeomFace;
}

Grid *gridSetNGeomFace(Grid *grid, int nGeomFace)
{
  grid->nGeomFace = nGeomFace;
  return grid;
}

Grid *gridAddGeomEdge(Grid *grid, int edge, int n0, int n1 )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return NULL;
  grid->geomEdge[0+(edge-1)*2] = n0;
  grid->geomEdge[1+(edge-1)*2] = n1;
  return grid;
}

int gridGeomEdgeStart( Grid *grid, int edge )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return EMPTY;
  if ( NULL == grid->geomEdge) return EMPTY;
  return grid->geomEdge[0+2*(edge-1)];
}

int gridGeomEdgeEnd( Grid *grid, int edge )
{
  if ( edge<1 || edge>grid->nGeomEdge ) return EMPTY;
  if ( NULL == grid->geomEdge) return EMPTY;
  return grid->geomEdge[1+2*(edge-1)];
}

int gridGeomEdgeSize( Grid *grid, int edge )
{
  int startNode;
  startNode = gridGeomEdgeStart( grid, edge );
  if ( startNode == EMPTY ) return EMPTY;
  return gridGeomCurveSize( grid, edge, startNode );
}

Grid *gridGeomEdge( Grid *grid, int edge, int *curve )
{
  int startNode;
  startNode = gridGeomEdgeStart( grid, edge );
  if ( startNode == EMPTY ) return NULL;
  return gridGeomCurve( grid, edge, startNode, curve );
}

int gridFrozenEdgeEndPoint( Grid *grid, int edgeId, int startNode )
{
  AdjIterator it;
  int node, lastnode, edge, n1, degree;
  bool found;

  degree =0;
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
	degree++;
	if ( node == grid->e2n[0+2*edge] ) {
	  n1 = grid->e2n[1+2*edge];
	}else{
	  n1 = grid->e2n[0+2*edge];	  
	}
	if ( n1 != lastnode &&
	     gridNodeFrozen(grid, n1) &&
	     gridNodeFrozen(grid, node)) { 
	  found = TRUE;
	  lastnode = node;
	  node = n1;
	}
      }
    }
  }

  if (degree==0)return EMPTY;
  return node;
}

bool gridGeometryNode(Grid *grid, int node)
{
  if (node < 0 ) return FALSE;
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

Grid *gridAddPrism(Grid *grid, int n0, int n1, int n2, int n3, int n4, int n5)
{
  int i;

  if (grid->nprism >= grid->maxprism) {
    grid->maxprism += 5000;
    if (grid->prism == NULL) {
      grid->prism = malloc(grid->maxprism*sizeof(Prism));
      grid->prismDeg = malloc(grid->maxnode*sizeof(int));
      for(i=0;i<grid->maxnode;i++) grid->prismDeg[i]=0;
    }else{
      grid->prism = realloc(grid->prism,grid->maxprism*sizeof(Prism));
    }
  }

  grid->prism[grid->nprism].nodes[0] = n0;
  grid->prism[grid->nprism].nodes[1] = n1;
  grid->prism[grid->nprism].nodes[2] = n2;
  grid->prism[grid->nprism].nodes[3] = n3;
  grid->prism[grid->nprism].nodes[4] = n4;
  grid->prism[grid->nprism].nodes[5] = n5;

  for(i=0;i<6;i++) grid->prismDeg[grid->prism[grid->nprism].nodes[i]]++;

  grid->nprism++;

  return grid;
}

Grid *gridPrism(Grid *grid, int prismIndex, int *nodes)
{
  int i;
  if (prismIndex<0 || prismIndex >= gridNPrism(grid) ) return NULL; 

  for (i=0;i<6;i++){
    nodes[i]=grid->prism[prismIndex].nodes[i];
  }

  return grid;
}


Grid *gridAddPyramid(Grid *grid, int n0, int n1, int n2, int n3, int n4)
{

  if (grid->npyramid >= grid->maxpyramid) {
    grid->maxpyramid += 5000;
    if (grid->pyramid == NULL) {
      grid->pyramid = malloc(grid->maxpyramid*sizeof(Pyramid));
    }else{
      grid->pyramid = realloc(grid->pyramid,grid->maxpyramid*sizeof(Pyramid));
    }
  }

  grid->pyramid[grid->npyramid].nodes[0] = n0;
  grid->pyramid[grid->npyramid].nodes[1] = n1;
  grid->pyramid[grid->npyramid].nodes[2] = n2;
  grid->pyramid[grid->npyramid].nodes[3] = n3;
  grid->pyramid[grid->npyramid].nodes[4] = n4;

  grid->npyramid++;

  return grid;
}

Grid *gridPyramid(Grid *grid, int pyramidIndex, int *nodes)
{
  int i;
  if (pyramidIndex<0 || pyramidIndex >= gridNPyramid(grid) ) return NULL; 

  for (i=0;i<5;i++){
    nodes[i]=grid->pyramid[pyramidIndex].nodes[i];
  }

  return grid;
}

Grid *gridAddQuad(Grid *grid, int n0, int n1, int n2, int n3, int faceId )
{

  if (grid->nquad >= grid->maxquad) {
    grid->maxquad += 5000;
    if (grid->quad == NULL) {
      grid->quad = malloc(grid->maxquad*sizeof(Quad));
    }else{
      grid->quad = realloc(grid->quad,grid->maxquad*sizeof(Quad));
    }
  }

  grid->quad[grid->nquad].nodes[0] = n0;
  grid->quad[grid->nquad].nodes[1] = n1;
  grid->quad[grid->nquad].nodes[2] = n2;
  grid->quad[grid->nquad].nodes[3] = n3;
  grid->quad[grid->nquad].faceId = faceId;

  grid->nquad++;

  return grid;
}

Grid *gridQuad(Grid *grid, int quadIndex, int *nodes, int *faceId )
{
  int i;
  if (quadIndex<0 || quadIndex >= gridNQuad(grid) ) return NULL; 

  for (i=0;i<4;i++){
    nodes[i]=grid->quad[quadIndex].nodes[i];
  }
  *faceId = grid->quad[quadIndex].faceId;

  return grid;
}

Grid *gridMap(Grid *grid, int node, double *map)
{
  if ( !gridValidNode(grid, node) ) return NULL;

  map[0] = grid->map[0+6*node];
  map[1] = grid->map[1+6*node];
  map[2] = grid->map[2+6*node];
  map[3] = grid->map[3+6*node];
  map[4] = grid->map[4+6*node];
  map[5] = grid->map[5+6*node];
  
  return grid;
}

Grid *gridSetMap(Grid *grid, int node,
		 double m11, double m12, double m13,
		             double m22, double m23,
		                         double m33)
{
  if ( !gridValidNode(grid, node) ) return NULL;

  grid->map[0+6*node] = m11;
  grid->map[1+6*node] = m12;
  grid->map[2+6*node] = m13;
  grid->map[3+6*node] = m22;
  grid->map[4+6*node] = m23;
  grid->map[5+6*node] = m33;

  return grid;
}

int gridStoredARDegree( Grid *grid )
{
  return grid->degAR;
}

Grid *gridClearStoredAR( Grid *grid )
{
  grid->degAR = 0;
  return grid;
}

Grid *gridAddStoredAR( Grid *grid, double AR, double *dARdX )
{
  if (grid->degAR == MAXDEG) return NULL;

  grid->AR[grid->degAR] = AR;
  grid->dARdX[0+3*grid->degAR] = dARdX[0];
  grid->dARdX[1+3*grid->degAR] = dARdX[1];
  grid->dARdX[2+3*grid->degAR] = dARdX[2];

  grid->degAR++;
  return grid;
}

double gridStoredAR( Grid *grid, int index )
{
  if ( index < 0 || index >= gridStoredARDegree( grid ) ) return DBL_MAX;
  return grid->AR[index];
}

Grid *gridStoredARDerivative( Grid *grid, int index, double *dARdX )
{
  if ( index < 0 || index >= gridStoredARDegree( grid ) ) return NULL;
  dARdX[0] = grid->dARdX[0+3*index];
  dARdX[1] = grid->dARdX[1+3*index];
  dARdX[2] = grid->dARdX[2+3*index];
  
  return grid;
}

Grid *gridFreezeLinesNodes(Grid *grid)
{
  int line, index, node;
  Lines *lines;

  lines = gridLines(grid);

  for(line=0;line<linesNumber(lines);line++){
    index = 0;
    node = linesNode(lines,line,index);
    while (node>=0) {
      gridFreezeNode(grid,node);
      index++;
      node = linesNode(lines,line,index);
    }
  }

  return grid;
}

Grid *gridReportLinesLocation(Grid *grid)
{
  int line, index, node;
  double xyz[3];
  Lines *lines;

  lines = gridLines(grid);

  for(line=0;line<linesNumber(lines);line++){
    index = 0;
    node = linesNode(lines,line,index);
    while (node>=0) {
      gridNodeXYZ(grid,node,xyz);
      printf("line %5d index %5d node %5d x %10.8f y %10.8f z %10.8f\n",
	     line,index,node,xyz[0],xyz[1],xyz[2]);
      index++;
      node = linesNode(lines,line,index);
    }
  }
  
  return grid;
}

Grid *gridCopyAboutY0(Grid *grid)
{
  int node, orignode, origface, origcell;
  int *o2n;
  double xyz[3];
  int nodes[4];
  int face, i, faceid;
  int cell;

  printf("gridCopyAboutY0: pack\n");
  if (NULL == gridPack(grid)) {
    printf("gridCopyAboutY0: gridPack failed.\n");
    return NULL;
  }

  orignode = gridNNode(grid);

  o2n = malloc(sizeof(int) * orignode); 

  printf("gridCopyAboutY0: copy nodes, y = -y\n");

  for ( node = 0 ; node < orignode ; node++){
    gridNodeXYZ(grid,node,xyz);
    o2n[node] = gridAddNode(grid,xyz[0],-xyz[1],xyz[2]);
  }

  printf("gridCopyAboutY0: remove duplicate nodes.\n");

  origface = gridNFace(grid);
  for ( face = 0 ; face < origface ; face++ ){
    gridFace(grid,face,nodes,&faceid);
    if ( faceid == 5 || faceid == 6 ) {
      for (i=0;i<3;i++){
	node = nodes[i];
	if (o2n[node] >= orignode){
	  gridRemoveNode(grid,o2n[node]);
	  o2n[node] = node;
	}
      }
    }
  }

  printf("gridCopyAboutY0: copy faces, swap node 0 and 1\n");

  for ( face = 0 ; face < origface ; face++ ){
    gridFace(grid,face,nodes,&faceid);
    gridAddFace(grid,o2n[nodes[1]],o2n[nodes[0]],o2n[nodes[2]],faceid);
  }


  printf("gridCopyAboutY0: copy cells, swap node 0 and 1\n");

  origcell = gridNCell(grid);
  for ( cell = 0 ; cell < origcell ; cell++ ){
    gridCell(grid,cell,nodes);
    gridAddCell(grid,o2n[nodes[1]],o2n[nodes[0]],o2n[nodes[2]],o2n[nodes[3]]);
  }

  printf("gridCopyAboutY0: remove sym face\n");

  gridThawAll(grid);
  gridDeleteThawedFaces(grid,5);
  gridDeleteThawedFaces(grid,6);

  free(o2n);
    
  return grid;
}
