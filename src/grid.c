
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
#ifdef __APPLE__       /* Not needed on Mac OS X */
#include <float.h>
#else
#include <malloc.h>
#include <values.h>
#endif
#include "sort.h"
#include "gridmath.h"
#include "grid.h"

Grid* gridCreate(int maxnode, int maxcell, int maxface, int maxedge)
{
  return gridImport(maxnode, 0, 
		    maxface, 0, 
		    maxcell, 0,
		    maxedge, 
		    NULL, NULL, NULL, NULL);
}

Grid *gridImport(int maxnode, int nnode, 
		 int maxface, int nface, 
		 int maxcell, int ncell,
		 int maxedge, 
		 double *xyz, int *f2n, int *faceId, int *c2n )
{
  int i;
  Grid *grid;

  grid = (Grid *)malloc(sizeof(Grid));

  grid->parent = NULL;
  grid->child  = NULL;

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
  grid->partId=0;
  grid->globalNNode=0;
  grid->globalNCell=0;
  grid->nGeomNode = 0;
  grid->nGeomEdge = 0;
  grid->nGeomFace = 0;
  grid->geomEdge = NULL;

  if (NULL == xyz ){
    grid->xyz = (double *)malloc(3 * grid->maxnode * sizeof(double));
  }else{
    grid->xyz = xyz;
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

  grid->map = (double *)malloc(grid->maxnode * 6 * sizeof(double));
  for (i=0;i < grid->maxnode; i++ ) {
    grid->map[0+6*i] = 1.0;
    grid->map[1+6*i] = 0.0;
    grid->map[2+6*i] = 0.0;
    grid->map[3+6*i] = 1.0;
    grid->map[4+6*i] = 0.0;
    grid->map[5+6*i] = 1.0;
  }

  grid->frozen = (GridBool *)malloc(grid->maxnode * sizeof(GridBool));
  for (i=0;i < grid->maxnode; i++ ) grid->frozen[i] = FALSE;

  grid->geomNode = (int *)malloc(grid->maxnode * sizeof(int));
  for (i=0;i < grid->maxnode; i++ ) grid->geomNode[i] = EMPTY;

  grid->child_reference  = NULL;
  grid->nodeGlobal  = NULL;
  grid->part = NULL;
  grid->sortedGlobal = NULL;
  grid->sortedLocal = NULL;
  grid->maxUnusedNodeGlobal = 0;
  grid->nUnusedNodeGlobal = 0;
  grid->unusedNodeGlobal  = NULL;

  grid->aux = NULL;
  grid->naux = 0;

  /* cells */
  if ( NULL == c2n ) {
    grid->c2n = (int *)malloc(4 * grid->maxcell * sizeof(int));
  } else {
    grid->c2n = c2n;
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
  grid->cellAdj = adjCreate(grid->maxnode,grid->maxcell*4,5000*4);

  for ( i=0 ; i < grid->ncell ; i++ ) {
    adjRegister(grid->cellAdj,grid->c2n[0+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[1+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[2+4*i],i);
    adjRegister(grid->cellAdj,grid->c2n[3+4*i],i);
  }

  grid->cellGlobal  = NULL;
  grid->maxUnusedCellGlobal = 0;
  grid->nUnusedCellGlobal = 0;
  grid->unusedCellGlobal  = NULL;

  grid->constrain_surface_node = FALSE;

  if (NULL == f2n) {
    grid->f2n    = (int *)malloc(3 * grid->maxface * sizeof(int));
  }else{
    grid->f2n    = f2n;
  }

  if (NULL == f2n) { /* this should be (NULL == faceId) ? */
    grid->faceId = (int *)malloc(1 * grid->maxface * sizeof(int));
  }else{
    grid->faceId = faceId;
  }

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
  grid->faceU  = (double *)malloc(3 * grid->maxface * sizeof(double));
  grid->faceV  = (double *)malloc(3 * grid->maxface * sizeof(double));
  for (i=0;i < grid->maxface; i++ ) {
    grid->faceU[0+3*i] = DBL_MAX;
    grid->faceU[1+3*i] = DBL_MAX;
    grid->faceU[2+3*i] = DBL_MAX;
    grid->faceV[0+3*i] = DBL_MAX;
    grid->faceV[1+3*i] = DBL_MAX;
    grid->faceV[2+3*i] = DBL_MAX;
  }

  grid->faceAdj = adjCreate(grid->maxnode,grid->maxface*3,5000*3);

  /* addface */

  for ( i=0 ; i < grid->nface ; i++ ) {
    adjRegister(grid->faceAdj,grid->f2n[0+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[1+3*i],i);
    adjRegister(grid->faceAdj,grid->f2n[2+3*i],i);
  }

  /*  edge */
  grid->e2n    = (int *)malloc(2 * grid->maxedge * sizeof(int));
  grid->edgeId = (int *)malloc(1 * grid->maxedge * sizeof(int));
  grid->edgeT  = (double *)malloc(2 * grid->maxedge * sizeof(double));
  for (i=0;i < grid->maxedge; i++ ) {
    grid->e2n[0+2*i] = EMPTY; 
    grid->e2n[1+2*i] = i+1; 
    grid->edgeId[i] = EMPTY; 
    grid->edgeT[0+2*i] = DBL_MAX; 
    grid->edgeT[1+2*i] = DBL_MAX; 
  }
  grid->e2n[1+2*(grid->maxedge-1)] = EMPTY; 
  grid->blanke2n = 0;

  grid->edgeAdj = adjCreate(grid->maxnode,grid->maxedge*2,5000*2);

  grid->nequ = 0;
  grid->ngem = 0;

  grid->nconn = 0;
  grid->cell2conn = NULL;
  grid->conn2node = NULL;

  grid->costDegree = 0;

  grid->costFunction = gridCOST_FCN_MEAN_RATIO;
  grid->costConstraint = gridCOST_CNST_VOLUME;
  grid->min_allowed_insert_cost = 0.01;
  grid->min_allowed_surface_smooth_cost = 0.01;
  grid->min_allowed_swap_cost = -0.5;
  grid->min_allowed_swap_cost_improvement = 1.0e-12;

  grid->tecplotGeomFile = NULL;
  grid->tecplotScalarFile = NULL;

  grid->packFunc = NULL;
  grid->packData = NULL;

  grid->renumberFunc = NULL;
  grid->renumberData = NULL;

  grid->reallocFunc = NULL;
  grid->reallocData = NULL;

  grid->freeNotificationFunc = NULL;
  grid->freeNotificationData = NULL;

  grid->lines = linesCreate();

  grid->phase = gridALL_PHASE;

  grid->model = -1;

  return  grid;
}

Grid *gridDup( Grid *grid )
{
  Grid *child;
  int node;
  double xyz[3], map[6];
  int cell;
  int nodes[4];
  int face, faceId;
  double uvs[3][2];
  int edge, edgeId;
  double ts[2];

  if (grid != gridSortNodeGridEx(grid)) {
    printf("%s: %d: gridDup: gridSortNodeGridEx failed.\n",
	   __FILE__, __LINE__ );
    return NULL;
  }

  child = gridCreate( gridNNode(grid), gridNCell(grid),
		      gridNFace(grid), gridNEdge(grid) );

  for ( node = 0 ; node < gridNNode(grid) ; node++ ) {
    gridNodeXYZ( grid,  node, xyz );
    gridAddNode( child, xyz[0], xyz[1], xyz[2] );
    gridMap( grid, node, map );
    gridSetMap( child, node, map[0], map[1], map[2], map[3], map[4], map[5] );
  }

  for ( cell = 0 ; cell < gridNCell(grid) ; cell++ ) {
    gridCell( grid, cell, nodes );
    gridAddCell( child, nodes[0], nodes[1], nodes[2], nodes[3] );
  }

  for ( face = 0 ; face < gridNFace(grid) ; face++ ) {
    gridFace( grid, face, nodes, &faceId );
    gridNodeUV( grid, nodes[0], faceId, uvs[0] );
    gridNodeUV( grid, nodes[1], faceId, uvs[1] );
    gridNodeUV( grid, nodes[2], faceId, uvs[2] );
    gridAddFaceUV( child, 
		   nodes[0], uvs[0][0],  uvs[0][1], 
		   nodes[1], uvs[1][0],  uvs[1][1], 
		   nodes[2], uvs[2][0],  uvs[2][1], 
		   faceId);
  }

  for ( edge = 0 ; edge < gridNEdge(grid) ; edge++ ) {
    gridEdge( grid, edge, nodes, &edgeId );
    gridNodeT( grid, nodes[0], edgeId, &(ts[0]) );
    gridNodeT( grid, nodes[1], edgeId, &(ts[1]) );
    gridAddEdge( child, nodes[0], nodes[1], edgeId, ts[0], ts[1] );
  }

  gridSetNGeomNode( child, gridNGeomNode( grid ) );
  gridSetNGeomEdge( child, gridNGeomEdge( grid ) );
  gridSetNGeomFace( child, gridNGeomFace( grid ) );

  for ( edgeId = 1 ; edgeId <= gridNGeomEdge( grid ) ; edgeId++ ) {
    gridAddGeomEdge( child, edgeId, 
		     gridGeomEdgeStart( grid, edgeId ), 
		     gridGeomEdgeEnd(   grid, edgeId ) );
  }

  return child;
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

Grid *gridImportNGP( char *filename )
{
  FILE *file;
  int nnode, maxface, ncell;
  int node, prefix;
  double xyz[3];
  int nodes[3];
  int face;
  int lc,rc, bc;
  int *c2n;
  int cell;
  Grid *grid;
 
  file = fopen(filename,"r");\
  if (NULL == file) return NULL;

  if ( 3 != fscanf(file,"%d %d %d",&nnode,&ncell,&maxface) )
    {
      printf("ERROR: gridImportNGP: %s: %d: header read\n",
	     __FILE__, __LINE__ );
      return NULL;
    }

  printf("nnode %d ncell %d maxface %d\n",nnode,ncell,maxface);

  grid = gridCreate(nnode, ncell, maxface, 0);
  if (NULL == grid) return NULL;

  for ( node = 0; node < nnode ; node++ )
    {
      if (4 != fscanf(file,"%d %lf %lf %lf",
		      &prefix,&(xyz[0]),&(xyz[1]),&(xyz[2]) ) )
	{
	  printf("ERROR: gridImportNGP: %s: %d: node %d read\n",
		 __FILE__, __LINE__, (node+1) );
	  gridFree(grid); return NULL;
	}
      if ( (node+1) != prefix )
	{
	  printf("ERROR: gridImportNGP: %s: %d: node %d prefix %d mismatch\n",
		 __FILE__, __LINE__, (node+1), prefix );
	  gridFree(grid); return NULL;
	}
      if ( node != gridAddNode( grid, xyz[0], xyz[1], xyz[2] ) )
	{
	  printf("ERROR: gridImportNGP: %s: %d: gridAddNode %d failed\n",
		 __FILE__, __LINE__, node );
	  gridFree(grid); return NULL;
	}
    }

  c2n = (int *)malloc( 4*ncell*sizeof(int) );
  for ( node = 0 ; node < 4*ncell ; node++ ) c2n[node] = EMPTY;

  for ( face = 0; face < maxface ; face++ )
    {
      if (4 != fscanf(file,"%d %d %d %d",
		      &prefix,&(nodes[0]),&(nodes[1]),&(nodes[2]) ) )
	{
	  printf("ERROR: gridImportNGP: %s: %d: face %d read line 1\n",
		 __FILE__, __LINE__, face );
	  gridFree(grid); return NULL;
	}
      if ( (face+1) != prefix )
	{
	  printf("ERROR: gridImportNGP: %s: %d: face %d prefix %d mismatch\n",
		 __FILE__, __LINE__, (face+1), prefix );
	  gridFree(grid); return NULL;
	}
      nodes[0]--;
      nodes[1]--;
      nodes[2]--;
      if (4 != fscanf(file,"%d %d %d %d",
		      &prefix, &lc, &rc, &bc ) )
	{
	  printf("ERROR: gridImportNGP: %s: %d: face %d read line 2\n",
		 __FILE__, __LINE__, face );
	  gridFree(grid); return NULL;
	}
      if ( (face+1) != prefix )
	{
	  printf("ERROR: gridImportNGP: %s: %d: face %d prefix %d mismatch\n",
		 __FILE__, __LINE__, (face+1), prefix );
	  gridFree(grid); return NULL;
	}
      if ( lc > 0 ) 
	{
	  cell = lc-1;
	  if ( EMPTY == c2n[0+4*cell] )
	    {
	      for (node=0;node<3;node++) c2n[node+4*cell]=nodes[node];
	    }
	  else
	    {
	      for (node=0;node<3;node++)
		{
		  if ( nodes[node] != c2n[0+4*cell] &&
		       nodes[node] != c2n[1+4*cell] &&
		       nodes[node] != c2n[2+4*cell] ) 
		    c2n[3+4*cell] = nodes[node];
		}
	    }
	}
      else
	{
	  gridAddFace( grid, nodes[0], nodes[1], nodes[2], ABS(lc) );
	}
      if ( rc > 0 ) 
	{
	  cell = rc-1;
	  if ( EMPTY == c2n[0+4*cell] )
	    {
	      for (node=0;node<3;node++) c2n[node+4*cell]=nodes[node];
	    }
	  else
	    {
	      for (node=0;node<3;node++)
		{
		  if ( nodes[node] != c2n[0+4*cell] &&
		       nodes[node] != c2n[1+4*cell] &&
		       nodes[node] != c2n[2+4*cell] ) 
		    c2n[3+4*cell] = nodes[node];
		}
	    }	  
	}
      else
	{
	  gridAddFace( grid, nodes[0], nodes[1], nodes[2], ABS(rc) );
	}

    }

  for (cell=0;cell<ncell;cell++)
    {
      gridAddCell(grid, 
		  c2n[0+4*cell],c2n[1+4*cell],c2n[2+4*cell],c2n[3+4*cell]);
    }

  free(c2n);

  return grid;
}

Grid *gridExportNGP( Grid *grid, char *filename )
{
  FILE *file;
  int *c2f;
  int face;
  int nface;
  int node;
  int cell, cell_face;
  int face_nodes[3];
  int cell_nodes[4];
  int other_cell;
  int other_nodes[4];
  int other_face;
  int other_face_nodes[3];
 
  int cellfacenode[4][3] = {{1, 3, 2}, {0,2,3}, {0,3,1}, {0,1,2}};

  if ( FALSE )
    {
      if (NULL == gridPack(grid)) {
	printf("gridExportNGP: gridPack failed.\n");
	return NULL;
      }
    }
  else
    {
      printf("gridExportNGP grid Pack skipped\n");
    }

  c2f = (int *)malloc( 4 * gridNCell(grid) * sizeof(int) );
  for ( face = 0 ; face < 4 * gridNCell(grid) ; face++ )
    c2f[face] = EMPTY;

  nface = 0;
  for ( cell = 0; cell < gridNCell(grid) ; cell++ )
    for (cell_face = 0; cell_face < 4; cell_face++)
      if ( EMPTY == c2f[cell_face+4*cell] )
      {
	c2f[cell_face+4*cell] = nface; 
	if (grid != gridCell(grid,cell,cell_nodes))
	  { 
	    printf("ERROR: gridExportNGP: %s: %d: gridCell %d failed\n",
		   __FILE__, __LINE__, cell );
	    free(c2f); return NULL;
	  }
	face_nodes[0] = cell_nodes[cellfacenode[cell_face][0]];
	face_nodes[1] = cell_nodes[cellfacenode[cell_face][1]];
	face_nodes[2] = cell_nodes[cellfacenode[cell_face][2]];
	other_cell = gridFindOtherCellWith3Nodes(grid,
						 face_nodes[1],
						 face_nodes[0],
						 face_nodes[2],
						 cell );
	if ( EMPTY != other_cell )
	  {
	    gridCell(grid,other_cell,other_nodes);
	    for (other_face = 0; other_face < 4; other_face++)
	      {
		other_face_nodes[0] = other_nodes[cellfacenode[other_face][0]];
		other_face_nodes[1] = other_nodes[cellfacenode[other_face][1]];
		other_face_nodes[2] = other_nodes[cellfacenode[other_face][2]];
		if ( (face_nodes[0] == other_face_nodes[1]) &&
		     (face_nodes[1] == other_face_nodes[0]) &&
		     (face_nodes[2] == other_face_nodes[2]) ) 
		  c2f[other_face+4*other_cell] = nface;
		if ( (face_nodes[0] == other_face_nodes[2]) &&
		     (face_nodes[1] == other_face_nodes[1]) &&
		     (face_nodes[2] == other_face_nodes[0]) ) 
		  c2f[other_face+4*other_cell] = nface;
		if ( (face_nodes[0] == other_face_nodes[0]) &&
		     (face_nodes[1] == other_face_nodes[2]) &&
		     (face_nodes[2] == other_face_nodes[1]) ) 
		  c2f[other_face+4*other_cell] = nface;
	      }
	  }
	nface++;
      }

  if (NULL != filename) {
    file = fopen(filename,"w");
  }else{
    file = fopen("grid.ngp","w");
  }
  if (NULL == file) {
    printf("gridExportNGP: file open failed.\n");
    return NULL;
  }

  fprintf(file,"%d %d %d\n",gridNNode(grid),gridNCell(grid), nface);

  for( node=0; node<gridNNode(grid) ; node++ ) 
    fprintf(file,"%d %25.15e %25.15e %25.15e\n", (node+1),
	    grid->xyz[0+3*node], grid->xyz[1+3*node], grid->xyz[2+3*node]);

  fclose(file);

  free(c2f);

  return grid;
}

Grid *gridImportFAST( char *filename )
{
  FILE *file;
  int i, nnode, nface, maxcell, ncell;
  double *xyz;
  int *f2n, *faceId;
  int *c2n;

  file = fopen(filename,"r");
  if ( 3 != fscanf(file,"%d %d %d",&nnode,&nface,&ncell) ) return NULL;

  xyz = (double *)malloc(3*nnode*sizeof(double));

  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[0+3*i]) ) return NULL;
  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[1+3*i]) ) return NULL;
  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[2+3*i]) ) return NULL;

  f2n = (int *)malloc(3*nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    if ( 1 != fscanf(file,"%d",&f2n[0+3*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&f2n[1+3*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&f2n[2+3*i]) ) return NULL;
    f2n[0+3*i]--;
    f2n[1+3*i]--;
    f2n[2+3*i]--;
  }

  faceId = (int *)malloc(nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    if ( 1 != fscanf(file,"%d",&faceId[i]) ) return NULL;
  }

  maxcell = ncell*2;

  c2n = (int *)malloc(4*maxcell*sizeof(int));

  for( i=0; i<ncell ; i++ ) {
    if ( 1 != fscanf(file,"%d",&c2n[0+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[1+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[2+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[3+4*i]) ) return NULL;
    c2n[0+4*i]--;
    c2n[1+4*i]--;
    c2n[2+4*i]--;
    c2n[3+4*i]--;
  }

  fclose(file);

  return gridImport( nnode, nnode, nface, nface, maxcell, ncell, 0,
		     xyz, f2n, faceId, c2n );
}

Grid *gridExportFAST( Grid *grid, char *filename )
{
  FILE *file;
  int i;

  if (NULL == gridPack(grid)) {
    printf("gridExportFAST: gridPack failed.\n");
    return NULL;
  }

  if (NULL != filename) {
    file = fopen(filename,"w");
  }else{
    file = fopen("grid.fgrid","w");
  }

  fprintf(file,"%10d %10d %10d\n",grid->nnode,grid->nface,grid->ncell);

  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[0+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[1+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[2+3*i]);

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%10d %10d %10d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%4d\n",grid->faceId[i]);
  }

  for( i=0; i<grid->ncell ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);

  }

  fclose(file);

  return grid;
}

Grid *gridImportRef( char *filename )
{
  FILE *file;
  int i, nnode, nface, ncell;
  double *xyz;
  int *f2n, *faceId;
  int *c2n;
  Grid *grid;
  
  int nGeomNode, nGeomEdge, nGeomFace;
  int edge_id, edge_nnodes, *edge_nodes;
  int edge_id_verification;
  double *t;
  int node;
  int lines, uv_nnodes;
  double u, v;
  int face_id;

  file = fopen(filename,"r");
  if ( 3 != fscanf(file,"%d %d %d",&nnode,&nface,&ncell) ) return NULL;

  xyz = (double *)malloc(3*nnode*sizeof(double));

  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[0+3*i]) ) return NULL;
  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[1+3*i]) ) return NULL;
  for( i=0; i<nnode ; i++ ) 
    if ( 1 != fscanf(file,"%lf",&xyz[2+3*i]) ) return NULL;

  f2n = (int *)malloc(3*nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    if ( 1 != fscanf(file,"%d",&f2n[0+3*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&f2n[1+3*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&f2n[2+3*i]) ) return NULL;
    (f2n[0+3*i])--;
    (f2n[1+3*i])--;
    (f2n[2+3*i])--;
  }

  faceId = (int *)malloc(nface*sizeof(int));

  for( i=0; i<nface ; i++ ) {
    if ( 1 != fscanf(file,"%d",&faceId[i]) ) return NULL;
  }

  c2n = (int *)malloc(4*ncell*sizeof(int));

  for( i=0; i<ncell ; i++ ) {
    if ( 1 != fscanf(file,"%d",&c2n[0+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[1+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[2+4*i]) ) return NULL;
    if ( 1 != fscanf(file,"%d",&c2n[3+4*i]) ) return NULL;
    (c2n[0+4*i])--;
    (c2n[1+4*i])--;
    (c2n[2+4*i])--;
    (c2n[3+4*i])--;
  }

  grid = gridImport( nnode, nnode, nface, nface, ncell, ncell, 0,
		     xyz, f2n, faceId, c2n );

  if ( 3 != fscanf(file,"%d %d %d",&nGeomNode,&nGeomEdge,&nGeomFace) ) 
    return NULL;

  gridSetNGeomNode( grid, nGeomNode );
  gridSetNGeomEdge( grid, nGeomEdge );
  gridSetNGeomFace( grid, nGeomFace );

  for( edge_id=1; edge_id<=gridNGeomEdge(grid) ; edge_id++ ) {
    if ( 2 != fscanf( file,"%d %d", &edge_nnodes, &edge_id_verification) ) 
      return NULL;
    edge_nodes = (int *)malloc( edge_nnodes * sizeof(int) );
    t = (double *)malloc( edge_nnodes * sizeof(double) );
    for( i=0; i<edge_nnodes ; i++ ) {
      if ( 2 != fscanf( file,"%d %lf", &(edge_nodes[i]), &(t[i])) ) 
      return NULL;
    }
    gridAddGeomEdge( grid, edge_id, edge_nodes[0], edge_nodes[edge_nnodes-1]);
    for( i=0; i<(edge_nnodes-1) ; i++ ) {
      gridAddEdge(grid, edge_nodes[i], edge_nodes[i+1], edge_id,
		  t[i], t[i+1]);
    }

    free(edge_nodes);
    free(t);
  }
  
  if ( 1 != fscanf( file,"%d\n",&uv_nnodes) ) return NULL;
  for( lines=0; lines<uv_nnodes ; lines++ ) {
    if ( 4 != fscanf(file,"%d %d %lf %lf\n", 
	   &face_id, &node, &u, &v) ) return NULL;
    gridSetNodeUV(grid, node, face_id, u, v);
  }

  fclose(file);

  return grid;
}

Grid *gridExportRef( Grid *grid, char *filename )
{
  FILE *file;
  int i;
  int edge_id, edge_nnodes, *edge_nodes;
  double *t;
  int face, node;
  double uv[2];

  if (NULL == gridSortNodeGridEx( grid ) ) {
    printf("gridExportRef: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL != filename) {
    file = fopen(filename,"w");
  }else{
    file = fopen("grid.ref","w");
  }

  fprintf(file,"%10d %10d %10d\n",grid->nnode,grid->nface,grid->ncell);

  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[0+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[1+3*i]);
  for( i=0; i<grid->nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[2+3*i]);

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%10d %10d %10d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%4d\n",grid->faceId[i]);
  }

  for( i=0; i<grid->ncell ; i++ ) {
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);

  }

  fprintf( file,"%10d %10d %10d\n",
	   gridNGeomNode(grid),gridNGeomEdge(grid),gridNGeomFace(grid));

  for( edge_id=1; edge_id<=gridNGeomEdge(grid) ; edge_id++ ) {
    edge_nnodes = gridGeomEdgeSize( grid, edge_id );
    edge_nodes = (int *)malloc( edge_nnodes * sizeof(int) );
    t = (double *)malloc( edge_nnodes * sizeof(double) );
    gridGeomEdge( grid, edge_id, edge_nodes );
    gridGeomCurveT( grid, edge_id, edge_nodes[0], t );
    fprintf( file,"%10d %10d\n", edge_nnodes, edge_id);
    for( i=0; i<edge_nnodes ; i++ ) {
      fprintf( file,"%10d %24.16e\n", edge_nodes[i], t[i]);
    }
    free(edge_nodes);
    free(t);
  }

  fprintf( file,"%10d\n",3*grid->nface);
  for( face=0; face<grid->nface ; face++ ) {
    for( node=0; node<3 ; node++ ) {
      gridNodeUV(grid,grid->f2n[node+3*face], grid->faceId[face], uv);
      fprintf(file,"%10d %10d %24.16e %24.16e\n", 
	      grid->faceId[face], grid->f2n[node+3*face], uv[0], uv[1]);
    }
  }

  fclose(file);

  return grid;
}

Grid *gridExportFASTSurface( Grid *grid, char *filename )
{
  FILE *file;
  int i;
  int nnode;

  if (NULL == gridSortNodeGridEx(grid)) {
    printf("gridExportFAST: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL != filename) {
    file = fopen(filename,"w");
  }else{
    file = fopen("grid_surface.fgrid","w");
  }

  nnode = 0;
  for( i=0; i<grid->nface ; i++ ) {
    nnode = MAX(nnode,grid->f2n[0+3*i]);
    nnode = MAX(nnode,grid->f2n[1+3*i]);
    nnode = MAX(nnode,grid->f2n[2+3*i]);
  }
  nnode++;

  fprintf(file,"%10d %10d %10d\n",nnode,grid->nface,0);

  for( i=0; i<nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[0+3*i]);
  for( i=0; i<nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[1+3*i]);
  for( i=0; i<nnode ; i++ ) fprintf(file,"%25.15e\n",grid->xyz[2+3*i]);

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%10d %10d %10d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  for( i=0; i<grid->nface ; i++ ) {
    fprintf(file,"%4d\n",grid->faceId[i]);
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

Grid *gridImportGRI( char *filename )
{
  FILE *file;
  int i, nnode, nface, ncell;
  int cell_group_size;
  double *xyz;
  int *f2n, *faceId;
  int nface_groups;
  int group;
  int **face_group;
  int *face_group_size;
  int *c2n;
  GridBool verbose;

  verbose = TRUE;

  file = fopen(filename,"r");
  if ( 2 != fscanf(file,"%d %d",&nnode,&ncell) ) return NULL;

  if (verbose) printf("gri size: %d nodes %d cells.\n",nnode,ncell);

  if (verbose) printf("reading xyz...\n");
  
  xyz = (double *)malloc(3*nnode*sizeof(double));

  for( i=0; i<nnode ; i++ ) 
    if ( 3 != fscanf(file,"%lf %lf %lf",&xyz[0+3*i],&xyz[1+3*i],&xyz[2+3*i]) ) 
      return NULL;

  if ( 1 != fscanf(file,"%d",&nface_groups) ) return NULL;

  if (verbose) printf("gri size: %d face IDs.\n",nface_groups);

  face_group = (int **)malloc(nface_groups*sizeof(int *));
  face_group_size = (int *)malloc(nface_groups*sizeof(int));

  nface = 0;
  for (group=0;group<nface_groups;group++) {

    if (verbose) printf("reading faces in group %d...\n",group+1);
  
    if ( 1 != fscanf(file,"%d",&(face_group_size[group])) ) return NULL;
    nface += face_group_size[group];

    if (verbose) printf("gri size: %d group %d faces.\n",
			group+1,face_group_size[group]);

    f2n = (int *)malloc(3*face_group_size[group]*sizeof(int));
    face_group[group] = f2n;

    for( i=0; i<face_group_size[group] ; i++ ) {
      if ( 3 != fscanf(file,"%d %d %d",
		       &(f2n[0+3*i]),&(f2n[1+3*i]),&(f2n[2+3*i])) ) return NULL;
      f2n[0+3*i]--; f2n[1+3*i]--; f2n[2+3*i]--;
    }

  }

  f2n = (int *)malloc(3*nface*sizeof(int));
  faceId = (int *)malloc(nface*sizeof(int));

  nface = 0;
  for (group=0;group<nface_groups;group++) {
    for( i=0; i<face_group_size[group] ; i++ ) {
      f2n[0+3*nface] = face_group[group][0+3*i];
      f2n[1+3*nface] = face_group[group][1+3*i];
      f2n[2+3*nface] = face_group[group][2+3*i];
      faceId[nface] = group+1;
      nface++;
    }
  }

  /* PXGrid.c logic
  fgets(line, PX_MAXSTRLEN, fgri);
  if (fgets(line, PX_MAXSTRLEN, fgri) == NULL)
    return PXError(PX_FILE_READ_ERROR);

  ierr = sscanf(line, "%d %d", &nelem, &ntype);
  if (ierr != 2){
    ierr = sscanf(line, "%d", &nelem);
    if (ierr != 1) return PXError(PX_FILE_READ_ERROR);
    ntype = 1; // no specified type implies Q1
  }
   */

  if ( 1 != fscanf(file,"%d",&cell_group_size) ) return NULL;

  if (verbose) printf("gri size: %d cells (assume Q1 type).\n",
		      cell_group_size);

  if (verbose) printf("reading cells...\n");
  
  c2n = (int *)malloc(4*ncell*sizeof(int));

  for( i=0; i<ncell ; i++ ) {
    if ( 4 != fscanf(file,"%d %d %d %d",
		     &c2n[0+4*i],&c2n[1+4*i],&c2n[2+4*i],&c2n[3+4*i]) ) 
      return NULL;
    c2n[0+4*i]--; c2n[1+4*i]--; c2n[2+4*i]--; c2n[3+4*i]--;
  }

  fclose(file);

  return gridImport( nnode, nnode, nface, nface, ncell, ncell, 0,
		     xyz, f2n, faceId, c2n );
}

Grid *gridExportGRI( Grid *grid, char *filename )
{
  FILE *file;
  int i;
  int face, nFaceId, faceId, nface;

  if (NULL == gridPack(grid)) {
    printf("gridExportFAST: gridPack failed.\n");
    return NULL;
  }

  if (NULL != filename) {
    file = fopen(filename,"w");
  }else{
    file = fopen("grid.gri","w");
  }

  fprintf(file,"%10d %10d\n",grid->nnode,grid->ncell);

  for( i=0; i<grid->nnode ; i++ ) 
    fprintf(file,"%25.15e %25.15e %25.15e\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i]);

  nFaceId = 0;
  for ( face=0 ; face < grid->nface ; face++ ) 
    nFaceId = MAX(nFaceId, grid->faceId[face]);

  fprintf(file,"%10d\n",nFaceId);

  for (faceId=1;faceId<=nFaceId;faceId++) {
    nface = 0;
    for ( face=0 ; face < grid->nface ; face++ ) 
      if (faceId == grid->faceId[face]) nface++;
    fprintf(file,"%10d\n",nface);
    for ( face=0 ; face < grid->nface ; face++ ) 
      if (faceId == grid->faceId[face]) 
	fprintf(file,"%10d %10d %10d\n",
		grid->f2n[0+3*face]+1,
		grid->f2n[1+3*face]+1,
		grid->f2n[2+3*face]+1);
  }

  fprintf(file,"%10d %10d\n",grid->ncell,1);

  for( i=0; i<grid->ncell ; i++ )
    fprintf(file,"%10d %10d %10d %10d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);

  return grid;
}

Grid *gridImportAdapt( Grid *grid, char *filename )
{
  int i, nnode;
  FILE *file;

  file = fopen(filename,"r");

  if ( 1 != fscanf(file,"%d",&nnode) ) return NULL;

  printf("metric nnodes %10d %10d\n",nnode,gridNNode(grid));

  for( i=0; i<grid->nnode ; i++ ) 
    if ( 6 != fscanf(file,"%lf %lf %lf %lf %lf %lf\n",
       &grid->map[0+6*i], &grid->map[1+6*i], &grid->map[2+6*i],
                          &grid->map[3+6*i], &grid->map[4+6*i],
                                             &grid->map[5+6*i]) ) return NULL;
  fclose(file);
  return grid;
}

Grid *gridAttachPacker(Grid *grid, 
		       void (*packFunc)(void *packData, 
					int nnode, int maxnode, int *nodeo2n,
					int ncell, int maxcell, int *cello2n,
					int nface, int maxface, int *faceo2n,
					int nedge, int maxedge, int *edgeo2n),
		       void *packData )
{
  grid->packFunc = packFunc;
  grid->packData = packData;
  return grid;
}

Grid *gridDetachPacker(Grid *grid )
{
  grid->packFunc = NULL;
  grid->packData = NULL;
  return grid;
}

Grid *gridAttachNodeSorter(Grid *grid, 
			   void (*renumberFunc)(void *renumberData, 
						int maxnode, int *o2n), 
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

Grid *gridAttachReallocator(Grid *grid, 
			    void (*reallocFunc)(void *reallocData, 
						int reallocType, 
						int lastSize, int newSize), 
			    void *reallocData )
{
  grid->reallocFunc = reallocFunc;
  grid->reallocData = reallocData;
  return grid;
}

Grid *gridDetachReallocator(Grid *grid )
{
  grid->reallocFunc = NULL;
  grid->reallocData = NULL;
  return grid;
}

Grid *gridAttachFreeNotifier(Grid *grid, void (*freeNotificationFunc)
			     (void *freeNotificationData),
			     void *freeNotificationData)
{
  grid->freeNotificationFunc = freeNotificationFunc;
  grid->freeNotificationData = freeNotificationData;
  return grid;
}

Grid *gridDetachFreeNotifier(Grid *grid)
{
  grid->freeNotificationFunc = NULL;
  grid->freeNotificationData = NULL;
  return grid;
}


void gridFree(Grid *grid)
{
  if (NULL != grid->freeNotificationFunc) 
    (*grid->freeNotificationFunc)( grid->freeNotificationData );
  if (NULL != grid->lines) linesFree(grid->lines);

  gridCloseTecplotGeomFile(grid);
  gridCloseTecplotScalarFile(grid);

  if (NULL != grid->conn2node) free(grid->conn2node);
  if (NULL != grid->cell2conn) free(grid->cell2conn);

  if (NULL != grid->geomEdge) free(grid->geomEdge);
  if (NULL != grid->geomNode) free(grid->geomNode);

  if (NULL != grid->prismDeg) free(grid->prismDeg);

  if (NULL != grid->quad) free(grid->quad);
  if (NULL != grid->pyramid) free(grid->pyramid);
  if (NULL != grid->prism) free(grid->prism);
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
  if (NULL != grid->unusedCellGlobal) free(grid->unusedCellGlobal);
  if (NULL != grid->cellGlobal) free(grid->cellGlobal);
  free(grid->c2n);
  if (NULL != grid->unusedNodeGlobal) free(grid->unusedNodeGlobal);
  if (NULL != grid->sortedLocal) free(grid->sortedLocal);
  if (NULL != grid->sortedGlobal) free(grid->sortedGlobal);
  if (NULL != grid->part) free(grid->part);
  if (NULL != grid->nodeGlobal) free(grid->nodeGlobal);
  if (NULL != grid->child_reference) free(grid->child_reference);
  if (NULL != grid->aux) free(grid->aux);
  free(grid->frozen);
  if (NULL != grid->map) free(grid->map);
  free(grid->xyz);
  if (NULL != grid->child) gridFree(grid->child);
  free(grid);
}

Grid *gridPack(Grid *grid)
{
  int i;
  int orignode, packnode, origcell, packcell;
  int origface, packface, origedge, packedge;
  int conn;
  int nFaceId;
  int iface, n0, n1, id, n[3];
  double t0, t1, u[3], v[3];
  int *nodeo2n, *cello2n, *faceo2n, *edgeo2n;
  int prismIndex, pyramidIndex, quadIndex;

  nodeo2n = (int *)malloc(grid->maxnode*sizeof(int));
  for (i=0;i<grid->maxnode;i++) nodeo2n[i] = EMPTY;
  cello2n = (int *)malloc(grid->maxcell*sizeof(int));
  for (i=0;i<grid->maxcell;i++) cello2n[i] = EMPTY;
  faceo2n = (int *)malloc(grid->maxface*sizeof(int));
  for (i=0;i<grid->maxface;i++) faceo2n[i] = EMPTY;
  edgeo2n = (int *)malloc(grid->maxedge*sizeof(int));
  for (i=0;i<grid->maxedge;i++) edgeo2n[i] = EMPTY;

  packnode = 0;
  for ( orignode=0 ; orignode < grid->maxnode ; orignode++ )
    if ( grid->xyz[0+3*orignode] != DBL_MAX) {
      nodeo2n[orignode] = packnode;
      grid->xyz[0+3*packnode] = grid->xyz[0+3*orignode];
      grid->xyz[1+3*packnode] = grid->xyz[1+3*orignode];
      grid->xyz[2+3*packnode] = grid->xyz[2+3*orignode];
      if ( NULL != grid->map ) { 
	grid->map[0+6*packnode] = grid->map[0+6*orignode];
	grid->map[1+6*packnode] = grid->map[1+6*orignode];
	grid->map[2+6*packnode] = grid->map[2+6*orignode];
	grid->map[3+6*packnode] = grid->map[3+6*orignode];
	grid->map[4+6*packnode] = grid->map[4+6*orignode];
	grid->map[5+6*packnode] = grid->map[5+6*orignode];
      }
      grid->frozen[packnode]  = grid->frozen[orignode];
      if (NULL != grid->child_reference) 
	grid->child_reference[packnode] = grid->child_reference[orignode];
      if (NULL != grid->nodeGlobal) 
	grid->nodeGlobal[packnode] = grid->nodeGlobal[orignode];
      if (NULL != grid->part) grid->part[packnode] = grid->part[orignode];
      for ( i=0;i<grid->naux;i++) 
	grid->aux[i+grid->naux*packnode] = grid->aux[i+grid->naux*orignode];
      grid->geomNode[packnode] = grid->geomNode[orignode];
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

  if (NULL != grid->sortedLocal)
    for (i=0 ; i < grid->nsorted; i++ )
      grid->sortedLocal[i] = nodeo2n[grid->sortedLocal[i]];

  packcell = 0;
  for ( origcell=0 ; origcell < grid->maxcell ; origcell++ )
    if ( grid->c2n[0+4*origcell] != EMPTY) {
      cello2n[origcell] = packcell;
      adjRemove( grid->cellAdj, grid->c2n[0+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[1+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[2+4*origcell], origcell );
      adjRemove( grid->cellAdj, grid->c2n[3+4*origcell], origcell );
      grid->c2n[0+4*packcell] = nodeo2n[grid->c2n[0+4*origcell]];
      grid->c2n[1+4*packcell] = nodeo2n[grid->c2n[1+4*origcell]];
      grid->c2n[2+4*packcell] = nodeo2n[grid->c2n[2+4*origcell]];
      grid->c2n[3+4*packcell] = nodeo2n[grid->c2n[3+4*origcell]];
      if (NULL != grid->cellGlobal) 
	grid->cellGlobal[packcell] = grid->cellGlobal[origcell];
      adjRegister( grid->cellAdj, grid->c2n[0+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[1+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[2+4*packcell], packcell );
      adjRegister( grid->cellAdj, grid->c2n[3+4*packcell], packcell );
      if (NULL != grid->cell2conn)
	for(conn=0;conn<6;conn++)
	  grid->cell2conn[conn+6*packcell] = grid->cell2conn[conn+6*origcell];
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

  for(conn=0;conn<gridNConn(grid);conn++) {
    if (EMPTY!=grid->conn2node[0+2*conn]) 
      grid->conn2node[0+2*conn] = nodeo2n[grid->conn2node[0+2*conn]];
    if (EMPTY!=grid->conn2node[1+2*conn]) 
      grid->conn2node[1+2*conn] = nodeo2n[grid->conn2node[1+2*conn]];
  }

  packface=0;

  /* gridSortNodeGridEx needs the faceIds in ascending order */

  nFaceId = 0;
  for ( origface=0 ; origface < grid->maxface ; origface++ ) 
    if (grid->f2n[0+3*origface] != EMPTY)
      nFaceId = MAX(nFaceId, grid->faceId[origface]);

  iface=0;
  for (iface=1;iface <= nFaceId;iface++){
    for ( origface=0 ; origface < grid->maxface ; origface++ ){ 
      if ( grid->faceId[origface]==iface && grid->f2n[0+3*origface] != EMPTY) {
	faceo2n[origface] = packface;
	if (origface == packface) {
	  for (i=0;i<3;i++){
	    adjRemove( grid->faceAdj, grid->f2n[i+3*packface], packface );
	    grid->f2n[i+3*packface] = nodeo2n[grid->f2n[i+3*packface]];
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
	    grid->f2n[i+3*packface] = nodeo2n[n[i]];
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
      edgeo2n[origedge] = packedge;
      n0 = nodeo2n[grid->e2n[0+2*origedge]];
      n1 = nodeo2n[grid->e2n[1+2*origedge]];
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
			nodeo2n[grid->prism[prismIndex].nodes[i]];
  }

  for (pyramidIndex=0;pyramidIndex<gridNPyramid(grid);pyramidIndex++){
    for (i=0;i<5;i++) grid->pyramid[pyramidIndex].nodes[i] =
			nodeo2n[grid->pyramid[pyramidIndex].nodes[i]];
  }

  /* note, these should be counted as boundaries */
  for (quadIndex=0;quadIndex<gridNQuad(grid);quadIndex++){
    for (i=0;i<4;i++) grid->quad[quadIndex].nodes[i] =
			nodeo2n[grid->quad[quadIndex].nodes[i]];
  }

  if ( NULL != grid->packFunc ) 
    (*grid->packFunc)( grid->packData, 
		       grid->nnode, grid->maxnode, nodeo2n,
		       grid->ncell, grid->maxcell, cello2n,
		       grid->nface, grid->maxface, faceo2n,
		       grid->nedge, grid->maxedge, edgeo2n);


  if ( NULL != gridLines(grid) ) linesRenumber(gridLines(grid),nodeo2n);

  free(nodeo2n);
  free(cello2n);
  free(faceo2n);
  free(edgeo2n);

  return  grid;
}

Grid *gridSortNodeGridEx(Grid *grid)
{
  int i, newnode, edge, nCurveNode;
  int node, face;
  int *o2n, *curve;

  if (NULL == gridPack(grid)) {
    printf("gridSortNodeGridEx: gridPack failed.\n");
    return NULL;
  }

  o2n = (int *)malloc( grid->maxnode * sizeof(int) );
  for (i=0;i<grid->maxnode;i++) o2n[i] = EMPTY;

  /* geom nodes */
  newnode = 0;
  for (node=0;node<grid->nnode;node++) {
    if ( EMPTY != grid->geomNode[node] ) {
      o2n[node] = newnode;
      newnode++;
    }
  }

  /* edge stuff */
  for (edge=1; edge<=grid->nGeomEdge; edge++){

    nCurveNode = gridGeomEdgeSize( grid, edge );
    if (nCurveNode > 0) { /* do not malloc size 0 for off proc edges */
      curve = (int *)malloc( nCurveNode * sizeof(int) );
      gridGeomEdge( grid, edge, curve );

      for ( i=1; i<(nCurveNode-1); i++){ /* skip end points */
	if (o2n[curve[i]] != EMPTY) 
	  printf("gridSortNodeGridEx: %s: %d: newnode error %d\n",
		 __FILE__, __LINE__, o2n[curve[i]] );
	o2n[curve[i]] = newnode;
	newnode++;
      }

      free(curve);
    }
  }

  /* face stuff
   *   assuming that the bc faces are sorted for contiguous
   *   face node numbering by gridPack. */

  for ( face=0; face<grid->nface; face++ ){
    for ( i=0; i<3; i++ ){
      node = grid->f2n[i+3*face];
      if ( o2n[node] == EMPTY ) {
	o2n[node] = newnode;
	newnode++;
      }
    }
  }

  /* interior nodes */
  for ( node=0; node<grid->nnode; node++ ){
    if ( o2n[node] == EMPTY ) {
      o2n[node] = newnode;
      newnode++;
    }
  }

  if (newnode != grid->nnode) 
    printf("ERROR: gridSortNodeGridEx, newnode %d nnode %d, line %d of %s\n.",
	   newnode,grid->nnode,__LINE__, __FILE__);

  gridRenumber(grid, o2n);

  free(o2n);

  return grid;
}

Grid *gridSortNodeFUN3D(Grid *grid, int *nnodes0)
{
  int i, newnode;
  int node;
  int *o2n;

  if (NULL == gridPack(grid)) {
    printf("gridSortNodeFUN3D: gridPack failed.\n");
    return NULL;
  }

  o2n = (int *)malloc( grid->maxnode * sizeof(int) );
  for (i=0;i<grid->maxnode;i++) o2n[i] = EMPTY;

  newnode = 0;

  /* local nodes */
  for ( node=0; node<grid->nnode; node++ ){
    if ( gridNodeLocal(grid,node) ) {
      o2n[node] = newnode;
      newnode++;
    }
  }

  *nnodes0 = newnode;

  /* interior nodes */
  for ( node=0; node<grid->nnode; node++ ){
    if ( o2n[node] == EMPTY ) {
      o2n[node] = newnode;
      newnode++;
    }
  }
  if (newnode != grid->nnode) 
    printf("ERROR: gridSortNodeFUN3D, newnode %d nnode %d, line %d of %s\n.",
	   newnode,grid->nnode,__LINE__, __FILE__);

  gridRenumber(grid, o2n);

  free(o2n);

  return grid;
}

Grid *gridRenumber(Grid *grid, int *o2n)
{
  int i, edge;
  int ixyz, node, face, cell, inode;
  int conn;
  double *temp_xyz;
  GridBool *temp_frozen;
  int *temp_int;
  int prismIndex, pyramidIndex, quadIndex;

  temp_xyz = (double *)malloc( grid->nnode * sizeof(double) );
  for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
    for ( node = 0 ; node < grid->nnode ; node++ ){
      temp_xyz[o2n[node]] = grid->xyz[ixyz+3*node];
    }
    for ( node = 0 ; node < grid->nnode ; node++ ){
      grid->xyz[ixyz+3*node] = temp_xyz[node];
    }
  }
  if ( NULL != grid->map ) {
    for ( ixyz = 0; ixyz < 6 ; ixyz++ ){
      for ( node = 0 ; node < grid->nnode ; node++ ){
	temp_xyz[o2n[node]] = grid->map[ixyz+6*node];
      }
      for ( node = 0 ; node < grid->nnode ; node++ ){
	grid->map[ixyz+6*node] = temp_xyz[node];
      }
    }
  }
  for ( ixyz = 0; ixyz < grid->naux ; ixyz++ ){
    for ( node = 0 ; node < grid->nnode ; node++ ){
      temp_xyz[o2n[node]] = grid->aux[ixyz+grid->naux*node];
    }
    for ( node = 0 ; node < grid->nnode ; node++ ){
      grid->aux[ixyz+grid->naux*node] = temp_xyz[node];
    }
  }
  free(temp_xyz);

  temp_frozen = (GridBool *)malloc( grid->nnode * sizeof(GridBool) );
  for ( node = 0 ; node < grid->nnode ; node++ )
    temp_frozen[o2n[node]] = grid->frozen[node];
  for ( node = 0 ; node < grid->nnode ; node++ )
    grid->frozen[node] = temp_frozen[node];
  free(temp_frozen);

  temp_int = (int *)malloc( grid->nnode * sizeof(int) );
  for ( node = 0 ; node < grid->nnode ; node++ )
    temp_int[o2n[node]] = grid->geomNode[node];
  for ( node = 0 ; node < grid->nnode ; node++ )
    grid->geomNode[node] = temp_int[node];
  if (NULL != grid->child_reference) {
    for ( node = 0 ; node < grid->nnode ; node++ )
      temp_int[o2n[node]] = grid->child_reference[node];
    for ( node = 0 ; node < grid->nnode ; node++ )
      grid->child_reference[node] = temp_int[node];
  }
  if (NULL != grid->nodeGlobal) {
    for ( node = 0 ; node < grid->nnode ; node++ )
      temp_int[o2n[node]] = grid->nodeGlobal[node];
    for ( node = 0 ; node < grid->nnode ; node++ )
      grid->nodeGlobal[node] = temp_int[node];
  }
  if (NULL != grid->part) {
    for ( node = 0 ; node < grid->nnode ; node++ )
      temp_int[o2n[node]] = grid->part[node];
    for ( node = 0 ; node < grid->nnode ; node++ )
      grid->part[node] = temp_int[node];
  }
  free(temp_int);

  if (NULL != grid->sortedLocal)
    for (i=0 ; i < grid->nsorted; i++ )
      grid->sortedLocal[i] = o2n[grid->sortedLocal[i]];


  for ( cell = 0; cell < grid->ncell ; cell++ ){
    for ( inode = 0 ; inode < 4 ; inode++ ){
      adjRemove( grid->cellAdj, grid->c2n[inode+4*cell], cell );
      grid->c2n[inode+4*cell] = o2n[grid->c2n[inode+4*cell]];
      adjRegister( grid->cellAdj, grid->c2n[inode+4*cell], cell );
    }
  }

  for(conn=0;conn<gridNConn(grid);conn++) {
    if (EMPTY!=grid->conn2node[0+2*conn]) 
      grid->conn2node[0+2*conn] = o2n[grid->conn2node[0+2*conn]];
    if (EMPTY!=grid->conn2node[1+2*conn]) 
      grid->conn2node[1+2*conn] = o2n[grid->conn2node[1+2*conn]];
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

  /* note, these should be counted as boundaries */
  for (quadIndex=0;quadIndex<gridNQuad(grid);quadIndex++){
    for (i=0;i<4;i++) grid->quad[quadIndex].nodes[i] =
			o2n[grid->quad[quadIndex].nodes[i]];
  }

  if ( NULL != grid->renumberFunc ) 
    (*grid->renumberFunc)( grid->renumberData, grid->maxnode, o2n );

  if ( NULL != gridLines(grid) ) linesRenumber(gridLines(grid),o2n);

  return grid;
}

Grid *gridWriteTecplotSurfaceGeom(Grid *grid, char *filename)
{
  int i, nfacenode;
  if ( grid !=  gridSortNodeGridEx(grid) ) {
    printf("gridWriteTecplotSurfaceGeom: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("gridGeom.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  nfacenode=0;
  for(i=0;i<3*grid->nface;i++){
    nfacenode = MAX(nfacenode, grid->f2n[i]);
  }
  nfacenode++;

  fprintf(grid->tecplotGeomFile,
	  "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  nfacenode, grid->nface);

  for ( i=0; i<nfacenode ; i++ ){
    AdjIterator it;
    int nodes[3], face, faceId;
    it = adjFirst(gridFaceAdj(grid),i);
    face = adjItem(it);
    gridFace(grid, face, nodes, &faceId);
    fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e %d\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i],faceId);
  }

  fprintf(grid->tecplotGeomFile, "\n");

  for ( i=0; i<grid->nface ; i++ ){
    fprintf(grid->tecplotGeomFile, " %9d %9d %9d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridWriteTecplotVolumeGeom(Grid *grid, char *filename)
{
  int i, nfacenode;
  if ( grid !=  gridSortNodeGridEx(grid) ) {
    printf("gridWriteTecplotVolumeGeom: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("gridGeom.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  nfacenode=0;
  for(i=0;i<3*grid->nface;i++){
    nfacenode = MAX(nfacenode, grid->f2n[i]);
  }
  nfacenode++;

  fprintf(grid->tecplotGeomFile,
	  "zone t=volume, i=%d, j=%d, f=fepoint, et=tetrahedron\n",
	  grid->nnode, grid->ncell);

  for ( i=0; i<grid->nnode ; i++ ){
    fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e %d\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i],0);
  }

  fprintf(grid->tecplotGeomFile, "\n");

  for ( i=0; i<grid->ncell ; i++ ){
    fprintf(grid->tecplotGeomFile, " %9d %9d %9d %9d\n",
	    grid->c2n[0+4*i]+1,grid->c2n[1+4*i]+1,
	    grid->c2n[2+4*i]+1,grid->c2n[3+4*i]+1);
  }

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridWriteTecplotGeomFaceUV(Grid *grid, char *filename, int id )
{
  int face, nodes[3], faceId;
  int nface, nnode;
  int i, node;
  double *xyz;
  int *f2n;
  int *g2l;
  Grid *status;
  char zone[256];

  nface = gridTrianglesOnFaceId(grid,id);
  if (0==nface) return NULL;

  f2n = (int *)malloc(sizeof(int)*3*nface);
  nface = 0;
  for (face=0; face<gridMaxFace(grid); face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      if (id == faceId) {
	f2n[0+3*nface] = nodes[0];
	f2n[1+3*nface] = nodes[1];
	f2n[2+3*nface] = nodes[2];
	nface++;
      }
    }
  }

  g2l = (int *)malloc(sizeof(int)*gridMaxNode(grid));
  for (node=0; node<gridMaxNode(grid); node++) g2l[node]=EMPTY;
  nnode = 0;
  for (face=0; face<nface; face++) {
    for (i=0;i<3;i++) {
      if (EMPTY == g2l[f2n[i+3*face]]) {
	g2l[f2n[i+3*face]] = nnode;
	nnode++;
      }
    }
  }

  xyz = (double *)malloc(sizeof(double)*3*nnode);
  for (node=0; node<3*nnode; node++) xyz[node]=0.0;
  for (face=0; face<nface; face++) {
    for (i=0;i<3;i++) {
      gridNodeUV(grid,f2n[i+3*face],id,&(xyz[3*g2l[f2n[i+3*face]]]));
      f2n[i+3*face] = g2l[f2n[i+3*face]];
    }
  }
  free(g2l);

  sprintf(zone, "face_%04d", id);

  status = gridWriteTecplotTriangleZone(grid, filename, zone,
					nnode, xyz,
					nface, f2n);
  free(f2n);
  free(xyz);

  return status;
}

Grid *gridWriteTecplotGeomFaceXYZ(Grid *grid, char *filename, int id )
{
  int face, nodes[3], faceId;
  int nface, nnode;
  int i, node;
  double *xyz;
  int *f2n;
  int *g2l;
  Grid *status;
  char zone[256];

  nface = gridTrianglesOnFaceId(grid,id);
  if (0==nface) return NULL;

  f2n = (int *)malloc(sizeof(int)*3*nface);
  nface = 0;
  for (face=0; face<gridMaxFace(grid); face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      if (id == faceId) {
	f2n[0+3*nface] = nodes[0];
	f2n[1+3*nface] = nodes[1];
	f2n[2+3*nface] = nodes[2];
	nface++;
      }
    }
  }

  g2l = (int *)malloc(sizeof(int)*gridMaxNode(grid));
  for (node=0; node<gridMaxNode(grid); node++) g2l[node]=EMPTY;
  nnode = 0;
  for (face=0; face<nface; face++) {
    for (i=0;i<3;i++) {
      if (EMPTY == g2l[f2n[i+3*face]]) {
	g2l[f2n[i+3*face]] = nnode;
	nnode++;
      }
    }
  }

  xyz = (double *)malloc(sizeof(double)*3*nnode);
  for (node=0; node<3*nnode; node++) xyz[node]=0.0;
  for (face=0; face<nface; face++) {
    for (i=0;i<3;i++) {
      gridNodeXYZ(grid,f2n[i+3*face],&(xyz[3*g2l[f2n[i+3*face]]]));
      f2n[i+3*face] = g2l[f2n[i+3*face]];
    }
  }
  free(g2l);

  sprintf(zone, "face_%04d", id);

  status = gridWriteTecplotTriangleZone(grid, filename, zone,
					nnode, xyz,
					nface, f2n);
  free(f2n);
  free(xyz);

  return status;
}

Grid *gridWriteTecplotTriangleZone(Grid *grid, char *filename, char *zone,
				   int nnode, double *xyz,
				   int nface, int *f2n)
{
  int i;

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("gridGeom.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  if ( NULL == zone ) {
    fprintf(grid->tecplotGeomFile,
	    "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	    nnode, nface);
  }else{
    fprintf(grid->tecplotGeomFile,
	    "zone t=%s, i=%d, j=%d, f=fepoint, et=triangle\n",
	    zone, nnode, nface);
  }

  for ( i=0; i<nnode ; i++ ){
    fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e %d\n",
	    xyz[0+3*i],xyz[1+3*i],xyz[2+3*i],0);
  }

  for ( i=0; i<nface ; i++ ){
    fprintf(grid->tecplotGeomFile, " %9d %9d %9d\n",
	    f2n[0+3*i]+1,f2n[1+3*i]+1,f2n[2+3*i]+1);
  }

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridWriteTecplotComment(Grid *grid, char *comment)
{
  if (NULL == grid->tecplotGeomFile) return NULL;
  fprintf(grid->tecplotGeomFile, "# %s\n", comment);
  fflush(grid->tecplotGeomFile);
  return grid;
}

Grid *gridWriteTecplotCellGeom(Grid *grid, int *nodes, double *scalar, 
			       char *filename)
{
  int i;
  double xyz[3];

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("grid.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  fprintf(grid->tecplotGeomFile, "zone t=cell, n=%d, e=%d, f=fepoint, et=%s\n",
	  4, 1, "tetrahedron");

  for ( i=0; i<4 ; i++ ){
    gridNodeXYZ(grid,nodes[i],xyz);
    if ( NULL == scalar ) {
      fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e %d\n",
	      xyz[0],xyz[1],xyz[2],0);
    }else{
      fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e%23.15e\n",
	      xyz[0],xyz[1],xyz[2],scalar[i]);
    }
  }

  fprintf(grid->tecplotGeomFile, "1 2 3 4\n");

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridWriteTecplotNodeOrbit(Grid *grid, int node, char *filename )
{
  AdjIterator it;
  int i, cell, nodes[4];
  double xyz[3];

  if ( !gridValidNode(grid, node) ) return NULL;

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("grid_orbit.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  fprintf(grid->tecplotGeomFile,
	  "zone t=orbit, n=%d, e=%d, f=fepoint, et=%s\n",
	  4*gridCellDegree(grid,node), gridCellDegree(grid,node), 
	  "tetrahedron");

  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    cell = adjItem(it);
    if ( grid != gridCell(grid, cell, nodes ) ) return NULL;
    for ( i = 0 ; i < 4 ; i++ ) {
      gridNodeXYZ(grid,nodes[i],xyz);
      fprintf(grid->tecplotGeomFile,
	      "%23.15e%23.15e%23.15e%5.1f\n", 
	      xyz[0],xyz[1],xyz[2],(double)gridGeometryFace(grid,nodes[i]));
    }
  } 

  i = 1; /* tecplot uses 1-base indexes */
  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    fprintf(grid->tecplotGeomFile, "%6d%6d%6d%6d\n", i+0, i+1, i+2, i+3 );
    i += 4;
  } 

  fflush(grid->tecplotGeomFile);

  return grid;
}
Grid *gridWriteTecplotEquator(Grid *grid, int n0, int n1, char *filename )
{
  int i, last;
  double xyz[3];

  if (grid != gridEquator(grid, n0, n1 )) return NULL;

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("grid_equator.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  fprintf(grid->tecplotGeomFile,
	  "zone t=equator, n=%d, e=%d, f=fepoint, et=%s\n",
	  2 + gridNEqu(grid), gridNGem(grid), "tetrahedron");

  gridNodeXYZ(grid,n0,xyz);
  fprintf(grid->tecplotGeomFile,
	  "%23.15e%23.15e%23.15e%5.1f\n", xyz[0],xyz[1],xyz[2],0.0);
  gridNodeXYZ(grid,n1,xyz);
  fprintf(grid->tecplotGeomFile,
	  "%23.15e%23.15e%23.15e%5.1f\n", xyz[0],xyz[1],xyz[2],1.0);

  for ( i=0; i<gridNEqu(grid) ; i++ ){
    gridNodeXYZ(grid,gridEqu(grid,i),xyz);
    fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e%5.1f\n",xyz[0],xyz[1],xyz[2],0.5);
  }
 

  for ( i=0; i<gridNGem(grid) ; i++ ){
    last = 4+i;
    if (last>(gridNEqu(grid)+2)) last = 3;
    fprintf(grid->tecplotGeomFile, "%6d%6d%6d%6d\n", 1, 2, 3+i, last );
  }

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridWriteTecplotEquatorFaces(Grid *grid, int n0, int n1, char *filename )
{
  int i;
  double xyz[3];

  if (grid != gridEquator(grid, n0, n1 )) return NULL;

  if (NULL == grid->tecplotGeomFile) {
    if (NULL == filename) {
      grid->tecplotGeomFile = fopen("grid_equator.t","w");
    }else{
      grid->tecplotGeomFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotGeomFile, "title=\"tecplot refine geometry file\"\n");
    fprintf(grid->tecplotGeomFile, "variables=\"X\",\"Y\",\"Z\",\"Face\"\n");
  }

  fprintf(grid->tecplotGeomFile,
	  "zone t=equator, n=%d, e=%d, f=fepoint, et=%s\n",
	  gridNEqu(grid)+2, gridNEqu(grid), "triangle");

  gridNodeXYZ(grid,n0,xyz);
  fprintf(grid->tecplotGeomFile,
	  "%23.15e%23.15e%23.15e%5.1f\n", xyz[0],xyz[1],xyz[2],0.0);
  gridNodeXYZ(grid,n1,xyz);
  fprintf(grid->tecplotGeomFile,
	  "%23.15e%23.15e%23.15e%5.1f\n", xyz[0],xyz[1],xyz[2],1.0);

  for ( i=0; i<gridNEqu(grid) ; i++ ){
    gridNodeXYZ(grid,gridEqu(grid,i),xyz);
    fprintf(grid->tecplotGeomFile, "%23.15e%23.15e%23.15e%5.1f\n",xyz[0],xyz[1],xyz[2],0.5);
  }
  

  for ( i=0; i<gridNEqu(grid) ; i++ ){
    fprintf(grid->tecplotGeomFile, "%6d%6d%6d\n", 1, 2, 3+i );
  }

  fflush(grid->tecplotGeomFile);

  return grid;
}

Grid *gridCloseTecplotGeomFile(Grid *grid)
{
  if (NULL == grid->tecplotGeomFile) return NULL;

  fclose(grid->tecplotGeomFile);
  grid->tecplotGeomFile = NULL;
  return grid;
}

Grid *gridWriteTecplotSurfaceScalar(Grid *grid, char *filename, double *scalar)
{
  int i, nfacenode;
  if ( grid !=  gridSortNodeGridEx(grid) ) {
    printf("gridWriteTecplotSurfaceScalar: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL == grid->tecplotScalarFile) {
    if (NULL == filename) {
      grid->tecplotScalarFile = fopen("gridScalar.t","w");
    }else{
      grid->tecplotScalarFile = fopen(filename,"w");
    } 
    fprintf(grid->tecplotScalarFile, "title=\"tecplot refine scalar file\"\n");
    fprintf(grid->tecplotScalarFile, "variables=\"X\",\"Y\",\"Z\",\"S\"\n");
  }

  nfacenode=0;
  for(i=0;i<3*grid->nface;i++){
    nfacenode = MAX(nfacenode, grid->f2n[i]);
  }
  nfacenode++;

  fprintf(grid->tecplotScalarFile, "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  nfacenode, grid->nface);

  for ( i=0; i<nfacenode ; i++ ){
    fprintf(grid->tecplotScalarFile, "%23.15e%23.15e%23.15e%23.15e\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i],scalar[i]);
  }

  fprintf(grid->tecplotScalarFile, "\n");

  for ( i=0; i<grid->nface ; i++ ){
    fprintf(grid->tecplotScalarFile, " %9d %9d %9d\n",
	    grid->f2n[0+3*i]+1,grid->f2n[1+3*i]+1,grid->f2n[2+3*i]+1);
  }

  fflush(grid->tecplotScalarFile);

  return grid;
}

Grid *gridCloseTecplotScalarFile(Grid *grid)
{
  if (NULL == grid->tecplotScalarFile) return NULL;

  fclose(grid->tecplotScalarFile);
  grid->tecplotScalarFile = NULL;
  return grid;
}

Grid *gridWriteVTK(Grid *grid, char *filename)
{
  FILE *vtk;
  int i;
  if ( grid !=  gridSortNodeGridEx(grid) ) {
    printf("gridWriteVTK: gridSortNodeGridEx failed.\n");
    return NULL;
  }

  if (NULL == filename) {
    vtk = fopen("grid.vtu","w");
  }else{
    vtk = fopen(filename,"w");
  } 
  fprintf(vtk, "<?xml version=\"1.0\"?>\n");
  fprintf(vtk, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(vtk, "  <UnstructuredGrid>\n");

  fprintf(vtk, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
	  gridNNode(grid), gridNCell(grid));

  fprintf(vtk, "      <Points Scalars=\"my_scalars\">\n");
  fprintf(vtk, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for ( i=0; i<gridNNode(grid) ; i++ ){
    fprintf(vtk, "           %23.15e%23.15e%23.15e\n",
	    grid->xyz[0+3*i],grid->xyz[1+3*i],grid->xyz[2+3*i]);
  }
  fprintf(vtk, "        </DataArray>\n");
  fprintf(vtk, "      </Points>\n");
  fprintf(vtk, "      <Cells>\n");

  fprintf(vtk, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n","connectivity");
  for ( i=0; i<gridNCell(grid) ; i++ ){
    fprintf(vtk, "           %10d %10d %10d %10d\n",
	    grid->c2n[0+4*i],grid->c2n[1+4*i],
	    grid->c2n[2+4*i],grid->c2n[3+4*i]);
  }
  fprintf(vtk, "        </DataArray>\n");

  fprintf(vtk, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n","offsets");
  fprintf(vtk, "          ");
  for ( i=0; i<gridNCell(grid) ; i++ ) fprintf(vtk, " %d",(i+1)*4);
  fprintf(vtk, "\n");
  fprintf(vtk, "        </DataArray>\n");


  fprintf(vtk, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n","types");
  fprintf(vtk, "          ");
  for ( i=0; i<gridNCell(grid) ; i++ ) fprintf(vtk, " %d",10);
  fprintf(vtk, "\n");
  fprintf(vtk, "        </DataArray>\n");

  fprintf(vtk, "      </Cells>\n");
  fprintf(vtk, "    </Piece>\n");
  fprintf(vtk, "  </UnstructuredGrid>\n");
  fprintf(vtk, "</VTKFile>\n");

  fclose(vtk);

  return grid;
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

int gridPartId(Grid *grid)
{
  return grid->partId;
}

int gridNAux(Grid *grid)
{
  return grid->naux;
}

Grid *gridSetNAux(Grid *grid, int naux )
{
  grid->naux = MAX(0,naux);
  if (NULL != grid->aux) free(grid->aux);
  if (grid->naux>0) 
    grid->aux = (double *)malloc( grid->maxnode * grid->naux * sizeof(double) );
  return grid;
}

double gridAux(Grid *grid, int node, int aux)
{
  if (aux<0  || aux >= grid->naux) return 0.0;
  if (node<0 || node >= grid->maxnode) return 0.0;
  
  return grid->aux[aux+grid->naux*node];
}

Grid *gridSetAux(Grid *grid, int node, int aux, double value )
{
  if (aux<0  || aux >= grid->naux) return NULL;
  if (node<0 || node >= grid->maxnode) return NULL;
  grid->aux[aux+grid->naux*node] = value;
  return grid;
}

Grid *gridInterpolateAux2(Grid *grid, int node0, int node1, double ratio, 
			  int target)
{
  int i;
  int n;

  if ( !gridValidNode(grid, node0) ) return NULL;
  if ( !gridValidNode(grid, node1) ) return NULL;
  if ( !gridValidNode(grid, target) ) return NULL;

  n = grid->naux;

  for (i=0;i<n;i++) {
    grid->aux[i+n*target] = 
      ratio*(grid->aux[i+n*node1]) + (1.0-ratio)*(grid->aux[i+n*node0]); 
  }
  return grid;
}



Grid *gridSetAuxToAverageOfNodes2(Grid *grid, int avgNode,
				  int n0, int n1 )
{
  int aux;
  if ( !gridValidNode(grid,n0) || 
       !gridValidNode(grid,n1) ) return NULL;

  for (aux=0;aux<gridNAux(grid);aux++) {
    gridSetAux(grid, avgNode, aux, 
	       1.0/2.0*( gridAux(grid, n0, aux) +
			 gridAux(grid, n1, aux) ) );
  }

  return grid;
}

Grid *gridSetAuxToAverageOfNodes3(Grid *grid, int avgNode,
				  int n0, int n1, int n2 )
{
  int aux;
  if ( !gridValidNode(grid,n0) || 
       !gridValidNode(grid,n1) || 
       !gridValidNode(grid,n2) ) return NULL;

  for (aux=0;aux<gridNAux(grid);aux++) {
    gridSetAux(grid, avgNode, aux, 
	       1.0/3.0*( gridAux(grid, n0, aux) +
			 gridAux(grid, n1, aux) +
			 gridAux(grid, n2, aux) ) );
  }

  return grid;
}

Grid *gridSetAuxToAverageOfNodes4(Grid *grid, int avgNode,
				  int n0, int n1, int n2, int n3 )
{
  int aux;
  if ( !gridValidNode(grid,n0) || 
       !gridValidNode(grid,n1) || 
       !gridValidNode(grid,n2) || 
       !gridValidNode(grid,n3) ) return NULL;

  for (aux=0;aux<gridNAux(grid);aux++) {
    gridSetAux(grid, avgNode, aux, 
	       1.0/4.0*( gridAux(grid, n0, aux) +
			 gridAux(grid, n1, aux) +
			 gridAux(grid, n2, aux) +
			 gridAux(grid, n3, aux) ) );
  }

  return grid;
}

Grid *gridSetPartId(Grid *grid, int partId )
{
  grid->partId = partId;
  return grid;
}

int gridGlobalNNode(Grid *grid)
{
  return grid->globalNNode;
} 

Grid *gridSetGlobalNNode(Grid *grid, int nglobal )
{
  grid->globalNNode = nglobal;
  return grid;
}

int gridGlobalNCell(Grid *grid)
{
  return grid->globalNCell;
}

Grid *gridSetGlobalNCell(Grid *grid, int nglobal )
{
  grid->globalNCell = nglobal;
  return grid;
}

int gridNUnusedNodeGlobal(Grid *grid )
{
  return grid->nUnusedNodeGlobal;
}

int gridNUnusedCellGlobal(Grid *grid )
{
  return grid->nUnusedCellGlobal;
}

Grid *gridGetUnusedNodeGlobal(Grid *grid, int *unused )
{
  int i;
  for (i=0;i<grid->nUnusedNodeGlobal;i++) 
    unused[i] = grid->unusedNodeGlobal[i];
  return grid;
}

Grid *gridGetUnusedCellGlobal(Grid *grid, int *unused )
{
  int i;
  for (i=0;i<grid->nUnusedCellGlobal;i++) 
    unused[i] = grid->unusedCellGlobal[i];
  return grid;
}

Grid *gridJoinUnusedNodeGlobal(Grid *grid, int global )
{
  int insertpoint, index;

  if (EMPTY == global) return grid;

  if (NULL == grid->unusedNodeGlobal) {
    grid->maxUnusedNodeGlobal = 500;
    grid->unusedNodeGlobal =
                         (int *)malloc(grid->maxUnusedNodeGlobal * sizeof(int));
  }
  if ((grid->nUnusedNodeGlobal+1) >= grid->maxUnusedNodeGlobal) {
    grid->maxUnusedNodeGlobal += 500;
    grid->unusedNodeGlobal =
                         (int *)realloc( grid->unusedNodeGlobal, 
				      grid->maxUnusedNodeGlobal*sizeof(int));
  }
  
  insertpoint = 0;
  if (grid->nUnusedNodeGlobal > 0) {
    for (index=grid->nUnusedNodeGlobal-1; index>=0; index--) {
      if (grid->unusedNodeGlobal[index] < global) {
	insertpoint = index+1;
	break;
      }
    }
    if ( grid->nUnusedNodeGlobal != insertpoint &&
	 grid->unusedNodeGlobal[insertpoint] == global) return grid;
    for(index=grid->nUnusedNodeGlobal;index>insertpoint;index--)
      grid->unusedNodeGlobal[index] = grid->unusedNodeGlobal[index-1];
  }

  grid->unusedNodeGlobal[insertpoint] = global;
  grid->nUnusedNodeGlobal++;
  
  return grid;
}

Grid *gridJoinUnusedCellGlobal(Grid *grid, int global )
{
  int insertpoint, index;

  if (EMPTY == global) return grid;

  if (NULL == grid->unusedCellGlobal) {
    grid->maxUnusedCellGlobal = 500;
    grid->unusedCellGlobal =
                         (int *)malloc(grid->maxUnusedCellGlobal * sizeof(int));
  }
  if ((grid->nUnusedCellGlobal+1) >= grid->maxUnusedCellGlobal) {
    grid->maxUnusedCellGlobal += 500;
    grid->unusedCellGlobal =
                         (int *)realloc( grid->unusedCellGlobal, 
				      grid->maxUnusedCellGlobal*sizeof(int));
  }
  
  insertpoint = 0;
  if (grid->nUnusedCellGlobal > 0) {
    for (index=grid->nUnusedCellGlobal-1; index>=0; index--) {
      if (grid->unusedCellGlobal[index] < global) {
	insertpoint = index+1;
	break;
      }
    }
    if ( grid->nUnusedCellGlobal != insertpoint &&
	 grid->unusedCellGlobal[insertpoint] == global) return grid;
    for(index=grid->nUnusedCellGlobal;index>insertpoint;index--)
      grid->unusedCellGlobal[index] = grid->unusedCellGlobal[index-1];
  }

  grid->unusedCellGlobal[insertpoint] = global;
  grid->nUnusedCellGlobal++;
  
  return grid;
}

Grid *gridEliminateUnusedNodeGlobal(Grid *grid )
{
  int node;
  int sort, offset;

  if ( 0 == gridNUnusedNodeGlobal(grid) ) return grid;

  if ( (NULL == grid->sortedLocal) || (grid->nnode != grid->nsorted) ) {
    if ( grid != gridCreateSortedGlobal(grid ) ) {
      printf("%s: %d: gridEliminateUnusedNodeGlobal: gridCreateSortedGlobal NULL.",
	     __FILE__,__LINE__);
      return NULL;
    }
  }

  /* for safety, see if other proc has already deleted this global */
  for ( offset = 0 ; offset < grid->nUnusedNodeGlobal ; offset++ ) {
    node = gridGlobal2Local(grid,grid->unusedNodeGlobal[offset]);
    if (EMPTY!=node) {
      printf("WARNING: %s: %d: Something is hosed. Found an eliminated node.%d\n",
	     __FILE__,__LINE__,node);      
      gridRemoveNodeWithOutGlobal(grid,node);
    }
  }

  offset = 0;
  for (sort=0;sort<grid->nsorted;sort++) {
    while ( (offset < grid->nUnusedNodeGlobal ) &&
	    (grid->unusedNodeGlobal[offset] < grid->sortedGlobal[sort] ) ) {
      offset++;
    }
    if ( offset < grid->nUnusedNodeGlobal &&
	 grid->unusedNodeGlobal[offset]==grid->sortedGlobal[sort]) {
      printf("ERROR: %s: %d: Global Node %d exists in sortedGlobal.%d\n",
	     __FILE__,__LINE__,grid->unusedNodeGlobal[offset],gridPartId(grid));
    }
    node = grid->sortedLocal[sort];
    grid->nodeGlobal[node] -= offset;
    grid->sortedGlobal[sort] -= offset;
  }
  grid->globalNNode -= grid->nUnusedNodeGlobal;
  grid->nUnusedNodeGlobal = 0;

  return grid;
}

Grid *gridEliminateUnusedCellGlobal(Grid *grid )
{
  int ncell, cell;
  int *pack;
  int *sortedGlobal, *sortedLocal;
  int sort, offset;

  if ( 0 == gridNUnusedCellGlobal(grid) ) return grid;

  pack         = (int *)malloc(grid->maxcell * sizeof(int));
  sortedGlobal = (int *)malloc(grid->maxcell * sizeof(int));
  sortedLocal  = (int *)malloc(grid->maxcell * sizeof(int));
  
  ncell = 0;
  for (cell=0;cell<grid->maxcell;cell++)
    if (gridCellValid(grid,cell)) {
      sortedGlobal[ncell] = grid->cellGlobal[cell];
      pack[ncell] = cell;
      ncell++;
    }

  if (ncell != grid->ncell)
    printf("%s: %d: gridEliminateUnusedCellGlobal: ncell error %d %d.",
	   __FILE__,__LINE__,ncell,grid->ncell);

  sortHeap(ncell,sortedGlobal,sortedLocal);

  offset = 0;
  for (sort=0;sort<ncell;sort++) {
    cell = pack[sortedLocal[sort]];
    while ( (offset < grid->nUnusedCellGlobal ) &&
	    (grid->unusedCellGlobal[offset] < grid->cellGlobal[cell] ) ) {
      offset++;
    }
    grid->cellGlobal[cell] -= offset;
  }
  grid->globalNCell -= grid->nUnusedCellGlobal;
  grid->nUnusedCellGlobal = 0;

  free(pack);
  free(sortedGlobal);
  free(sortedLocal);

  return grid;
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

Grid *gridRemoveGemAndQueue(Grid *grid, Queue *queue) {
  int i;
  for ( i = 0 ; i < grid->ngem ; i++ ) 
    gridRemoveCellAndQueue( grid, queue, grid->gem[i] );
  return grid;
}

int gridCellDegree(Grid *grid, int id)
{
  return adjDegree(grid->cellAdj, id);
}

int gridCellGlobal(Grid *grid, int cell )
{
  if ( !gridCellValid(grid, cell) ) return EMPTY;
  if (NULL == grid->cellGlobal) return EMPTY;
  return grid->cellGlobal[cell];
}

Grid *gridSetCellGlobal(Grid *grid, int cell, int global )
{
  if ( !gridCellValid(grid, cell) ) return NULL;
  if (NULL == grid->cellGlobal) 
    grid->cellGlobal = (int *)malloc(grid->maxcell*sizeof(int));
  grid->cellGlobal[cell] = global;
  return grid;
}

Grid *gridGlobalShiftCell(Grid *grid, int oldncellg, int newncellg, 
			  int celloffset )
{
  int cell;
  gridSetGlobalNCell(grid,newncellg);
  if (NULL == grid->cellGlobal) return NULL;  
  for (cell=0;cell<grid->maxcell;cell++)
    if ( gridCellValid(grid,cell) && (grid->cellGlobal[cell] >= oldncellg) ) 
      grid->cellGlobal[cell] += celloffset;
  for (cell=0;cell<grid->nUnusedCellGlobal;cell++)
    if ( grid->unusedCellGlobal[cell] >= oldncellg )  
      grid->unusedCellGlobal[cell] += celloffset;
  return grid;
}

int gridGetNextCellGlobal(Grid *grid)
{
  int global, i;

  if (NULL == grid->cellGlobal) {
    global = EMPTY;
  } else {
    if (grid->nUnusedCellGlobal > 0) {
      global = grid->unusedCellGlobal[0];
      for (i=1;i<grid->nUnusedCellGlobal;i++)
	grid->unusedCellGlobal[i-1]=grid->unusedCellGlobal[i];
      grid->nUnusedCellGlobal--;
    }else{
      global = grid->globalNCell;
      grid->globalNCell++;
    }
  }

  return global;
}

int gridAddCell(Grid *grid, int n0, int n1, int n2, int n3)
{
  int global;
  global = gridGetNextCellGlobal(grid);
  return gridAddCellWithGlobal(grid,n0,n1,n2,n3,global);
}

int gridAddCellAndQueue(Grid *grid, Queue *queue, 
			int n0, int n1, int n2, int n3)
{
  int global;
  int inode;
  int nodes[4], globalnodes[4], nodeParts[4];
  double xyz[1000];
  int dim, aux;
  
  global = gridGetNextCellGlobal(grid);

  nodes[0] = n0; nodes[1] = n1; nodes[2] = n2; nodes[3] = n3;

  if ( NULL != queue && gridCellHasGhostNode(grid, nodes) ) {
    dim = 3 + 6 + gridNAux(grid);
    if (dim>250) printf( "ERROR: %s: %d: undersized static xyz.\n", 
			 __FILE__, __LINE__);
    for ( inode = 0 ; inode < 4 ; inode++ ) {
      globalnodes[inode] = gridNodeGlobal(grid,nodes[inode]);
      nodeParts[inode] = gridNodePart(grid,nodes[inode]);
      gridNodeXYZ(grid,nodes[inode],&xyz[dim*inode]);
      gridMap(grid,nodes[inode],&xyz[3+dim*inode]);
      for ( aux = 0 ; aux < gridNAux(grid) ; aux++ ) 
	xyz[aux+9+dim*inode] = gridAux(grid, nodes[inode], aux);
    }
    queueAddCell(queue,globalnodes,global,nodeParts,xyz);
  }
  
  if ( gridCellHasLocalNode(grid,nodes) )
    return gridAddCellWithGlobal(grid, n0, n1, n2, n3, global);

  return EMPTY;
}

int gridAddCellWithGlobal(Grid *grid, int n0, int n1, int n2, int n3, 
			  int global )
{
  int cell, origSize, chunkSize;
  int cellId;
  if ( grid->blankc2n == EMPTY ) {
    origSize = grid->maxcell;
    chunkSize = MAX(5000,origSize/6);
    grid->maxcell += chunkSize;
    grid->c2n = (int *)realloc( grid->c2n, 4 * grid->maxcell * sizeof(int) );
    for (cell=origSize;cell < grid->maxcell; cell++ ) {
      grid->c2n[0+4*cell] = EMPTY; 
      grid->c2n[1+4*cell] = cell+1; 
    }
    grid->c2n[1+4*(grid->maxcell-1)] = EMPTY; 
    grid->blankc2n = origSize;
    if (NULL != grid->cellGlobal) 
      grid->cellGlobal =
                   (int *)realloc(grid->cellGlobal,grid->maxcell * sizeof(int));
    if (NULL != grid->reallocFunc)
      (*grid->reallocFunc)( grid->reallocData, gridREALLOC_CELL, 
			    origSize, grid->maxcell);
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
  
  if (NULL != grid->cellGlobal) grid->cellGlobal[cellId] = global;

  return cellId;
}

Grid *gridRemoveCell(Grid *grid, int cellId )
{
  Grid *result;

  result = gridRemoveCellWithOutGlobal(grid, cellId );

  if (grid == result && NULL != grid->cellGlobal) 
    gridJoinUnusedCellGlobal(grid,grid->cellGlobal[cellId]);

  return result;
}

Grid *gridRemoveCellAndQueue(Grid *grid, Queue *queue, int cellId )
{
  int inode, *nodes, globalnodes[4], nodeParts[4];

  if ( !gridCellValid(grid, cellId) ) return NULL;

  nodes = &(grid->c2n[4*cellId]);
  if (NULL!=queue && gridCellHasGhostNode(grid, nodes)) {
    for ( inode = 0 ; inode < 4 ; inode++ ) { 
      globalnodes[inode] = gridNodeGlobal(grid, nodes[inode]);
      nodeParts[inode] = gridNodePart(grid, nodes[inode]);
    }
    queueRemoveCell(queue,globalnodes,nodeParts);
  }
  
  return gridRemoveCell( grid, cellId );
}

Grid *gridRemoveCellWithOutGlobal(Grid *grid, int cellId )
{
  if ( !gridCellValid(grid, cellId) ) return NULL;
  
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
  GridBool frozen;

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

GridBool gridCellValid(Grid *grid, int cellId )
{
  if ( cellId < 0 || cellId >= grid->maxcell ) return FALSE;
  return (EMPTY != grid->c2n[4*cellId]);
}

GridBool gridCellEdge(Grid *grid, int n0, int n1 )
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

GridBool gridCellFace(Grid *grid, int n0, int n1, int n2 )
{
  return (EMPTY != gridFindOtherCellWith3Nodes(grid, n0, n1, n2, EMPTY ) );
}

int gridFindOtherCellWith3Nodes(Grid *grid, int n0, int n1, int n2,
				int currentCell )
{
  AdjIterator it;
  int cell, nodes[4];

  if (n0==n1 || n0==n2 || n1==n2 ) return EMPTY;

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

int gridFindCell(Grid *grid, int *nodes )
{
  AdjIterator it;
  int cell, n[4];

  for ( it=adjFirst(grid->cellAdj,nodes[0]); adjValid(it); it=adjNext(it) ) {
    cell = adjItem(it);
    gridCell( grid, cell, n );
    if ( ( n[0]==nodes[1]||n[1]==nodes[1]||n[2]==nodes[1]||n[3]==nodes[1] ) &&
	 ( n[0]==nodes[2]||n[1]==nodes[2]||n[2]==nodes[2]||n[3]==nodes[2] ) &&
	 ( n[0]==nodes[3]||n[1]==nodes[3]||n[2]==nodes[3]||n[3]==nodes[3] ) )
      return cell; 
  }

  return EMPTY;
}

int gridFindEnclosingCell(Grid *grid, int starting_guess, 
			  double *target, double *bary )
{
  int current_cell, last_cell;
  int tries;

  if (EMPTY == starting_guess) return EMPTY;

  current_cell = starting_guess;
  for (tries=1;tries<=1000;tries++) {
    int nodes[4];
    double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
    int other_cell;
    double tol;

    if (grid != gridCell( grid, current_cell, nodes )) {
      printf("%s: %d: gridFindClosestBoundaryCell %s %d\n",
	     __FILE__,__LINE__,"invalid cell",current_cell);
      return EMPTY;
    }

    gridNodeXYZ(grid, nodes[0], xyz0);
    gridNodeXYZ(grid, nodes[1], xyz1);
    gridNodeXYZ(grid, nodes[2], xyz2);
    gridNodeXYZ(grid, nodes[3], xyz3);
    
    gridBarycentricCoordinate( xyz0, xyz1, xyz2, xyz3, target, bary );

    tol = -1.0e-13;
    if ( bary[0] >= tol && bary[1] >= tol &&
	 bary[2] >= tol && bary[3] >= tol ) {
      return current_cell;
    }
    
    last_cell = current_cell;
    current_cell = EMPTY;

    if  ( bary[0] < bary[1] && bary[0] < bary[2] && bary[0] < bary[3] ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[1], nodes[2], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[1] < bary[0] && bary[1] < bary[2] && bary[1] < bary[3] ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[2], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[2] < bary[1] && bary[2] < bary[1] && bary[2] < bary[3] ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[3] < bary[0] && bary[3] < bary[1] && bary[3] < bary[2] ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[2],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[0] < tol ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[1], nodes[2], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[1] < tol ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[2], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[2] < tol ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[3],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( bary[3] < tol ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[2],
					       last_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if (EMPTY == current_cell) {
      current_cell = gridFindClosestBoundaryCell(grid, last_cell, target, bary);
      return current_cell;
    }
  }

  return EMPTY;
}

int gridFindClosestBoundaryCell(Grid *grid, int starting_guess, 
				double *target, double *bary )
{
  int current_face;
  int current_cell;
  int faceId;
  int nodes[4];
  double trib[3];
  double projected_target[3];
  double xyz0[3], xyz1[3], xyz2[3], xyz3[3]; 

  if (grid != gridCell( grid, starting_guess, nodes )) {
    printf("%s: %d: gridFindClosestBoundaryCell bad starting_guess.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }

  current_face = EMPTY;

  if (EMPTY == current_face) {
    current_face = gridFindFace( grid, nodes[0], nodes[1], nodes[2] );
  }
  if (EMPTY == current_face) {
    current_face = gridFindFace( grid, nodes[0], nodes[1], nodes[3] );
  }
  if (EMPTY == current_face) {
    current_face = gridFindFace( grid, nodes[0], nodes[2], nodes[3] );
  }
  if (EMPTY == current_face) {
    current_face = gridFindFace( grid, nodes[1], nodes[2], nodes[3] );
  }

  if (EMPTY == current_face) {
    printf("%s: %d: gridFindClosestBoundaryCell starting_guess has no face.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }

  current_face = gridFindClosestBoundaryFace(grid, current_face, 
					     target, trib );
  if (EMPTY == current_face) {
    return EMPTY;
  }

  current_cell = gridFindCellWithFace( grid, current_face );

  if (EMPTY == current_cell) {
    printf("%s: %d: gridFindClosestBoundaryCell: can't find closest cell.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }

  projected_target[0] = target[0];
  projected_target[1] = target[1];
  projected_target[2] = target[2];
  gridFace( grid, current_face, nodes, &faceId );
  gridNodeXYZ(grid, nodes[0], xyz0);
  gridNodeXYZ(grid, nodes[1], xyz1);
  gridNodeXYZ(grid, nodes[2], xyz2);
  gridProjectToTriangle(projected_target, xyz0, xyz1, xyz2  );

  gridCell( grid, current_cell, nodes );
  gridNodeXYZ(grid, nodes[0], xyz0);
  gridNodeXYZ(grid, nodes[1], xyz1);
  gridNodeXYZ(grid, nodes[2], xyz2);
  gridNodeXYZ(grid, nodes[3], xyz3);

  gridBarycentricCoordinate( xyz0, xyz1, xyz2, xyz3, projected_target, bary );

  return current_cell;
}

int gridFindClosestBoundaryFace(Grid *grid, int starting_guess, 
				double *target, double *trib )
{
  int current_face;
  int tries;

  if (EMPTY == starting_guess) return EMPTY;

  current_face = starting_guess;
  for (tries=1;tries<=100;tries++) {
    int nodes[3];
    int faceId;
    double xyz0[3], xyz1[3], xyz2[3];
    double projected_target[3];
    int last_face, other_face;
    double tol;

    if (grid != gridFace( grid, current_face, nodes, &faceId )) {
      printf("%s: %d: gridFindClosestBoundaryFace %s\n",
	     __FILE__,__LINE__," invalid face");
      return EMPTY;
    }

    gridNodeXYZ(grid, nodes[0], xyz0);
    gridNodeXYZ(grid, nodes[1], xyz1);
    gridNodeXYZ(grid, nodes[2], xyz2);

    projected_target[0] = target[0];
    projected_target[1] = target[1];
    projected_target[2] = target[2];

    gridProjectToTriangle(projected_target, xyz0, xyz1, xyz2  );
    
    gridBarycentricCoordinateTri( xyz0, xyz1, xyz2, projected_target, trib );

    tol = -1.0e-13;
    if ( trib[0] >= tol && trib[1] >= tol && trib[2] >= tol ) {
      return current_face;
    }
    
    last_face = current_face;
    current_face = EMPTY;

    if  ( trib[0] < tol ) {
      other_face = gridFindFaceWithNodesUnless(grid,
					       nodes[1], nodes[2],
					       last_face );
      if (EMPTY != other_face ) {
	current_face = other_face;
	continue;
      }
    }

    if  ( trib[1] < tol ) {
      other_face = gridFindFaceWithNodesUnless(grid,
					       nodes[0], nodes[2],
					       last_face );
      if (EMPTY != other_face ) {
	current_face = other_face;
	continue;
      }
    }

    if  ( trib[2] < tol ) {
      other_face = gridFindFaceWithNodesUnless(grid,
					       nodes[0], nodes[1],
					       last_face );
      if (EMPTY != other_face ) {
	current_face = other_face;
	continue;
      }
    }

    if (EMPTY == current_face) {
      printf("%s: %d: gridFindClosestBoundaryFace %s\n",
	     __FILE__,__LINE__," cant find next face for search");
      return EMPTY;
    }
  }

  return EMPTY;
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

static void gridAddConnToCell2Conn(Grid *grid, int conn, int node0, int node1)
{
  int iconn;
  int cell, nodes[4];
  int conn2node0[6] = {0, 0, 0, 1, 1, 2};
  int conn2node1[6] = {1, 2, 3, 2, 3, 3};
  int local0, local1;
  AdjIterator it;

  for ( it = adjFirst(gridCellAdj(grid),node0); 
	adjValid(it); 
	it = adjNext(it) ) {
    cell = adjItem(it);
    gridCell( grid, cell, nodes );
    for(iconn=0;iconn<6;iconn++) {
      local0 = nodes[conn2node0[iconn]];
      local1 = nodes[conn2node1[iconn]];
      if ( MIN(local0, local1) == MIN(node0,node1) &&
	   MAX(local0, local1) == MAX(node0,node1) ) {
	grid->cell2conn[iconn+6*cell] = conn;
      } 
    }
  }
}

Grid *gridCreateConn(Grid *grid) 
{
  int cell, conn;
  int nodes[4];
  int node0, node1;
  int conn2node0[6] = {0, 0, 0, 1, 1, 2};
  int conn2node1[6] = {1, 2, 3, 2, 3, 3};

  if (NULL != grid->cell2conn || NULL != grid->conn2node) return NULL; 

  grid->cell2conn = (int*)malloc(6*gridMaxCell(grid)*sizeof(int));
  for(cell=0;cell<6*gridMaxCell(grid);cell++) grid->cell2conn[cell] = EMPTY;
  grid->nconn = 0;
  for(cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid == gridCell(grid,cell,nodes)) {
      for(conn=0;conn<6;conn++) {
	if ( EMPTY == grid->cell2conn[conn+6*cell] ) {
	  grid->cell2conn[conn+6*cell] = grid->nconn;
	  gridAddConnToCell2Conn(grid, grid->nconn,
				 nodes[conn2node0[conn]], 
				 nodes[conn2node1[conn]]);
	  grid->nconn++;
	}
      }
    }
  } 

  grid->conn2node = (int*)malloc(2*gridNConn(grid)*sizeof(int));
  for(conn=0;conn<2*gridNConn(grid);conn++) grid->conn2node[conn] = EMPTY;

  for(cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid == gridCell(grid,cell,nodes)) {
      for(conn=0;conn<6;conn++) {
	if (EMPTY!=grid->cell2conn[conn+6*cell]) {
	  node0 = MIN(nodes[conn2node0[conn]],nodes[conn2node1[conn]]);
	  node1 = MAX(nodes[conn2node0[conn]],nodes[conn2node1[conn]]);
	  grid->conn2node[0+2*grid->cell2conn[conn+6*cell]] = node0;
	  grid->conn2node[1+2*grid->cell2conn[conn+6*cell]] = node1;
	}else{
	  printf("ERROR: %s: %d: cell2conn EMPTY.\n",__FILE__,__LINE__);
	}
      }
    }
  }

  return grid;
}

int gridCell2Conn(Grid *grid, int cell, int index )
{
  if (!gridCellValid(grid, cell)) return EMPTY;
  if (index<0||index>5) return EMPTY;
  if (NULL==grid->cell2conn) gridCreateConn(grid);
  return grid->cell2conn[index+6*cell];
}

Grid *gridConn2Node(Grid *grid, int conn, int *nodes )
{
  if (NULL==grid->conn2node) gridCreateConn(grid);
  if (conn<0||conn>=gridNConn(grid)) return NULL;
  nodes[0] = grid->conn2node[0+2*conn];
  nodes[1] = grid->conn2node[1+2*conn];
  return grid;
}

int gridFindConn(Grid *grid, int node0, int node1 )
{
  AdjIterator it;
  int cell, index, conn, conn0, conn1;
  if (NULL==grid->conn2node) return EMPTY;
  if (!gridValidNode(grid,node0) || !gridValidNode(grid,node1) ) return EMPTY;
  for ( it = adjFirst(gridCellAdj(grid),node0);
	adjValid(it);
	it = adjNext(it) ){
    cell = adjItem(it);
    for ( index=0 ; index < 6 ; index++ ) {
      conn = grid->cell2conn[index+6*cell];
      conn0 = grid->conn2node[0+2*conn];
      conn1 = grid->conn2node[1+2*conn];
      if ( node0==conn0 && node1==conn1 ) return conn;
      if ( node0==conn1 && node1==conn0 ) return conn;
    }
  }
  return EMPTY;
}

Grid *gridEraseConn(Grid* grid)
{
  grid->nconn=0;
  if (NULL != grid->cell2conn) free(grid->cell2conn);
  grid->cell2conn = NULL;
  if (NULL != grid->conn2node) free(grid->conn2node);
  grid->conn2node = NULL;
  return grid;
}

Grid *gridConstrainSurfaceNode(Grid *grid)
{
  grid->constrain_surface_node = TRUE;
  return grid;
}

Grid *gridUnconstrainSurfaceNode(Grid *grid)
{
  grid->constrain_surface_node = FALSE;
  return grid;
}

GridBool gridSurfaceNodeConstrained(Grid *grid)
{
  return grid->constrain_surface_node;
}

int gridAddFace(Grid *grid, int n0, int n1, int n2, int faceId )
{
 return gridAddFaceUV(grid, 
		      n0, DBL_MAX, DBL_MAX,
		      n1, DBL_MAX, DBL_MAX,
		      n2, DBL_MAX, DBL_MAX,
		      faceId );
}

int gridAddFaceUVAndQueue(Grid *grid, Queue *queue,
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId )
{
  int g0, g1, g2;
  int p0, p1, p2;

  if ( NULL != queue && gridFaceHasGhostNode(grid,n0,n1,n2) ) {
    g0 = gridNodeGlobal(grid,n0);
    g1 = gridNodeGlobal(grid,n1);
    g2 = gridNodeGlobal(grid,n2);
    p0 = gridNodePart(grid,n0);
    p1 = gridNodePart(grid,n1);
    p2 = gridNodePart(grid,n2);
    queueAddFaceScalar(queue, 
		       g0, p0, u0, v0,
		       g1, p1, u1, v1,
		       g2, p2, u2, v2, faceId );   
  }
  
  if ( gridFaceHasLocalNode(grid,n0,n1,n2) )
    return gridAddFaceUV( grid, 
			  n0, u0, v0,
			  n1, u1, v1,
			  n2, u2, v2, faceId );
  
  return EMPTY;
}

int gridAddFaceUV(Grid *grid, 
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId )
{
  int face, i, origSize, chunkSize;
  if ( grid->blankf2n == EMPTY ) {
    origSize = grid->maxface;
    chunkSize = MAX(5000,origSize/6);
    grid->maxface += chunkSize;
    grid->f2n    =
                (int    *)realloc(grid->f2n,    3*grid->maxface*sizeof(int));
    grid->faceU  =
                (double *)realloc(grid->faceU,  3*grid->maxface*sizeof(double));
    grid->faceV  =
                (double *)realloc(grid->faceV,  3*grid->maxface*sizeof(double));
    grid->faceId =
                (int    *)realloc(grid->faceId, 1*grid->maxface*sizeof(int));
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
    if (NULL != grid->reallocFunc)
      (*grid->reallocFunc)( grid->reallocData, gridREALLOC_FACE, 
			    origSize, grid->maxface);
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

Grid *gridRemoveFaceAndQueue(Grid *grid, Queue *queue, int face )
{
  int inode, globalnodes[3], nodeParts[3];
  if (face >= grid->maxface || face < 0) return NULL;
  if (EMPTY == grid->f2n[3*face]) return NULL;

  if (NULL!=queue && gridFaceHasGhostNode(grid,
					  grid->f2n[0+3*face],
					  grid->f2n[1+3*face],
					  grid->f2n[2+3*face])) {
    for ( inode = 0 ; inode < 3 ; inode++ ) {
      globalnodes[inode] = gridNodeGlobal(grid, grid->f2n[inode+3*face]);
      nodeParts[inode] = gridNodePart(grid, grid->f2n[inode+3*face]);
    }
    queueRemoveFace(queue,globalnodes,nodeParts);
  }
  
  return gridRemoveFace(grid, face );
}

int gridFindFace(Grid *grid, int n0, int n1, int n2 )
{
  AdjIterator it0, it1, it2;
  Adj *adj=grid->faceAdj;

  if ( n0 == n1 || n0 == n2 || n1 == n2 ) return EMPTY;

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

int gridFindFaceWithNodesUnless(Grid *grid, int n0, int n1, int unless_face )
{
  AdjIterator it0, it1;
  Adj *adj=grid->faceAdj;

  if ( n0 == n1 ) return EMPTY;

  for ( it0 = adjFirst(adj,n0); adjValid(it0); it0 = adjNext(it0) )
    for ( it1 = adjFirst(adj,n1); adjValid(it1); it1 = adjNext(it1) )
      if ( adjItem(it0) == adjItem(it1) )
	if ( adjItem(it0) != unless_face ) return adjItem(it0);

  return EMPTY;
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

GridBool gridReconnectionOfAllFacesOK(Grid *grid, int oldNode, int newNode )
{
  AdjIterator it;
  int face, i, nodes[3], faceId;
  if (oldNode < 0 || oldNode >= grid->maxnode ) return FALSE;
  if (newNode < 0 || newNode >= grid->maxnode ) return FALSE;
  if (newNode == oldNode) return TRUE;

  if (!gridGeometryFace(grid, oldNode)) return TRUE;
  if (!gridGeometryFace(grid, newNode)) return TRUE;

  for (it = adjFirst(grid->faceAdj,oldNode); adjValid(it); it = adjNext(it)) {
    face = adjItem(it);
    gridFace(grid, face, nodes, &faceId );
    for (i=0;i<3;i++){
      if (oldNode == nodes[i]) nodes[i]=newNode;
    }
    /* verify both ways in case gridFindFace ever checks face orientation */
    if (EMPTY != gridFindFace(grid, nodes[0], nodes[1], nodes[2]))return FALSE;
    if (EMPTY != gridFindFace(grid, nodes[1], nodes[0], nodes[2]))return FALSE;
  }

  return TRUE;
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
  GridBool found;
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

int gridNodeFaceIdDegree(Grid *grid, int node)
{
  int id[MAXFACEIDDEG];
  int face, faceId, search, ids;
  GridBool found;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return EMPTY;

  ids = 0;
  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    faceId = grid->faceId[face];
    found = FALSE;
    for (search=0;!found && (search<ids);search++) {
      found = ( id[search] == faceId );
    }
    if (!found) {
      if (ids >= MAXFACEIDDEG) {
	printf("%s: %d: need more MAXFACEIDDEG.\n",__FILE__,__LINE__);
	return EMPTY;
      }
      id[ids] = faceId;
      ids++;
    }
  }

  return ids;
}

Grid *gridNodeFaceId(Grid *grid, int node, int maxId, int *ids_arg, int *id )
{
  int face, faceId, search, insertpoint, index, ids;
  GridBool found;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return NULL;

  *ids_arg = 0;
  ids = 0;
  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    faceId = grid->faceId[face];
    found = FALSE;
    for (search=0;!found && (search<ids);search++) {
      found = ( id[search] == faceId );
    }
    if (!found) {
      if (ids >= maxId) {
	printf("%s: %d: need more maxId.\n",__FILE__,__LINE__);
	return NULL;
      }
      insertpoint = 0;
      for (index=ids-1; index>=0; index--) {
	if (id[index] < faceId) {
	  insertpoint = index+1;
	  break;
	}
      }
      for(index=ids;index>insertpoint;index--)
	id[index] = id[index-1];
      ids++;
      id[insertpoint] = faceId;
    }
  }

  *ids_arg = ids;
  return grid;
}

int gridTrianglesOnFaceId(Grid *grid, int id )
{
  int face, count;
  int nodes[3];
  int faceId;

  count = 0;
  for(face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace( grid, face, nodes, &faceId )) {
      if (faceId==id) count++;
    }
  }
  return count;
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
  GridBool found;
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

int gridNodeEdgeIdDegree(Grid *grid, int node)
{
  int id[MAXEDGEIDDEG];
  int edge, edgeId, search, ids;
  GridBool found;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return EMPTY;

  ids = 0;
  for ( it = adjFirst(grid->edgeAdj,node); adjValid(it); it = adjNext(it) ){
    edge = adjItem(it);
    edgeId = grid->edgeId[edge];
    found = FALSE;
    for (search=0;!found && (search<ids);search++) {
      found = ( id[search] == edgeId );
    }
    if (!found) {
      if (ids >= MAXEDGEIDDEG) {
	printf("%s: %d: need more MAXEDGEIDDEG.\n",__FILE__,__LINE__);
	return EMPTY;
      }
      id[ids] = edgeId;
      ids++;
    }
  }

  return ids;
}

Grid *gridNodeEdgeId(Grid *grid, int node, int maxId, int *ids_arg, int *id )
{
  int edge, edgeId, search, insertpoint, index, ids;
  GridBool found;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return NULL;

  *ids_arg = 0;
  ids = 0;
  for ( it = adjFirst(grid->edgeAdj,node); adjValid(it); it = adjNext(it) ){
    edge = adjItem(it);
    edgeId = grid->edgeId[edge];
    found = FALSE;
    for (search=0;!found && (search<ids);search++) {
      found = ( id[search] == edgeId );
    }
    if (!found) {
      if (ids >= maxId) {
	printf("%s: %d: need more maxId.\n",__FILE__,__LINE__);
	return NULL;
      }
      insertpoint = 0;
      for (index=ids-1; index>=0; index--) {
	if (id[index] < edgeId) {
	  insertpoint = index+1;
	  break;
	}
      }
      for(index=ids;index>insertpoint;index--)
	id[index] = id[index-1];
      ids++;
      id[insertpoint] = edgeId;
    }
  }

  *ids_arg = ids;
  return grid;
}

int gridAddEdge(Grid *grid, int n0, int n1, 
		int edgeId, double t0, double t1 )
{
  int edge, i, origSize, chunkSize;
  if ( grid->blanke2n == EMPTY ) {
    origSize = grid->maxedge;
    chunkSize = MAX(5000,origSize/6);
    grid->maxedge += chunkSize;
    grid->e2n    =
                (int    *)realloc(grid->e2n,    2*grid->maxedge*sizeof(int));
    grid->edgeId =
                (int    *)realloc(grid->edgeId, 1*grid->maxedge*sizeof(int));
    grid->edgeT  =
                (double *)realloc(grid->edgeT,  2*grid->maxedge*sizeof(double));
    for (i=origSize;i < grid->maxedge; i++ ) {
      grid->e2n[0+2*i] = EMPTY; 
      grid->e2n[1+2*i] = i+1; 
      grid->edgeId[i] = EMPTY; 
      grid->edgeT[0+2*i] = DBL_MAX; 
      grid->edgeT[1+2*i] = DBL_MAX; 
    }
    grid->e2n[1+2*(grid->maxedge-1)] = EMPTY; 
    grid->blanke2n = origSize;
    if (NULL != grid->reallocFunc)
      (*grid->reallocFunc)( grid->reallocData, gridREALLOC_EDGE, 
			    origSize, grid->maxedge);
  }
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

int gridAddEdgeAndQueue(Grid *grid, Queue *queue, int n0, int n1, 
			int edgeId, double t0, double t1 )
{
  int globalNodes[2];
  int nodeParts[2];
  double ts[2];

  if ( NULL != queue && gridEdgeHasGhostNode(grid,n0,n1) ) {
    globalNodes[0] = gridNodeGlobal(grid,n0);
    globalNodes[1] = gridNodeGlobal(grid,n1);
    nodeParts[0] = gridNodePart(grid,n0);
    nodeParts[1] = gridNodePart(grid,n1);
    ts[0] = t0;
    ts[1] = t1;
    queueAddEdge(queue, globalNodes,  edgeId, nodeParts, ts);
  }
  
  if ( gridEdgeHasLocalNode(grid,n0,n1) )
    return gridAddEdge( grid, n0, n1, edgeId, t0, t1);
  
  return EMPTY;
}

int gridAddEdgeInGlobal(Grid *grid, int g0, int g1, 
			int edgeId, double t0, double t1 )
{
  int n0, n1;
  n0 = gridGlobal2Local(grid, g0);
  n1 = gridGlobal2Local(grid, g1);
  if ( EMPTY == n0 || EMPTY == n1 ) return EMPTY;
  if ( !gridNodeLocal( grid, n0 ) && !gridNodeLocal( grid, n1 ) ) return EMPTY;
  return gridAddEdge(grid, n0, n1, edgeId, t0, t1);
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

Grid *gridRemoveEdgeAndQueue(Grid *grid, Queue *queue, int edge )
{
  int n0, n1, globalnodes[2], nodeParts[2];
  if (edge >= grid->maxedge || edge < 0) return NULL;
  if (EMPTY == grid->e2n[2*edge]) return NULL;
  
  n0 = grid->e2n[0+2*edge];
  n1 = grid->e2n[1+2*edge];

  if (NULL!=queue && gridEdgeHasGhostNode(grid,n0,n1)) {
    globalnodes[0] = gridNodeGlobal(grid, n0);
    globalnodes[1] = gridNodeGlobal(grid, n1);
    nodeParts[0] = gridNodePart(grid, n0);
    nodeParts[1] = gridNodePart(grid, n1);
    queueRemoveEdge(queue,globalnodes,nodeParts);
  }
  
  return gridRemoveEdge(grid, edge );
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
  GridBool found;

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
  GridBool found;

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
  GridBool found;

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

GridBool gridNodeFrozen( Grid *grid, int node )
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

GridBool gridGemIsAllLocal(Grid *grid)
{
  int gem, cell;
  for ( gem = 0 ; gem < grid->ngem ; gem++ ) {
    cell = grid->gem[gem];
    if ( gridCellHasGhostNode(grid, &(grid->c2n[4*cell])) ) return FALSE;
  }
  return TRUE;
}

GridBool gridNodeNearGhost(Grid *grid, int node )
{
  AdjIterator it;

  for ( it = adjFirst(grid->cellAdj,node); 
	adjValid(it); 
	it = adjNext(it)) {
    if (gridCellHasGhostNode(grid, &(grid->c2n[4*adjItem(it)]))) return TRUE;
  }

  return FALSE;
}

Grid *gridFaceOppositeCellNode(Grid *grid, int *nodes, int node, int *face )
{
  if (nodes[0] == node) {
    face[0] = nodes[1];
    face[1] = nodes[3];
    face[2] = nodes[2];
    return grid;
  }

  if (nodes[1] == node) {
    face[0] = nodes[0];
    face[1] = nodes[2];
    face[2] = nodes[3];
    return grid;
  }

  if (nodes[2] == node) {
    face[0] = nodes[0];
    face[1] = nodes[3];
    face[2] = nodes[1];
    return grid;
  }

  if (nodes[3] == node) {
    face[0] = nodes[0];
    face[1] = nodes[1];
    face[2] = nodes[2];
    return grid;
  }
  
  return NULL;
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
  int igem, iequ, nodes[4];
  GridBool gap, found;
  grid->nequ = 0;

  if ( NULL == gridMakeGem( grid, n0, n1) ) return NULL;

  if ( grid->ngem == (MAXDEG-1) ) return NULL;
  if ( grid->ngem == 0 ) return grid;

  nodes[0] = n0;
  nodes[1] = n1;

  gap = FALSE;

  gridOrient( grid, &grid->c2n[4*grid->gem[0]], nodes );

  grid->equ[0]= nodes[3];

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

GridBool gridContinuousEquator(Grid *grid)
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

int gridGetNextNodeGlobal(Grid *grid)
{
  int global, i;

  if (NULL == grid->nodeGlobal) {
    global = EMPTY;
  } else {
    if (grid->nUnusedNodeGlobal > 0) {
      global = grid->unusedNodeGlobal[0];
      for (i=1;i<grid->nUnusedNodeGlobal;i++)
	grid->unusedNodeGlobal[i-1]=grid->unusedNodeGlobal[i];
      grid->nUnusedNodeGlobal--;
    }else{
      global = grid->globalNNode;
      grid->globalNNode++;
    }
  }

  return global;
}

int gridAddNode(Grid *grid, double x, double y, double z )
{
  int global;
  global = gridGetNextNodeGlobal(grid);
  return gridAddNodeWithGlobal(grid,x,y,z,global);
}

int gridAddNodeWithGlobal(Grid *grid, double x, double y, double z, int global )
{
  int node, i, origSize, chunkSize;

  if (EMPTY == grid->blanknode) {
    origSize = grid->maxnode;
    chunkSize = MAX(5000,origSize/6);
    grid->maxnode += chunkSize;
    grid->xyz =
               (double *)realloc(grid->xyz, 3 * grid->maxnode * sizeof(double));
    for (i=origSize;i < grid->maxnode; i++ ) {
      grid->xyz[0+3*i] = DBL_MAX;
      grid->xyz[1+3*i] = (double)(i+1);
    }
    grid->xyz[1+3*(grid->maxnode-1)] = (double)(EMPTY);
    grid->blanknode = origSize;

    if (NULL != grid->map) grid->map =
      (double *)realloc(grid->map, grid->maxnode * 6 * sizeof(double));
    if ( grid->naux > 0 ) grid->aux = 
      (double *)realloc(grid->aux, grid->maxnode * grid->naux * sizeof(double));
    grid->frozen =
             (GridBool *)realloc(grid->frozen,grid->maxnode * sizeof(GridBool));

    grid->geomNode = (int *)realloc(grid->geomNode, grid->maxnode*sizeof(int));
    for (i=origSize;i < grid->maxnode; i++ ) grid->geomNode[i] = EMPTY;

    if (NULL != grid->child_reference) 
      grid->child_reference =
                   (int *)realloc(grid->child_reference,grid->maxnode * sizeof(int));
    if (NULL != grid->nodeGlobal) 
      grid->nodeGlobal =
                   (int *)realloc(grid->nodeGlobal,grid->maxnode * sizeof(int));
    if (NULL != grid->part) 
      grid->part = (int *)realloc(grid->part,grid->maxnode * sizeof(int));
    if (NULL != grid->sortedGlobal) 
      grid->sortedGlobal =
                   (int *)realloc(grid->sortedGlobal,grid->maxnode*sizeof(int));
    if (NULL != grid->sortedLocal) 
      grid->sortedLocal =
                    (int *)realloc(grid->sortedLocal,grid->maxnode*sizeof(int));

    adjRealloc(grid->cellAdj,grid->maxnode);
    adjRealloc(grid->faceAdj,grid->maxnode);
    adjRealloc(grid->edgeAdj,grid->maxnode);

    if (NULL != grid->prismDeg) {
      grid->prismDeg = (int *)realloc(grid->prismDeg,grid->maxnode*sizeof(int));
      for (i=origSize;i < grid->maxnode; i++ ) grid->prismDeg[i] = 0;
    }
    if (NULL != grid->reallocFunc)
      (*grid->reallocFunc)( grid->reallocData, gridREALLOC_NODE, 
			    origSize, grid->maxnode);
  }

  node = grid->blanknode;
  grid->blanknode = (int)grid->xyz[1+3*node];
  grid->nnode++;

  grid->xyz[0+3*node] = x;
  grid->xyz[1+3*node] = y;
  grid->xyz[2+3*node] = z;
  grid->frozen[node] = FALSE;
  if (NULL != grid->map) {
    grid->map[0+6*node] = 1.0;
    grid->map[1+6*node] = 0.0;
    grid->map[2+6*node] = 0.0;
    grid->map[3+6*node] = 1.0;
    grid->map[4+6*node] = 0.0;
    grid->map[5+6*node] = 1.0;
  }
  if (NULL != grid->child_reference) grid->child_reference[node] = EMPTY;
  if (0 <= global) gridSetNodeGlobal(grid, node, global );
  if (NULL != grid->part) grid->part[node] = gridPartId(grid);

  return node;
}

Grid *gridRemoveNode(Grid *grid, int node )
{
  Grid *result;
  int global;

  global = gridNodeGlobal(grid,node);

  result = gridRemoveNodeWithOutGlobal(grid, node );

  if (grid == result && global != EMPTY ) 
    gridJoinUnusedNodeGlobal(grid,grid->nodeGlobal[node]);

  return result;
}

Grid *gridRemoveNodeWithOutGlobal(Grid *grid, int node )
{
  int index, removepoint;
  if (!gridValidNode(grid,node)) return NULL;
  grid->nnode--;
  grid->xyz[0+3*node] = DBL_MAX;
  grid->xyz[1+3*node] = (double)grid->blanknode;
  grid->blanknode = node;
  if (NULL != grid->sortedLocal) {
    removepoint = 
      sortSearch(grid->nsorted,grid->sortedGlobal,grid->nodeGlobal[node]);
    if (EMPTY != removepoint) {
      grid->nsorted--;
      for(index=removepoint;index<grid->nsorted;index++)
	grid->sortedGlobal[index] = grid->sortedGlobal[index+1];
      for(index=removepoint;index<grid->nsorted;index++)
	grid->sortedLocal[index] = grid->sortedLocal[index+1];
    }
  }

  return grid;
}

Grid *gridTranslate(Grid *grid, double dx, double dy, double dz)
{
  int node;
  for (node=0;node<gridMaxNode(grid);node++) {
    if (gridValidNode(grid,node)) {
      grid->xyz[0+3*node] += dx;
      grid->xyz[1+3*node] += dy;
      grid->xyz[2+3*node] += dz;
    }
  }
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

int gridNodeGlobal(Grid *grid, int node )
{
  if (!gridValidNode(grid,node)) return EMPTY;
  if (NULL == grid->nodeGlobal) return EMPTY;
  return grid->nodeGlobal[node];
}

Grid *gridCreateSortedGlobal(Grid *grid ) 
{
  int local, nnode;
  int *pack;

  if (NULL == grid->nodeGlobal) return NULL;

  if (NULL != grid->sortedLocal) free(grid->sortedLocal);
  if (NULL != grid->sortedGlobal) free(grid->sortedGlobal);
  grid->sortedLocal  = (int *)malloc(grid->maxnode * sizeof(int));
  grid->sortedGlobal = (int *)malloc(grid->maxnode * sizeof(int));

  pack = (int *)malloc(grid->maxnode * sizeof(int));

  nnode = 0;
  for (local=0;local<grid->maxnode;local++)
    if (gridValidNode(grid,local)) {
      grid->sortedGlobal[nnode] = grid->nodeGlobal[local];
      pack[nnode] = local;
      nnode++;
    }

  if (nnode != grid->nnode) {
    printf("%s: %d: gridCreateSortedGlobal: nnode error %d %d.",
	   __FILE__,__LINE__,nnode,grid->nnode);
    return NULL;
  }

  grid->nsorted = grid->nnode;
  sortHeap(grid->nsorted,grid->sortedGlobal,grid->sortedLocal);

  for (local=0;local<grid->nsorted;local++) {
    grid->sortedLocal[local] =  pack[grid->sortedLocal[local]];
    grid->sortedGlobal[local] = grid->nodeGlobal[grid->sortedLocal[local]];
  }
  free(pack);

  return grid;
}

int gridGlobal2Local(Grid *grid, int global )
{
  int position;

  if (NULL == grid->nodeGlobal) return EMPTY;

  if (NULL == grid->sortedLocal) gridCreateSortedGlobal(grid);

  position = sortSearch(grid->nsorted,grid->sortedGlobal,global);

  if (EMPTY == position) return EMPTY;

  return grid->sortedLocal[position];
}

Grid *gridSetNodeGlobal(Grid *grid, int node, int global )
{
  int index, insertpoint;
  if (!gridValidNode(grid,node)) return NULL;
  if (NULL == grid->nodeGlobal) 
    grid->nodeGlobal = (int *)malloc(grid->maxnode*sizeof(int));
  grid->nodeGlobal[node] = global;

  if (NULL != grid->sortedLocal) {
    if (grid->nsorted >= grid->maxnode) {
      printf("%s: %d: Error: gridSetNodeGlobal called twice on a node.\n",
	     __FILE__,__LINE__);
      return NULL;
    }
    insertpoint = 0;
    for (index=grid->nsorted-1; index>=0; index--) {
      if (grid->sortedGlobal[index] < global) {
	insertpoint = index+1;
	break;
      }
    }
    for(index=grid->nsorted;index>insertpoint;index--)
      grid->sortedGlobal[index] = grid->sortedGlobal[index-1];
    for(index=grid->nsorted;index>insertpoint;index--)
      grid->sortedLocal[index] = grid->sortedLocal[index-1];
    grid->nsorted++;
    grid->sortedGlobal[insertpoint] = global;
    grid->sortedLocal[insertpoint] = node;
  }

  return grid;
}

Grid *gridGlobalShiftNode(Grid *grid, int oldnnodeg, int newnnodeg, 
			  int nodeoffset )
{
  int node;

  gridSetGlobalNNode(grid,newnnodeg);

  if (NULL == grid->nodeGlobal) return NULL;  

  for (node=0;node<grid->maxnode;node++)
    if ( gridValidNode(grid,node) && (grid->nodeGlobal[node] >= oldnnodeg) ) 
      grid->nodeGlobal[node] += nodeoffset;

  for (node=0;node<grid->nUnusedNodeGlobal;node++)
    if ( grid->unusedNodeGlobal[node] >= oldnnodeg )  
      grid->unusedNodeGlobal[node] += nodeoffset;

  if (NULL == grid->sortedGlobal) return grid;

  for (node=0;node<grid->nsorted;node++)
    if ( grid->sortedGlobal[node] >= oldnnodeg ) 
      grid->sortedGlobal[node] += nodeoffset;

  return grid;
}

Grid *gridRenumberGlobalNodes(Grid *grid, int nnode, int *n2o)
{
  int oldglobal, newglobal;
  int oldposition, newposition;
  int oldlocal, newlocal;
  int fix;

  if ( NULL == grid->sortedLocal ) gridCreateSortedGlobal(grid);

  for (newglobal=0;newglobal<nnode;newglobal++){
    oldglobal = n2o[newglobal];
    if (newglobal != oldglobal) {
      newposition = sortSearch(grid->nsorted,grid->sortedGlobal,newglobal);
      newlocal = EMPTY;
      if (EMPTY!=newposition) newlocal = grid->sortedLocal[newposition];

      oldposition = sortSearch(grid->nsorted,grid->sortedGlobal,oldglobal);
      oldlocal = EMPTY;
      if (EMPTY!=oldposition) oldlocal = grid->sortedLocal[oldposition];

      if ( EMPTY != oldlocal || EMPTY != newlocal ) {
	if ( EMPTY != oldlocal ) grid->nodeGlobal[oldlocal] = newglobal;
	if ( EMPTY != newlocal ) grid->nodeGlobal[newlocal] = oldglobal;

	if ( EMPTY != oldlocal && EMPTY != newlocal ) {
	  grid->sortedLocal[newposition] = oldlocal;
	  grid->sortedLocal[oldposition] = newlocal;
	}else{
	  gridCreateSortedGlobal(grid);
	}
      }

      /* fix n2o list */
      n2o[newglobal] = newglobal;
      for (fix=newglobal+1;fix<nnode;fix++) {
	if ( newglobal == n2o[fix] ) {
	  n2o[fix] = oldglobal;
	  break;
	}
      }
    }
  }

  return grid;
}

int gridNodePart(Grid *grid, int node )
{
  if (!gridValidNode(grid,node)) return EMPTY;
  if (NULL == grid->part) return EMPTY;
  return grid->part[node];
}

Grid *gridSetNodePart(Grid *grid, int node, int part )
{
  int i;
  if (!gridValidNode(grid,node)) return NULL;
  if (NULL == grid->part) {
    grid->part = (int *)malloc(gridMaxNode(grid)*sizeof(int));
    for (i=0;i<gridMaxNode(grid);i++) grid->part[i] = gridPartId(grid);
  }
  grid->part[node] = part;
  return grid;
}

GridBool gridNodeLocal(Grid *grid, int node )
{
  if (!gridValidNode(grid,node)) return FALSE;
  if (NULL == grid->part) return TRUE;
  return (grid->partId==grid->part[node]);
}

GridBool gridNodeGhost(Grid *grid, int node )
{
  if (!gridValidNode(grid,node)) return FALSE;
  if (NULL == grid->part) return FALSE;
  return (grid->partId!=grid->part[node]);
}

Grid *gridDeleteNodesNotUsed(Grid *grid){
  int node, maxnode;
  GridBool prism;

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
  int global, node;
  grid->nGeomNode = nGeomNode;
  for(node=0;node<grid->maxnode;node++) grid->geomNode[node]=EMPTY;
  if (NULL == grid->nodeGlobal) {
    for(node=0;node<grid->nGeomNode;node++) grid->geomNode[node]=node;
  }else{
    for(global=0;global<grid->nGeomNode;global++) {
      node = gridGlobal2Local(grid, global );
      if (node!=EMPTY) grid->geomNode[node]=global;
    }
  }
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
  grid->geomEdge = (int *)malloc(2*nGeomEdge*sizeof(int));
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
  GridBool found;

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

GridBool gridGeometryNode(Grid *grid, int node)
{
  if (node < 0 || node >= grid->maxnode ) return FALSE;
  return (EMPTY!=grid->geomNode[node]);
}

GridBool gridGeometryEdge(Grid *grid, int node)
{

  if ( adjValid( adjFirst(grid->edgeAdj,node) ) ) 
    return TRUE;
  
  return FALSE;
}

GridBool gridGeometryFace(Grid *grid, int node)
{

  if ( adjValid( adjFirst(grid->faceAdj,node) ) )
    return TRUE;
  
  return FALSE;
}

GridBool gridGeometryBetweenFace(Grid *grid, int node)
{
  AdjIterator it;
  int faceId, firstFaceId, nodes[3];

  faceId = EMPTY;
  firstFaceId = EMPTY;
  gridFace(grid,adjItem(adjFirst(grid->faceAdj,node)),nodes,&firstFaceId);
    
  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ) {
    gridFace(grid,adjItem(it),nodes,&faceId);
    if (faceId !=firstFaceId) return TRUE;
  }

  return FALSE;
}

int gridParentGeometry(Grid *grid, int node0, int node1)
{
  int edgeId;
  int faceId, nodes[3];
  AdjIterator it0, it1;
  Adj *adj=grid->faceAdj;

  edgeId = gridEdgeId(grid, node0, node1 );
  if (edgeId>0) return (-edgeId);

  for ( it0 = adjFirst(adj,node0); adjValid(it0); it0 = adjNext(it0) )
    for ( it1 = adjFirst(adj,node1); adjValid(it1); it1 = adjNext(it1) )
      if ( adjItem(it0) == adjItem(it1) ) {
	if (grid==gridFace(grid, adjItem(it0), nodes, &faceId)) return faceId;
      }

  return 0;
}

Grid *gridAddPrism(Grid *grid, int n0, int n1, int n2, int n3, int n4, int n5)
{
  int i;
  int origSize;

  if (grid->nprism >= grid->maxprism) {
    origSize = grid->maxprism;
    grid->maxprism += 5000;
    if (grid->prism == NULL) {
      grid->prism = (Prism *)malloc(grid->maxprism*sizeof(Prism));
      grid->prismDeg = (int *)malloc(grid->maxnode*sizeof(int));
      for(i=0;i<grid->maxnode;i++) grid->prismDeg[i]=0;
    }else{
      grid->prism = (Prism *)realloc(grid->prism,grid->maxprism*sizeof(Prism));
      if (NULL != grid->reallocFunc)
	(*grid->reallocFunc)( grid->reallocData, 
			      gridREALLOC_PRISM, 
			      origSize, grid->maxprism);
    }
  }

  grid->prism[grid->nprism].nodes[0] = n0;
  grid->prism[grid->nprism].nodes[1] = n1;
  grid->prism[grid->nprism].nodes[2] = n2;
  grid->prism[grid->nprism].nodes[3] = n3;
  grid->prism[grid->nprism].nodes[4] = n4;
  grid->prism[grid->nprism].nodes[5] = n5;

  for(i=0;i<6;i++) {
    if ( (0 > grid->prism[grid->nprism].nodes[i]) || 
	 (grid->maxnode <= grid->prism[grid->nprism].nodes[i]) ) {
      printf("%s: %d: Could not register new prism, invalid node.\n",
	     __FILE__,__LINE__);
    } else {
      grid->prismDeg[grid->prism[grid->nprism].nodes[i]]++;      
    }
  }
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
  int origSize;

  if (grid->npyramid >= grid->maxpyramid) {
    origSize = grid->maxpyramid;
    grid->maxpyramid += 5000;
    if (grid->pyramid == NULL) {
      grid->pyramid = (Pyramid *)malloc(grid->maxpyramid*sizeof(Pyramid));
    }else{
      grid->pyramid =
             (Pyramid *)realloc(grid->pyramid,grid->maxpyramid*sizeof(Pyramid));
      if (NULL != grid->reallocFunc)
	(*grid->reallocFunc)( grid->reallocData, 
			      gridREALLOC_PYRAMID, 
			      origSize, grid->maxpyramid);
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
  int origSize;

  if (grid->nquad >= grid->maxquad) {
    origSize = grid->maxquad;
    grid->maxquad += 5000;
    if (grid->quad == NULL) {
      grid->quad = (Quad *)malloc(grid->maxquad*sizeof(Quad));
    }else{
      grid->quad = (Quad *)realloc(grid->quad,grid->maxquad*sizeof(Quad));
      if (NULL != grid->reallocFunc)
	(*grid->reallocFunc)( grid->reallocData, 
			      gridREALLOC_QUAD, 
			      origSize, grid->maxquad);
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
  int starting_guess, enclosing_cell;
  double xyz[3], bary[4];
  int nodes[4];
  double map0[6], map1[6], map2[6], map3[6];
  int i;
  Grid *child;
  
  if ( !gridValidNode(grid, node) ) return NULL;

  if ( NULL != grid->map ) {
    map[0] = grid->map[0+6*node];
    map[1] = grid->map[1+6*node];
    map[2] = grid->map[2+6*node];
    map[3] = grid->map[3+6*node];
    map[4] = grid->map[4+6*node];
    map[5] = grid->map[5+6*node];
  } else {
    if ( (NULL == grid->child) || ( NULL == grid->child_reference ) ) {
      printf("%s: %d: gridMap: everything is NULL?\n", __FILE__, __LINE__ );
      return NULL;
    }
    child = grid->child;
    starting_guess = grid->child_reference[node];
    gridNodeXYZ( grid, node, xyz );
    enclosing_cell = gridFindEnclosingCell( child, starting_guess, xyz, bary );
    if ( EMPTY == enclosing_cell ) {
      /*
      printf("%s: %d: gridMap: can not find enclosing_cell (EMPTY) for %d\n", 
	     __FILE__, __LINE__, node );
      */
      return NULL;
    }
    grid->child_reference[node] = enclosing_cell;
    gridCell( child, enclosing_cell, nodes );
    gridMap( child, nodes[0], map0 );
    gridMap( child, nodes[1], map1 );
    gridMap( child, nodes[2], map2 );
    gridMap( child, nodes[3], map3 );
    for(i=0;i<6;i++) {
      map[i] = bary[0]*map0[i] + bary[1]*map1[i] 
	     + bary[2]*map2[i] + bary[3]*map3[i];
    }
  }
  
  return grid;
}

Grid *gridSetMap(Grid *grid, int node,
		 double m11, double m12, double m13,
		             double m22, double m23,
		                         double m33)
{
  if ( !gridValidNode(grid, node) ) return NULL;

  /* this function is a no-op if the map array is NULL */
  if ( NULL != grid->map ) {
    grid->map[0+6*node] = m11;
    grid->map[1+6*node] = m12;
    grid->map[2+6*node] = m13;
    grid->map[3+6*node] = m22;
    grid->map[4+6*node] = m23;
    grid->map[5+6*node] = m33;
  }

  return grid;
}

Grid *gridInterpolateMap2(Grid *grid, int node0, int node1, double ratio, 
			  int target)
{
  int i;
  if ( !gridValidNode(grid, node0) ) return NULL;
  if ( !gridValidNode(grid, node1) ) return NULL;
  if ( !gridValidNode(grid, target) ) return NULL;

  if ( NULL != grid->map ) {
    for (i=0;i<6;i++) {
      grid->map[i+6*target] = 
	ratio*(grid->map[i+6*node1]) + (1.0-ratio)*(grid->map[i+6*node0]); 
    }
  }else{
    if ( NULL != grid->child_reference ) {
      double xyz[3], bary[4];
      gridNodeXYZ(grid,target,xyz);
      grid->child_reference[target] = 
	gridFindEnclosingCell(grid->child,grid->child_reference[node0],xyz,bary);
      if (EMPTY == grid->child_reference[target])
	grid->child_reference[target] = 
	  gridFindEnclosingCell(grid->child,grid->child_reference[node1],xyz,bary);
      if (EMPTY == grid->child_reference[target]) {
	printf("%s: %d: gridInterpolateMap2 cannot gridFindEnclosingCell\n",
	       __FILE__,__LINE__);
	return NULL;		      
      }
    }
  }
  return grid;
}

Grid *gridSetCostFunction(Grid *grid, int costFunction)
{
  if ( costFunction < gridCOST_FCN_MEAN_RATIO  || 
       costFunction > gridCOST_FCN_CONFORMITY ) return NULL;
  grid->costFunction = costFunction;
  return grid;
}

Grid *gridSetCostConstraint(Grid *grid, int costConstraint)
{
  if ( costConstraint < 0 || 
       costConstraint > ( gridCOST_CNST_VOLUME | 
			  gridCOST_CNST_VALID  |
			  gridCOST_CNST_AREAUV ) ) return NULL;
  grid->costConstraint = costConstraint;
  return grid;
}

Grid *gridSetMinInsertCost(Grid *grid, double min_cost )
{
  grid->min_allowed_insert_cost = min_cost;
  return grid;
}

double gridMinInsertCost(Grid *grid )
{
  return grid->min_allowed_insert_cost;
}

Grid *gridSetMinSurfaceSmoothCost(Grid *grid, double min_cost )
{
  grid->min_allowed_surface_smooth_cost = min_cost;
  return grid;
}

double gridMinSurfaceSmoothCost(Grid *grid )
{
  return grid->min_allowed_surface_smooth_cost;
}

Grid *gridSetMinSwapCost(Grid *grid, double min_cost )
{
  grid->min_allowed_swap_cost = min_cost;
  return grid;
}

double gridMinSwapCost(Grid *grid )
{
  return grid->min_allowed_swap_cost;
}

Grid *gridSetMinSwapCostImprovement(Grid *grid, double min_cost )
{
  grid->min_allowed_swap_cost_improvement = min_cost;
  return grid;
}

double gridMinSwapCostImprovement(Grid *grid )
{
  return grid->min_allowed_swap_cost_improvement;
}

int gridStoredCostDegree( Grid *grid )
{
  return grid->costDegree;
}

Grid *gridClearStoredCost( Grid *grid )
{
  grid->costDegree = 0;
  return grid;
}

Grid *gridStoreCost( Grid *grid, double cost, double *costDerivative )
{
  if (grid->costDegree == MAXDEG) return NULL;

  grid->storedCost[grid->costDegree] = cost;
  grid->storedCostDerivative[0+3*grid->costDegree] = costDerivative[0];
  grid->storedCostDerivative[1+3*grid->costDegree] = costDerivative[1];
  grid->storedCostDerivative[2+3*grid->costDegree] = costDerivative[2];

  grid->costDegree++;
  return grid;
}

Grid *gridUpdateStoredCost( Grid *grid, int index, 
			    double cost, double *costDerivative )
{
  if ( index < 0 || index >= grid->costDegree) return NULL;

  grid->storedCost[index] = cost;
  grid->storedCostDerivative[0+3*index] = costDerivative[0];
  grid->storedCostDerivative[1+3*index] = costDerivative[1];
  grid->storedCostDerivative[2+3*index] = costDerivative[2];

  return grid;
}

double gridStoredCost( Grid *grid, int index )
{
  if ( index < 0 || index >= gridStoredCostDegree( grid ) ) return DBL_MAX;
  return grid->storedCost[index];
}

Grid *gridStoredCostDerivative( Grid *grid, int index, double *costDerivative )
{
  if ( index < 0 || index >= gridStoredCostDegree( grid ) ) return NULL;
  costDerivative[0] = grid->storedCostDerivative[0+3*index];
  costDerivative[1] = grid->storedCostDerivative[1+3*index];
  costDerivative[2] = grid->storedCostDerivative[2+3*index];
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

int gridMirrorNodeAboutY0(Grid *grid, int node, int origGlobal, int mirrorAux )
{
  int newNode;
  double xyz[3];
  double map[6];
  int aux;
  double Y0;

  Y0 = 0.0;

  gridNodeXYZ(grid,node,xyz);
  if ( origGlobal > 0 ) {
    newNode = gridAddNodeWithGlobal(grid,xyz[0],Y0-(xyz[1]-Y0),xyz[2],
				    gridNodeGlobal(grid,node)+origGlobal);
    gridSetNodePart(grid,newNode,gridNodePart(grid,node));
  } else {
    newNode = gridAddNode(grid,xyz[0],Y0-xyz[1],xyz[2]);
  } 
  gridMap(grid, node, map);
  gridSetMap(grid, newNode, map[0], map[1], map[2], map[3], map[4], map[5]);
  for ( aux = 0; aux<gridNAux(grid);aux++)
    gridSetAux(grid, newNode, aux, gridAux(grid, node, aux));
  if (mirrorAux>EMPTY) 
    gridSetAux(grid, newNode, mirrorAux, -gridAux(grid, node, mirrorAux));
  if (gridNodeFrozen(grid,node)) gridFreezeNode(grid,newNode);
 
  return newNode;
}

Grid *gridCopyAboutY0(Grid *grid, int symmetryFaceId, int mirrorAux )
{
  int node, orignode, origNodeGlobal, origface, origcell, origCellGlobal;
  int *o2n;
  int nodes[4];
  int face, i, faceid;
  int cell;

  printf("gridCopyAboutY0: pack\n");
  if (NULL == gridPack(grid)) {
    printf("gridCopyAboutY0: gridPack failed.\n");
    return NULL;
  }

  orignode = gridNNode(grid);
  origNodeGlobal = gridGlobalNNode(grid);

  o2n = (int *)malloc(sizeof(int) * orignode); 

  printf("gridCopyAboutY0: copy nodes, y = -y\n");

  for ( node = 0 ; node < orignode ; node++){
    o2n[node] = gridMirrorNodeAboutY0(grid,node,origNodeGlobal,mirrorAux);
  }
  if ( origNodeGlobal > 0 ) gridSetGlobalNNode(grid,2*origNodeGlobal);

  printf("gridCopyAboutY0: remove duplicate nodes.\n");

  origface = gridNFace(grid);
  for ( face = 0 ; face < origface ; face++ ){
    gridFace(grid,face,nodes,&faceid);
    if ( faceid == symmetryFaceId ) {
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
  origCellGlobal = gridGlobalNCell(grid);
  for ( cell = 0 ; cell < origcell ; cell++ ){
    gridCell(grid,cell,nodes);
    if ( origCellGlobal > 0 ) {
      gridAddCellWithGlobal(grid,
			    o2n[nodes[1]],
			    o2n[nodes[0]],
			    o2n[nodes[2]],
			    o2n[nodes[3]],
			    origCellGlobal+gridCellGlobal(grid,cell));
    } else {
      gridAddCell(grid,
		  o2n[nodes[1]],
		  o2n[nodes[0]],
		  o2n[nodes[2]],
		  o2n[nodes[3]]);
    }
  }
  if ( origCellGlobal > 0 ) gridSetGlobalNCell(grid,2*origCellGlobal);
  printf("gridCopyAboutY0: remove sym face\n");

  gridThawAll(grid);
  gridDeleteThawedFaces(grid,symmetryFaceId);

  free(o2n);
    
  return grid;
}

Grid *gridReportZeroDegreeNodes(Grid *grid)
{
  int node;
  double xyz[3];
  for( node = 0 ; node < gridMaxNode(grid) ; node++ ) {
    if ( gridValidNode(grid, node) && 0 == gridCellDegree(grid, node) ) {
      printf("node %d has zero cell degree\n",node);
      if(grid==gridNodeXYZ(grid,node,xyz)) printf("at %f %f %f\n",
						  xyz[0],xyz[1],xyz[2]);
    }
  }
  return grid;
}

int gridPhase(Grid *grid){
  return grid->phase;
}

Grid *gridSetPhase(Grid *grid, int phase){
  grid->phase = phase;
  return grid;
}

Grid *gridCacheCurrentGridAndMap(Grid *grid){
  int node, cell;
  grid->child = gridDup( grid );
  free(grid->map); grid->map = NULL;
  grid->child_reference = (int *)malloc(gridMaxNode(grid) * sizeof(int));
  for (node=0;node<=gridMaxNode(grid); node++){
    if ( gridValidNode(grid,node) ) {
      cell = adjItem(adjFirst(gridCellAdj(grid),node));
      grid->child_reference[node] = cell;
    }
  }
  return grid;
}
