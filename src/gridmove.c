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
#include "gridmath.h"
#include "gridcad.h"
#include "gridmove.h"

GridMove *gridmoveCreate( Grid *grid )
{
  int i;
  GridMove *gm;
  gm = malloc(sizeof(GridMove));
  gm->grid = grid;

  gridAttachPacker( grid, gridmovePack, (void *)gm );
  gridAttachNodeSorter( grid, gridmoveSortNode, (void *)gm );
  gridAttachReallocator( grid, gridmoveReallocator, (void *)gm );
  gridAttachFreeNotifier( grid, gridmoveGridHasBeenFreed, (void *)gm );

  gm->displacement = malloc(3*gridMaxNode(grid)*sizeof(double));
  for (i=0;i<3*gridMaxNode(grid);i++) gm->displacement[i] = 0.0;

  gm->specified = malloc(gridMaxNode(grid)*sizeof(GridBool));
  for (i=0;i<gridMaxNode(grid);i++) gm->specified[i] = FALSE;

  gm->relaxationScheme = gridmoveEMPTY_SCHEME;

  gm->c2e = NULL;
  gm->springs = NULL;
  gm->xyz = NULL;
  gm->k = NULL;
  gm->source = NULL;

  gm->ksum = NULL;
  gm->kxyz = NULL;

  gm->rowStart = NULL;
  gm->compRow = NULL;
  gm->a = NULL;
  gm->lu = NULL;
  gm->lastxyz = NULL;

  return gm;
}

Grid *gridmoveGrid(GridMove *gm)
{
  return gm->grid;
}

GridMove *gridmoveFreeRelaxation(GridMove *gm)
{
  GridBool w;
  w = FALSE;

  if (NULL != gm->c2e)      { free(gm->c2e);      gm->c2e      = NULL; w=TRUE; }
  if (NULL != gm->springs)  { free(gm->springs);  gm->springs  = NULL; w=TRUE; }
  if (NULL != gm->xyz)      { free(gm->xyz);      gm->xyz      = NULL; w=TRUE; }
  if (NULL != gm->k)        { free(gm->k);        gm->k        = NULL; w=TRUE; }
  if (NULL != gm->source)   { free(gm->source);   gm->source   = NULL; w=TRUE; }

  if (NULL != gm->ksum)     { free(gm->ksum);     gm->ksum     = NULL; w=TRUE; }
  if (NULL != gm->kxyz)     { free(gm->kxyz);     gm->kxyz     = NULL; w=TRUE; }

  if (NULL != gm->rowStart) { free(gm->rowStart); gm->rowStart = NULL; w=TRUE; }
  if (NULL != gm->compRow)  { free(gm->compRow);  gm->compRow  = NULL; w=TRUE; }
  if (NULL != gm->a)        { free(gm->a);        gm->a        = NULL; w=TRUE; }
  if (NULL != gm->lu)       { free(gm->lu);       gm->lu       = NULL; w=TRUE; }
  if (NULL != gm->lastxyz)  { free(gm->lastxyz);  gm->lastxyz  = NULL; w=TRUE; }

  if (w) {
    printf("ERROR: A grid resize or renumber has called in the middle of a \n");
    printf("       grimove relaxation. You are most likely totally hosed.\n");
    printf("%s: %d: More code needs to be written.\n",__FILE__,__LINE__);
    return NULL;
  }

  return gm;
}

void gridmoveFree(GridMove *gm)
{
  if (NULL != gm->lastxyz) free(gm->lastxyz);
  if (NULL != gm->lu) free(gm->lu);
  if (NULL != gm->a) free(gm->a);
  if (NULL != gm->compRow) free(gm->compRow);
  if (NULL != gm->rowStart) free(gm->rowStart);

  if (NULL != gm->ksum) free(gm->ksum);
  if (NULL != gm->kxyz) free(gm->kxyz);

  if (NULL != gm->c2e) free(gm->c2e);
  if (NULL != gm->springs) free(gm->springs);
  if (NULL != gm->xyz) free(gm->xyz);
  if (NULL != gm->k) free(gm->k);
  if (NULL != gm->source) free(gm->source);

  free(gm->specified);
  free(gm->displacement);
  if (NULL != gm->grid) { 
    gridDetachPacker( gm->grid );
    gridDetachNodeSorter( gm->grid );
    gridDetachReallocator( gm->grid );
    gridDetachFreeNotifier( gm->grid );
  }
  free(gm);
}

void gridmovePack(void *voidGridMove, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n)
{
  GridMove *gm = (GridMove *)voidGridMove;
  int orignode, packnode;
  int ixyz;
  for ( orignode = 0 ; orignode < maxnode ; orignode++ ){
    packnode = nodeo2n[orignode];
    if (EMPTY!=packnode) {
      for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
	gm->displacement[ixyz+3*packnode] = gm->displacement[ixyz+3*orignode];
      }
    }
  }
  for ( packnode=nnode ; packnode < maxnode ; packnode++ ){ 
    for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
      gm->displacement[ixyz+3*packnode] = 0.0;
    }
  }

  for ( orignode = 0 ; orignode < maxnode ; orignode++ ){
    packnode = nodeo2n[orignode];
    if (EMPTY!=packnode) {
      gm->specified[packnode] = gm->specified[orignode];
    }
  }
  for ( packnode=nnode ; packnode < maxnode ; packnode++ ){ 
    gm->specified[packnode] = FALSE;
  }

  gridmoveFreeRelaxation(gm);
}

void gridmoveSortNode(void *voidGridMove, int maxnode, int *o2n)
{
  GridMove *gm = (GridMove *)voidGridMove;
  int node, ixyz;
  double *temp_xyz;
  GridBool *temp_bool;

  temp_xyz = malloc( maxnode * sizeof(double) );
  for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
    for ( node = 0 ; node < maxnode ; node++ )temp_xyz[node]=0.0;
    for ( node = 0 ; node < maxnode ; node++ ){
      if (EMPTY != o2n[node])
	temp_xyz[o2n[node]] = gm->displacement[ixyz+3*node];
    }
    for ( node = 0 ; node < maxnode ; node++ ){
      gm->displacement[ixyz+3*node] = temp_xyz[node];
    }
  }
  free(temp_xyz);

  temp_bool = malloc( maxnode * sizeof(GridBool) );
  for ( node = 0 ; node < maxnode ; node++ )temp_bool[node]=FALSE;
  for ( node = 0 ; node < maxnode ; node++ ){
    if (EMPTY != o2n[node])
      temp_bool[o2n[node]] = gm->specified[node];
  }
  for ( node = 0 ; node < maxnode ; node++ ){
    gm->specified[node] = temp_bool[node];
  }
  free(temp_bool);

  gridmoveFreeRelaxation(gm);
}

void gridmoveReallocator(void *voidGridMove, int reallocType, 
			 int lastSize, int newSize)
{
  GridMove *gm = (GridMove *)voidGridMove;
  int i;
  if (gridREALLOC_NODE == reallocType) {
    gm->displacement = realloc(gm->displacement, 3*newSize*sizeof(double));
    for (i=3*lastSize;i<3*newSize;i++) gm->displacement[i] = 0.0;
    gm->specified = realloc(gm->specified, newSize*sizeof(GridBool));
    for (i=lastSize;i<newSize;i++) gm->specified[i] = FALSE;
  }

  gridmoveFreeRelaxation(gm);
}

void gridmoveGridHasBeenFreed(void *voidGridMove )
{
  GridMove *gm = (GridMove *)voidGridMove;
  gm->grid = NULL;
}

GridMove *gridmoveDisplace(GridMove *gm, int node, double *displace)
{
  if (node < 0 || node >= gridMaxNode(gm->grid)) return NULL;
  gm->displacement[0+3*node] = displace[0];
  gm->displacement[1+3*node] = displace[1];
  gm->displacement[2+3*node] = displace[2];
  gm->specified[node] = TRUE;
  return gm;
}

GridMove *gridmoveDisplacement(GridMove *gm, int node, double *displacement)
{
  if (node < 0 || node >= gridMaxNode(gm->grid)) return NULL;
  displacement[0] = gm->displacement[0+3*node];
  displacement[1] = gm->displacement[1+3*node];
  displacement[2] = gm->displacement[2+3*node];
  return gm;
}

GridBool gridmoveSpecified(GridMove *gm, int node)
{
  if (node < 0 || node >= gridMaxNode(gm->grid)) return FALSE;
  return gm->specified[node];
}

GridMove *gridmoveCellFaceNormals(GridMove *gm, double *xyz, int *nodes, 
				  double normals[4][3])
{
  int edge, n0, n1, face;
  double e[6][3];
  int edge2node0[6] = {0, 0, 0, 1, 1, 2};
  int edge2node1[6] = {1, 2, 3, 2, 3, 3};
  int face2edge0[6] = {0, 2, 4, 1};
  int face2edge1[6] = {1, 0, 3, 2};

  for(edge=0;edge<6;edge++) {
    n0 = nodes[edge2node0[edge]];
    n1 = nodes[edge2node1[edge]];
    gridSubtractVector(&xyz[3*n1],&xyz[3*n0],e[edge]);
    gridVectorNormalize(e[edge]);
  }
  for(face=0;face<4;face++){
    gridCrossProduct(e[face2edge0[face]],e[face2edge1[face]],normals[face]);
    gridVectorNormalize(normals[face]);
  }
  return gm;
}

GridMove *gridmoveSpringConstant(GridMove *gm, double *xyz, int nsprings, 
				 double *k, int *springs, int *c2e)
{
  Grid *grid = gridmoveGrid(gm);
  int s, n0, n1;
  double dxyz[3];
  int cell, edge;
  int nodes[4];
  double n[4][3];
  int edge2face0[6] = {2, 1, 0, 1, 0, 0}; /* end points */
  int edge2face1[6] = {3, 2, 2, 3, 3, 1};
  double angle, invsinangle;

  /* touching
  int edge2face0[6] = {0, 0, 1, 0, 1, 2};
  int edge2face1[6] = {1, 3, 3, 2, 2, 3};
  */

  for(s=0;s<nsprings;s++) {
    n0 = springs[0+2*s];
    n1 = springs[1+2*s];
    gridSubtractVector(&xyz[3*n1],&xyz[3*n0],dxyz);
    k[s]= 1.0/gridVectorLength(dxyz);
  }

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid == gridCell(grid,cell,nodes)) {
      gridmoveCellFaceNormals( gm, xyz, nodes, n);
      for(edge=0;edge<6;edge++) {
	angle = acos(-gridDotProduct(n[edge2face0[edge]],n[edge2face1[edge]]));
	invsinangle = 1.0 / sin(angle);
	k[c2e[edge+6*cell]] += invsinangle*invsinangle;
      }
    }
  }

  return gm;
}

GridMove *gridmoveSpringRelaxationStartUp(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node;
  gm->c2e = malloc(6*gridMaxCell(grid)*sizeof(int));
  if (gm != gridmoveComputeC2E(gm, &(gm->nsprings), gm->c2e) ) { 
    free(gm->c2e); gm->c2e=NULL;
    return NULL; 
  }
  gm->springs = malloc(2*gm->nsprings*sizeof(int));
  if (gm != gridmoveComputeSpringsWithC2E(gm, gm->c2e, gm->springs)) {
    free(gm->c2e); gm->c2e=NULL;
    free(gm->springs); gm->springs=NULL;
    return NULL;
  }
  
  gm->xyz = malloc(3*gridMaxNode(grid)*sizeof(double));
  gm->k = malloc(gm->nsprings*sizeof(double));
  gm->source = malloc(3*gridMaxNode(grid)*sizeof(double));

  gm->ksum = malloc(gridMaxNode(grid)*sizeof(double));
  gm->kxyz = malloc(3*gridMaxNode(grid)*sizeof(double));
 
  for(node=0;node<gridMaxNode(grid);node++)
    gridNodeXYZ(grid,node,&(gm->xyz[3*node]));

  return gm;
}

GridMove *gridmoveSpringRelaxationStartStep(GridMove *gm, double position)
{
  Grid *grid = gridmoveGrid(gm);
  int node, s, i;
  int n0, n1;
  double res[3];

  gridmoveSpringConstant(gm, gm->xyz, 
			 gm->nsprings, gm->k, gm->springs,
			 gm->c2e);
  
  for(node=0;node<3*gridMaxNode(grid);node++) gm->source[node]=0.0;
  for(s=0;s<gm->nsprings;s++) {
    n0 = gm->springs[0+2*s];
    n1 = gm->springs[1+2*s];
    for(i=0;i<3;i++) {
      res[i] = gm->k[s] * ( gm->xyz[i+3*n0] - gm->xyz[i+3*n1] );
      gm->source[i+3*n0] += res[i];
      gm->source[i+3*n1] -= res[i];
    }
  }

  for(node=0;node<gridMaxNode(grid);node++) {
    if ( gridmoveSpecified(gm,node) && 
	 gridNodeLocal(grid,node) &&
	 grid == gridNodeXYZ(grid,node,&(gm->xyz[3*node]))) {
      for(i=0;i<3;i++) {
	gm->xyz[i+3*node] += position*gm->displacement[i+3*node];
      }
    }
  }

  return gm;
}

GridMove *gridmoveSpringRelaxationSubIteration(GridMove *gm, double *residual2)
{
  Grid *grid = gridmoveGrid(gm);
  int node, s, i;
  int n0, n1;
  double res[3];
  double residual;

  for(node=0;node<gridMaxNode(grid);node++) gm->ksum[node]=0.0;
  for(node=0;node<3*gridMaxNode(grid);node++) gm->kxyz[node]=0.0;
  for(s=0;s<gm->nsprings;s++) {
    n0 = gm->springs[0+2*s];
    n1 = gm->springs[1+2*s];
    gm->ksum[n0] += gm->k[s];
    gm->ksum[n1] += gm->k[s];
    for(i=0;i<3;i++) {
      gm->kxyz[i+3*n0] += gm->k[s] * gm->xyz[i+3*n1];
      gm->kxyz[i+3*n1] += gm->k[s] * gm->xyz[i+3*n0];
    }
  }

  residual = 0.0;
  for(node=0;node<gridMaxNode(grid);node++)
    if ( gridValidNode(grid,node) && 
	 gridNodeLocal(grid,node) &&
	 !gridmoveSpecified(gm,node)) {
      for(i=0;i<3;i++) {
	res[i] = gm->xyz[i+3*node] - 
	  ( gm->kxyz[i+3*node] + gm->source[i+3*node] ) / gm->ksum[node];
	gm->xyz[i+3*node] = (  gm->kxyz[i+3*node] + gm->source[i+3*node] ) 
	                    / gm->ksum[node];
      }
      residual += gridDotProduct(res,res);
    }

  *residual2 = residual;

  return gm;
}

GridMove *gridmoveSpringRelaxationShutDown(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node, i;
  double xyz0[3];

  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridValidNode(grid,node) && !gridmoveSpecified(gm,node)) {
      gridNodeXYZ(grid,node,xyz0);
      for(i=0;i<3;i++) gm->displacement[i+3*node] = gm->xyz[i+3*node] - xyz0[i];
    }
  }

  free(gm->ksum);    gm->ksum=NULL;
  free(gm->kxyz);    gm->kxyz=NULL;

  free(gm->c2e);     gm->c2e=NULL;
  free(gm->springs); gm->springs=NULL;
  free(gm->xyz);     gm->xyz=NULL;
  free(gm->k);       gm->k=NULL;
  free(gm->source);  gm->source=NULL;

  return gm;
}

GridMove *gridmoveSpringRelaxation(GridMove *gm, int nsteps, int subIterations)
{
  int step, iteration;
  double position;
  double rmsResidual;

  if (gm != gridmoveSpringRelaxationStartUp(gm)) return NULL;
  for(step=0;step<nsteps;step++) {
    position = (double)(step+1)/(double)nsteps;
    gridmoveSpringRelaxationStartStep(gm, position);    
    for(iteration=0;iteration<subIterations;iteration++) {
      gridmoveSpringRelaxationSubIteration(gm, &rmsResidual);
      /* printf("Iteration %4d Residual %23.15e\n",iteration,rmsResidual); */ 
    }
  }
  gridmoveSpringRelaxationShutDown(gm);

  return gm;
}

void grimoveAddEdgeToC2E(Grid *grid, int *c2e, int edge, int node0, int node1)
{
  int iedge;
  int cell, nodes[4];
  int edge2node0[6] = {0, 0, 0, 1, 1, 2};
  int edge2node1[6] = {1, 2, 3, 2, 3, 3};
  int local0, local1;
  AdjIterator it;

  for ( it = adjFirst(gridCellAdj(grid),node0); 
	adjValid(it); 
	it = adjNext(it) ) {
    cell = adjItem(it);
    gridCell( grid, cell, nodes );
    for(iedge=0;iedge<6;iedge++) {
      local0 = nodes[edge2node0[iedge]];
      local1 = nodes[edge2node1[iedge]];
      if ( MIN(local0, local1) == MIN(node0,node1) &&
	   MAX(local0, local1) == MAX(node0,node1) ) {
	c2e[iedge+6*cell] = edge;
      } 
    }
  }

}

GridMove *gridmoveComputeC2E(GridMove *gm, int *nedge, int *c2e)
{
  Grid *grid = gridmoveGrid(gm);
  int cell, edge, count;
  int nodes[4];
  int edge2node0[6] = {0, 0, 0, 1, 1, 2};
  int edge2node1[6] = {1, 2, 3, 2, 3, 3};

  for(cell=0;cell<6*gridMaxCell(grid);cell++) c2e[cell] = EMPTY;

  count = 0;
  for(cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid == gridCell(grid,cell,nodes)) {
      for(edge=0;edge<6;edge++) {
	if ( EMPTY == c2e[edge+6*cell] ) {
	  c2e[edge+6*cell] = count;
	  grimoveAddEdgeToC2E(grid, c2e, count,
			      nodes[edge2node0[edge]], nodes[edge2node1[edge]]);
	  count++;
	}
      }
    }
  }
  *nedge = count;
  return gm;
}

GridMove *gridmoveComputeSpringsWithC2E(GridMove *gm, int *c2e, int *springs)
{
  Grid *grid = gridmoveGrid(gm);
  int cell, edge;
  int nodes[4];
  int edge2node0[6] = {0, 0, 0, 1, 1, 2};
  int edge2node1[6] = {1, 2, 3, 2, 3, 3};
  int node0, node1;

  for(cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid == gridCell(grid,cell,nodes)) {
      for(edge=0;edge<6;edge++) {
	if (EMPTY!=c2e[edge+6*cell]) {
	  node0 = MIN(nodes[edge2node0[edge]],nodes[edge2node1[edge]]);
	  node1 = MAX(nodes[edge2node0[edge]],nodes[edge2node1[edge]]);
	  springs[0+2*c2e[edge+6*cell]] = node0;
	  springs[1+2*c2e[edge+6*cell]] = node1;
	}else{
	  printf("ERROR: %s: %d: c2e EMPTY.\n",__FILE__,__LINE__);
	}
      }
    }
  }

  return gm;
}

GridMove *gridmoveSprings(GridMove *gm, int *nsprings, int **springs)
{
  int *c2e;

  c2e = malloc(6*gridMaxCell(gridmoveGrid(gm))*sizeof(int));
  gridmoveComputeC2E(gm, nsprings, c2e);
  *springs = malloc(2*(*nsprings)*sizeof(int));
  gridmoveComputeSpringsWithC2E(gm, c2e, *springs);

  free(c2e);
  return gm;
}

GridMove *gridmoveApplyDisplacements(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node;
  double xyz[3], displacement[3];

  for ( node=0 ; node<gridMaxNode(grid) ; node++ ) {
    if ( gridNodeLocal(grid,node) && 
	 grid == gridNodeXYZ(grid,node,xyz)) {
      gridmoveDisplacement(gm,node,displacement);
      gridAddToVector(xyz,displacement);
      gridSetNodeXYZ(grid,node,xyz);
    }
  }

  return gm;
}

GridMove *gridmoveProjectionDisplacements(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node;
  double displacement[3];

  for ( node=0 ; node<gridMaxNode(grid) ; node++ ) {
    if ( gridNodeLocal(grid,node) && gridGeometryFace( grid, node ) ) {
      gridNodeProjectionDisplacement(grid,node,displacement);
      gridmoveDisplace(gm,node,displacement);
    }
  }

  return gm;
}

GridMove *gridmoveDataLeadingDimension( GridMove *gm, int *ndim )
{
  *ndim = 3;
  return gm;
}

GridMove *gridmoveInitializeMPITest(GridMove *gm )
{
  int node,i;
  Grid *grid = gridmoveGrid(gm);
  
  gm->xyz = malloc(3*gridMaxNode(grid)*sizeof(double));

  for(node=0;node<gridMaxNode(grid);node++) {
    for(i=0;i<3;i++) {
      gm->xyz[i+3*node] = -9.0 - 0.1*(double)i;
      if ( gridValidNode(grid,node) ) {
	gm->xyz[i+3*node] = -1.0 - 0.1*(double)i;
	if ( gridNodeLocal(grid,node) ) {
	  gm->xyz[i+3*node] = (double)gridNodeGlobal(grid,node) + 0.1*(double)i;
	}
      }
    }
  }

  return gm;
}

GridMove *gridmoveCompleteMPITest(GridMove *gm )
{
  int node,i;
  double correct, error;
  Grid *grid = gridmoveGrid(gm);

  for(node=0;node<gridMaxNode(grid);node++) {
    for(i=0;i<3;i++) {
      if ( gridValidNode(grid,node) ) {
	correct = (double)gridNodeGlobal(grid,node) + 0.1*(double)i;
	error = ABS(gm->xyz[i+3*node] - correct);
	if ( error > 1e-12) {
	  printf("%s: %d: ERROR mpi part %d node %d i %d error %e\n",
		 __FILE__,__LINE__,gridPartId,node,i,error);

	}
      }
    }
  }
  free(gm->xyz); gm->xyz = NULL;
  return gm;
}

GridMove *gridmoveLoadFortranNodeData( GridMove *gm, int nnode, 
				     int *nodes, double *data)
{
  int node, localnode, i;
  localnode = 0;
  node = 0;
  for (node=0;node<nnode;node++) {
    if (nodes[node] > 0) {
      localnode = nodes[node]-1;
      for(i=0;i<3;i++) data[i+3*node] = gm->xyz[i+3*localnode];
    }
  }
  return gm;
}

GridMove *gridmoveSetFortranNodeData( GridMove *gm, int nnode, 
				    int *nodes, double *data)
{
  int node, localnode, i;
  localnode = 0;
  node = 0;
  for (node=0;node<nnode;node++) {
    if (nodes[node] > 0) {
      localnode = nodes[node]-1;
      for(i=0;i<3;i++) gm->xyz[i+3*localnode] = data[i+3*node] ;
    }
  }
  return gm;
}

GridMove *gridmoveInitializeCompRow(GridMove *gm)
{
  int node, nnode;
  int node0, node1;
  int spring, nsprings, *springs;
  int row, fix, lowEntry, lowNode, findLower;

  if (NULL != gm->rowStart) { free(gm->rowStart); gm->rowStart = NULL; }
  if (NULL != gm->compRow) { free(gm->compRow); gm->compRow = NULL; }

  gridmoveSprings(gm, &nsprings, &springs);
  nnode = gridMaxNode(gridmoveGrid(gm));
  gm->rowStart = malloc( (1+nnode) *sizeof(int));
  for ( node = 0 ; node < nnode+1 ; node++ ) gm->rowStart[node] = 0;
  for ( spring = 0 ; spring < nsprings ; spring ++ ) {
    gm->rowStart[springs[0+2*spring]]++;
    gm->rowStart[springs[1+2*spring]]++;
  }
  for ( node = 0 ; node < nnode ; node++ ) 
    if (gm->rowStart[node]>0) gm->rowStart[node]++; /*diagonal*/
  for ( node = 0 ; node < nnode ; node++ ) 
    gm->rowStart[node+1] += gm->rowStart[node];
  gm->compRow = malloc( gm->rowStart[nnode] *sizeof(int));

  for ( node = nnode ; node > 0 ; node-- ) {
    if (gm->rowStart[node] > gm->rowStart[node-1]) {
      gm->rowStart[node]--;
      gm->compRow[gm->rowStart[node]] = node;
    }
  }
  if (gm->rowStart[0] > 0){
    gm->rowStart[node]--;
    gm->compRow[gm->rowStart[node]] = node;
  }

  for ( spring = 0 ; spring < nsprings ; spring ++ ) {
    node0 = springs[0+2*spring];
    node1 = springs[1+2*spring];
    gm->rowStart[node0]--;
    gm->compRow[gm->rowStart[node0]] = node1;
    node0 = springs[1+2*spring];
    node1 = springs[0+2*spring];
    gm->rowStart[node0]--;
    gm->compRow[gm->rowStart[node0]] = node1;
  }

  for ( row = 0 ; row < nnode ; row++ ) {
  
    for ( fix = gridmoveRowStart(gm, row) ;
	  fix < gridmoveRowStart(gm, row+1) ;
	  fix++ ) {
      lowEntry = fix;
      lowNode = gridmoveRowNode(gm, lowEntry);
      for ( findLower = fix+1 ;
	    findLower < gridmoveRowStart(gm, row+1) ;
	    findLower++ ) {
	node = gridmoveRowNode(gm, findLower);
	if ( node < lowNode) {
	  lowEntry = findLower;
	  lowNode = node;
	}
      }
      if (fix != lowEntry) {
	gm->compRow[lowEntry] = gm->compRow[fix];
	gm->compRow[fix] = lowNode;
      }
    }
  }

  free(springs);

  return gm;
}

int gridmoveRowStart(GridMove *gm, int row)
{
  if (row<0 || row>gridMaxNode(gridmoveGrid(gm)) ) return EMPTY;
  if (NULL == gm->rowStart) gridmoveInitializeCompRow(gm);
  return gm->rowStart[row];
}

int gridmoveNNZ(GridMove * gm)
{
  if (NULL == gm->rowStart) gridmoveInitializeCompRow(gm);
  return gm->rowStart[gridMaxNode(gridmoveGrid(gm))];
}

int gridmoveRowNode(GridMove *gm, int entry)
{
  if (entry<0 || entry>=gridmoveNNZ(gm) ) return EMPTY;
  if (NULL == gm->rowStart) gridmoveInitializeCompRow(gm);
  return gm->compRow[entry];
}

int gridmoveRowEntry(GridMove *gm, int row, int node)
{
  int entry;
  if (row<0 || row>=gridMaxNode(gridmoveGrid(gm)) ) return EMPTY;
  if (NULL == gm->rowStart) gridmoveInitializeCompRow(gm);
  for ( entry = gridmoveRowStart(gm, row) ;
	entry < gridmoveRowStart(gm, row+1) ;
	entry++ ) {
    if (gm->compRow[entry]==node) return entry;
    if (gm->compRow[entry]>node) return EMPTY;
  }
  return EMPTY;
}

GridMove *gridmoveElasticRelaxationStartUp(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node;
  
  gridmoveInitializeCompRow(gm);

  gm->a  = malloc(9*gridmoveNNZ(gm)*sizeof(double));

  gm->lu = malloc(9*gridMaxNode(grid)*sizeof(double));

  gm->lastxyz = malloc(3*gridMaxNode(grid)*sizeof(double));
  gm->xyz     = malloc(3*gridMaxNode(grid)*sizeof(double));
 
  for(node=0;node<3*gridMaxNode(grid);node++)
    gm->lastxyz[node] = 0.0;

  for(node=0;node<gridMaxNode(grid);node++)
    gridNodeXYZ(grid,node,&(gm->lastxyz[3*node]));

  for(node=0;node<3*gridMaxNode(grid);node++)
    gm->xyz[node] = 0.0;

  return gm;
}

GridMove *gridmoveElasticRelaxationStartStep(GridMove *gm, double position)
{
  Grid *grid = gridmoveGrid(gm);
  int i;
  int cell, nodes[4];
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double nx1, ny1, nz1;
  double nx2, ny2, nz2;
  double nx3, ny3, nz3;
  double nx4, ny4, nz4;
  double vol;

  int idiag1, idiag2, idiag3, idiag4;
  int ioff12, ioff13, ioff14;
  int ioff21, ioff23, ioff24;
  int ioff31, ioff32, ioff34;
  int ioff41, ioff42, ioff43;

  double c;
  double c1x, c1y, c1z;
  double c2x, c2y, c2z;
  double c3x, c3y, c3z;
  double c4x, c4y, c4z;

  double xn, yn, zn;

  double R, onePlusR;

  int node, entry;
  double xyz[3];


  for(node=0;node<3*gridMaxNode(grid);node++)
    gm->lastxyz[node] += gm->xyz[node];

  for(node=0;node<3*gridMaxNode(grid);node++)
    gm->xyz[node] = 0.0;

  for(node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode(grid,node) &&
	 gridmoveSpecified(gm,node) && 
	 gridNodeLocal(grid,node) ) {
      gridNodeXYZ(grid,node,xyz);
      for(i=0;i<3;i++)
	gm->xyz[i+3*node] = position*gm->displacement[i+3*node]+xyz[i] 
	  - gm->lastxyz[i+3*node];
    }
  }

  for(i=0;i<9*gridmoveNNZ(gm);i++) gm->a[i]=0.0;

  for(cell=0;cell<gridMaxCell(grid);cell++){
    if (grid==gridCell(grid,cell,nodes)) {

      x1 = gm->lastxyz[0+3*nodes[0]];
      y1 = gm->lastxyz[1+3*nodes[0]];
      z1 = gm->lastxyz[2+3*nodes[0]];
      x2 = gm->lastxyz[0+3*nodes[1]];
      y2 = gm->lastxyz[1+3*nodes[1]];
      z2 = gm->lastxyz[2+3*nodes[1]];
      x3 = gm->lastxyz[0+3*nodes[2]];
      y3 = gm->lastxyz[1+3*nodes[2]];
      z3 = gm->lastxyz[2+3*nodes[2]];
      x4 = gm->lastxyz[0+3*nodes[3]];
      y4 = gm->lastxyz[1+3*nodes[3]];
      z4 = gm->lastxyz[2+3*nodes[3]];

      nx1 = 0.5*((y2 - y4)*(z3 - z4) - (y3 - y4)*(z2 - z4));
      ny1 = 0.5*((z2 - z4)*(x3 - x4) - (z3 - z4)*(x2 - x4));
      nz1 = 0.5*((x2 - x4)*(y3 - y4) - (x3 - x4)*(y2 - y4));

      nx2 = 0.5*((y3 - y4)*(z1 - z4) - (y1 - y4)*(z3 - z4));
      ny2 = 0.5*((z3 - z4)*(x1 - x4) - (z1 - z4)*(x3 - x4));
      nz2 = 0.5*((x3 - x4)*(y1 - y4) - (x1 - x4)*(y3 - y4));

      nx3 = 0.5*((y1 - y4)*(z2 - z4) - (y2 - y4)*(z1 - z4));
      ny3 = 0.5*((z1 - z4)*(x2 - x4) - (z2 - z4)*(x1 - x4));
      nz3 = 0.5*((x1 - x4)*(y2 - y4) - (x2 - x4)*(y1 - y4));

      nx4 = -nx1 -nx2 -nx3;
      ny4 = -ny1 -ny2 -ny3;
      nz4 = -nz1 -nz2 -nz3;

      vol = ( ((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))*(x4-x1) -
	      ((x2-x1)*(z3-z1) - (x3-x1)*(z2-z1))*(y4-y1) +
	      ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))*(z4-z1) ) / 6.0;

      if(vol <= 0.0)printf("%s: %d: Negative Vol %f %f %f\n",
			   __FILE__,__LINE__,x1,y1,z1);

      idiag1 = 9*gridmoveRowEntry(gm,nodes[0],nodes[0]);
      ioff12 = 9*gridmoveRowEntry(gm,nodes[0],nodes[1]);
      ioff13 = 9*gridmoveRowEntry(gm,nodes[0],nodes[2]);
      ioff14 = 9*gridmoveRowEntry(gm,nodes[0],nodes[3]);

      ioff21 = 9*gridmoveRowEntry(gm,nodes[1],nodes[0]);
      idiag2 = 9*gridmoveRowEntry(gm,nodes[1],nodes[1]);
      ioff23 = 9*gridmoveRowEntry(gm,nodes[1],nodes[2]);
      ioff24 = 9*gridmoveRowEntry(gm,nodes[1],nodes[3]);

      ioff31 = 9*gridmoveRowEntry(gm,nodes[2],nodes[0]);
      ioff32 = 9*gridmoveRowEntry(gm,nodes[2],nodes[1]);
      idiag3 = 9*gridmoveRowEntry(gm,nodes[2],nodes[2]);
      ioff34 = 9*gridmoveRowEntry(gm,nodes[2],nodes[3]);

      ioff41 = 9*gridmoveRowEntry(gm,nodes[3],nodes[0]);
      ioff42 = 9*gridmoveRowEntry(gm,nodes[3],nodes[1]);
      ioff43 = 9*gridmoveRowEntry(gm,nodes[3],nodes[2]);
      idiag4 = 9*gridmoveRowEntry(gm,nodes[3],nodes[3]);

      c = 1.0/(3.0*vol);
      c1x = c*(nx2 + nx3 + nx4);
      c2x = c*(nx3 + nx4 + nx1);
      c3x = c*(nx4 + nx1 + nx2);
      c4x = c*(nx1 + nx2 + nx3);

      c1y = c*(ny2 + ny3 + ny4);
      c2y = c*(ny3 + ny4 + ny1);
      c3y = c*(ny4 + ny1 + ny2);
      c4y = c*(ny1 + ny2 + ny3);

      c1z = c*(nz2 + nz3 + nz4);
      c2z = c*(nz3 + nz4 + nz1);
      c3z = c*(nz4 + nz1 + nz2);
      c4z = c*(nz1 + nz2 + nz3);

      R = 1.0 / vol;
      onePlusR = 1.0 + R;

      /* node1 */
      xn = nx1/3.0;
      yn = ny1/3.0;
      zn = nz1/3.0;

      gm->a[0+idiag1] = gm->a[0+idiag1] + xn*onePlusR*c1x + yn*c1y+ zn*c1z;
      gm->a[1+idiag1] = gm->a[1+idiag1] + xn*R*c1y;
      gm->a[2+idiag1] = gm->a[2+idiag1] + xn*R*c1z;
      gm->a[3+idiag1] = gm->a[3+idiag1] + yn*R*c1x;
      gm->a[4+idiag1] = gm->a[4+idiag1] + yn*onePlusR*c1y + zn*c1z+ xn*c1x;
      gm->a[5+idiag1] = gm->a[5+idiag1] + yn*R*c1z;
      gm->a[6+idiag1] = gm->a[6+idiag1] + zn*R*c1x;
      gm->a[7+idiag1] = gm->a[7+idiag1] + zn*R*c1y;
      gm->a[8+idiag1] = gm->a[8+idiag1] + zn*onePlusR*c1z + xn*c1x+ yn*c1y;

      gm->a[0+ioff12] = gm->a[0+ioff12] + xn*onePlusR*c2x + yn*c2y+ zn*c2z;
      gm->a[1+ioff12] = gm->a[1+ioff12] + xn*R*c2y;
      gm->a[2+ioff12] = gm->a[2+ioff12] + xn*R*c2z;
      gm->a[3+ioff12] = gm->a[3+ioff12] + yn*R*c2x;
      gm->a[4+ioff12] = gm->a[4+ioff12] + yn*onePlusR*c2y + zn*c2z+ xn*c2x;
      gm->a[5+ioff12] = gm->a[5+ioff12] + yn*R*c2z;
      gm->a[6+ioff12] = gm->a[6+ioff12] + zn*R*c2x;
      gm->a[7+ioff12] = gm->a[7+ioff12] + zn*R*c2y;
      gm->a[8+ioff12] = gm->a[8+ioff12] + zn*onePlusR*c2z + xn*c2x+ yn*c2y;

      gm->a[0+ioff13] = gm->a[0+ioff13] + xn*onePlusR*c3x + yn*c3y+ zn*c3z;
      gm->a[1+ioff13] = gm->a[1+ioff13] + xn*R*c3y;
      gm->a[2+ioff13] = gm->a[2+ioff13] + xn*R*c3z;
      gm->a[3+ioff13] = gm->a[3+ioff13] + yn*R*c3x;
      gm->a[4+ioff13] = gm->a[4+ioff13] + yn*onePlusR*c3y + zn*c3z+ xn*c3x;
      gm->a[5+ioff13] = gm->a[5+ioff13] + yn*R*c3z;
      gm->a[6+ioff13] = gm->a[6+ioff13] + zn*R*c3x;
      gm->a[7+ioff13] = gm->a[7+ioff13] + zn*R*c3y;
      gm->a[8+ioff13] = gm->a[8+ioff13] + zn*onePlusR*c3z + xn*c3x+ yn*c3y;

      gm->a[0+ioff14] = gm->a[0+ioff14] + xn*onePlusR*c4x + yn*c4y+ zn*c4z;
      gm->a[1+ioff14] = gm->a[1+ioff14] + xn*R*c4y;
      gm->a[2+ioff14] = gm->a[2+ioff14] + xn*R*c4z;
      gm->a[3+ioff14] = gm->a[3+ioff14] + yn*R*c4x;
      gm->a[4+ioff14] = gm->a[4+ioff14] + yn*onePlusR*c4y + zn*c4z+ xn*c4x;
      gm->a[5+ioff14] = gm->a[5+ioff14] + yn*R*c4z;
      gm->a[6+ioff14] = gm->a[6+ioff14] + zn*R*c4x;
      gm->a[7+ioff14] = gm->a[7+ioff14] + zn*R*c4y;
      gm->a[8+ioff14] = gm->a[8+ioff14] + zn*onePlusR*c4z + xn*c4x+ yn*c4y;

      /* node2 */
      xn = nx2/3.0;
      yn = ny2/3.0;
      zn = nz2/3.0;

      gm->a[0+idiag2] = gm->a[0+idiag2] + xn*onePlusR*c2x + yn*c2y+ zn*c2z;
      gm->a[1+idiag2] = gm->a[1+idiag2] + xn*R*c2y;
      gm->a[2+idiag2] = gm->a[2+idiag2] + xn*R*c2z;
      gm->a[3+idiag2] = gm->a[3+idiag2] + yn*R*c2x;
      gm->a[4+idiag2] = gm->a[4+idiag2] + yn*onePlusR*c2y + zn*c2z+ xn*c2x;
      gm->a[5+idiag2] = gm->a[5+idiag2] + yn*R*c2z;
      gm->a[6+idiag2] = gm->a[6+idiag2] + zn*R*c2x;
      gm->a[7+idiag2] = gm->a[7+idiag2] + zn*R*c2y;
      gm->a[8+idiag2] = gm->a[8+idiag2] + zn*onePlusR*c2z + xn*c2x+ yn*c2y;

      gm->a[0+ioff23] = gm->a[0+ioff23] + xn*onePlusR*c3x + yn*c3y+ zn*c3z;
      gm->a[1+ioff23] = gm->a[1+ioff23] + xn*R*c3y;
      gm->a[2+ioff23] = gm->a[2+ioff23] + xn*R*c3z;
      gm->a[3+ioff23] = gm->a[3+ioff23] + yn*R*c3x;
      gm->a[4+ioff23] = gm->a[4+ioff23] + yn*onePlusR*c3y + zn*c3z+ xn*c3x;
      gm->a[5+ioff23] = gm->a[5+ioff23] + yn*R*c3z;
      gm->a[6+ioff23] = gm->a[6+ioff23] + zn*R*c3x;
      gm->a[7+ioff23] = gm->a[7+ioff23] + zn*R*c3y;
      gm->a[8+ioff23] = gm->a[8+ioff23] + zn*onePlusR*c3z + xn*c3x+ yn*c3y;

      gm->a[0+ioff24] = gm->a[0+ioff24] + xn*onePlusR*c4x + yn*c4y+ zn*c4z;
      gm->a[1+ioff24] = gm->a[1+ioff24] + xn*R*c4y;
      gm->a[2+ioff24] = gm->a[2+ioff24] + xn*R*c4z;
      gm->a[3+ioff24] = gm->a[3+ioff24] + yn*R*c4x;
      gm->a[4+ioff24] = gm->a[4+ioff24] + yn*onePlusR*c4y + zn*c4z+ xn*c4x;
      gm->a[5+ioff24] = gm->a[5+ioff24] + yn*R*c4z;
      gm->a[6+ioff24] = gm->a[6+ioff24] + zn*R*c4x;
      gm->a[7+ioff24] = gm->a[7+ioff24] + zn*R*c4y;
      gm->a[8+ioff24] = gm->a[8+ioff24] + zn*onePlusR*c4z + xn*c4x+ yn*c4y;

      gm->a[0+ioff21] = gm->a[0+ioff21] + xn*onePlusR*c1x + yn*c1y+ zn*c1z;
      gm->a[1+ioff21] = gm->a[1+ioff21] + xn*R*c1y;
      gm->a[2+ioff21] = gm->a[2+ioff21] + xn*R*c1z;
      gm->a[3+ioff21] = gm->a[3+ioff21] + yn*R*c1x;
      gm->a[4+ioff21] = gm->a[4+ioff21] + yn*onePlusR*c1y + zn*c1z+ xn*c1x;
      gm->a[5+ioff21] = gm->a[5+ioff21] + yn*R*c1z;
      gm->a[6+ioff21] = gm->a[6+ioff21] + zn*R*c1x;
      gm->a[7+ioff21] = gm->a[7+ioff21] + zn*R*c1y;
      gm->a[8+ioff21] = gm->a[8+ioff21] + zn*onePlusR*c1z + xn*c1x+ yn*c1y;
    
      /* node3 */
      xn = nx3/3.0;
      yn = ny3/3.0;
      zn = nz3/3.0;

      gm->a[0+idiag3] = gm->a[0+idiag3] + xn*onePlusR*c3x + yn*c3y+ zn*c3z;
      gm->a[1+idiag3] = gm->a[1+idiag3] + xn*R*c3y;
      gm->a[2+idiag3] = gm->a[2+idiag3] + xn*R*c3z;
      gm->a[3+idiag3] = gm->a[3+idiag3] + yn*R*c3x;
      gm->a[4+idiag3] = gm->a[4+idiag3] + yn*onePlusR*c3y + zn*c3z+ xn*c3x;
      gm->a[5+idiag3] = gm->a[5+idiag3] + yn*R*c3z;
      gm->a[6+idiag3] = gm->a[6+idiag3] + zn*R*c3x;
      gm->a[7+idiag3] = gm->a[7+idiag3] + zn*R*c3y;
      gm->a[8+idiag3] = gm->a[8+idiag3] + zn*onePlusR*c3z + xn*c3x+ yn*c3y;

      gm->a[0+ioff34] = gm->a[0+ioff34] + xn*onePlusR*c4x + yn*c4y+ zn*c4z;
      gm->a[1+ioff34] = gm->a[1+ioff34] + xn*R*c4y;
      gm->a[2+ioff34] = gm->a[2+ioff34] + xn*R*c4z;
      gm->a[3+ioff34] = gm->a[3+ioff34] + yn*R*c4x;
      gm->a[4+ioff34] = gm->a[4+ioff34] + yn*onePlusR*c4y + zn*c4z+ xn*c4x;
      gm->a[5+ioff34] = gm->a[5+ioff34] + yn*R*c4z;
      gm->a[6+ioff34] = gm->a[6+ioff34] + zn*R*c4x;
      gm->a[7+ioff34] = gm->a[7+ioff34] + zn*R*c4y;
      gm->a[8+ioff34] = gm->a[8+ioff34] + zn*onePlusR*c4z + xn*c4x+ yn*c4y;

      gm->a[0+ioff31] = gm->a[0+ioff31] + xn*onePlusR*c1x + yn*c1y+ zn*c1z;
      gm->a[1+ioff31] = gm->a[1+ioff31] + xn*R*c1y;
      gm->a[2+ioff31] = gm->a[2+ioff31] + xn*R*c1z;
      gm->a[3+ioff31] = gm->a[3+ioff31] + yn*R*c1x;
      gm->a[4+ioff31] = gm->a[4+ioff31] + yn*onePlusR*c1y + zn*c1z+ xn*c1x;
      gm->a[5+ioff31] = gm->a[5+ioff31] + yn*R*c1z;
      gm->a[6+ioff31] = gm->a[6+ioff31] + zn*R*c1x;
      gm->a[7+ioff31] = gm->a[7+ioff31] + zn*R*c1y;
      gm->a[8+ioff31] = gm->a[8+ioff31] + zn*onePlusR*c1z + xn*c1x+ yn*c1y;

      gm->a[0+ioff32] = gm->a[0+ioff32] + xn*onePlusR*c2x + yn*c2y+ zn*c2z;
      gm->a[1+ioff32] = gm->a[1+ioff32] + xn*R*c2y;
      gm->a[2+ioff32] = gm->a[2+ioff32] + xn*R*c2z;
      gm->a[3+ioff32] = gm->a[3+ioff32] + yn*R*c2x;
      gm->a[4+ioff32] = gm->a[4+ioff32] + yn*onePlusR*c2y + zn*c2z+ xn*c2x;
      gm->a[5+ioff32] = gm->a[5+ioff32] + yn*R*c2z;
      gm->a[6+ioff32] = gm->a[6+ioff32] + zn*R*c2x;
      gm->a[7+ioff32] = gm->a[7+ioff32] + zn*R*c2y;
      gm->a[8+ioff32] = gm->a[8+ioff32] + zn*onePlusR*c2z + xn*c2x+ yn*c2y;

      /* node4 */
      xn = nx4/3.0;
      yn = ny4/3.0;
      zn = nz4/3.0;

      gm->a[0+idiag4] = gm->a[0+idiag4] + xn*onePlusR*c4x + yn*c4y+ zn*c4z;
      gm->a[1+idiag4] = gm->a[1+idiag4] + xn*R*c4y;
      gm->a[2+idiag4] = gm->a[2+idiag4] + xn*R*c4z;
      gm->a[3+idiag4] = gm->a[3+idiag4] + yn*R*c4x;
      gm->a[4+idiag4] = gm->a[4+idiag4] + yn*onePlusR*c4y + zn*c4z+ xn*c4x;
      gm->a[5+idiag4] = gm->a[5+idiag4] + yn*R*c4z;
      gm->a[6+idiag4] = gm->a[6+idiag4] + zn*R*c4x;
      gm->a[7+idiag4] = gm->a[7+idiag4] + zn*R*c4y;
      gm->a[8+idiag4] = gm->a[8+idiag4] + zn*onePlusR*c4z + xn*c4x+ yn*c4y;

      gm->a[0+ioff41] = gm->a[0+ioff41] + xn*onePlusR*c1x + yn*c1y+ zn*c1z;
      gm->a[1+ioff41] = gm->a[1+ioff41] + xn*R*c1y;
      gm->a[2+ioff41] = gm->a[2+ioff41] + xn*R*c1z;
      gm->a[3+ioff41] = gm->a[3+ioff41] + yn*R*c1x;
      gm->a[4+ioff41] = gm->a[4+ioff41] + yn*onePlusR*c1y + zn*c1z+ xn*c1x;
      gm->a[5+ioff41] = gm->a[5+ioff41] + yn*R*c1z;
      gm->a[6+ioff41] = gm->a[6+ioff41] + zn*R*c1x;
      gm->a[7+ioff41] = gm->a[7+ioff41] + zn*R*c1y;
      gm->a[8+ioff41] = gm->a[8+ioff41] + zn*onePlusR*c1z + xn*c1x+ yn*c1y;

      gm->a[0+ioff42] = gm->a[0+ioff42] + xn*onePlusR*c2x + yn*c2y+ zn*c2z;
      gm->a[1+ioff42] = gm->a[1+ioff42] + xn*R*c2y;
      gm->a[2+ioff42] = gm->a[2+ioff42] + xn*R*c2z;
      gm->a[3+ioff42] = gm->a[3+ioff42] + yn*R*c2x;
      gm->a[4+ioff42] = gm->a[4+ioff42] + yn*onePlusR*c2y + zn*c2z+ xn*c2x;
      gm->a[5+ioff42] = gm->a[5+ioff42] + yn*R*c2z;
      gm->a[6+ioff42] = gm->a[6+ioff42] + zn*R*c2x;
      gm->a[7+ioff42] = gm->a[7+ioff42] + zn*R*c2y;
      gm->a[8+ioff42] = gm->a[8+ioff42] + zn*onePlusR*c2z + xn*c2x+ yn*c2y;

      gm->a[0+ioff43] = gm->a[0+ioff43] + xn*onePlusR*c3x + yn*c3y+ zn*c3z;
      gm->a[1+ioff43] = gm->a[1+ioff43] + xn*R*c3y;
      gm->a[2+ioff43] = gm->a[2+ioff43] + xn*R*c3z;
      gm->a[3+ioff43] = gm->a[3+ioff43] + yn*R*c3x;
      gm->a[4+ioff43] = gm->a[4+ioff43] + yn*onePlusR*c3y + zn*c3z+ xn*c3x;
      gm->a[5+ioff43] = gm->a[5+ioff43] + yn*R*c3z;
      gm->a[6+ioff43] = gm->a[6+ioff43] + zn*R*c3x;
      gm->a[7+ioff43] = gm->a[7+ioff43] + zn*R*c3y;
      gm->a[8+ioff43] = gm->a[8+ioff43] + zn*onePlusR*c3z + xn*c3x+ yn*c3y;
    }
  }

  for(node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode(grid,node) &&
	 ( gridmoveSpecified(gm,node) || gridNodeGhost(grid,node) ) ) {
      for ( entry = gridmoveRowStart(gm, node) ;
	    entry < gridmoveRowStart(gm, node+1) ;
	    entry++ ) {
	for(i=0;i<9;i++) gm->a[i+9*entry] = 0.0;
      }
      entry = gridmoveRowEntry(gm,node,node);
      gm->a[0+9*entry] = 1.0;
      gm->a[4+9*entry] = 1.0;
      gm->a[8+9*entry] = 1.0;
    }
  }

  for(node=0;node<gridMaxNode(grid);node++) {
    entry = gridmoveRowEntry(gm,node,node);
    if (EMPTY != entry) {
      gridLU3x3( &gm->a[9*entry], &gm->lu[9*node] );
    } else {
      for(i=0;i<9;i++) gm->lu[i+9*node] = 0.0;
      gm->lu[0+9*node] = 1.0;
      gm->lu[4+9*node] = 1.0;
      gm->lu[8+9*node] = 1.0;
    }    
  }

  return gm;
}

GridMove *gridmoveElasticRelaxationDumpA(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node, entry, i;

  for(node=0;node<gridMaxNode(grid);node++) {
    for ( entry = gridmoveRowStart(gm, node) ;
	  entry < gridmoveRowStart(gm, node+1) ;
	  entry++ ) {
      for(i=0;i<9;i++) {
	printf(" node%6d off%6d entry%6d elem%2d %20.15f\n",
	       node, gridmoveRowNode(gm,entry), entry, i, gm->a[i+9*entry]);
      }
    }
  }

  return gm;
}

GridMove *gridmoveElasticRelaxationSubIteration(GridMove *gm, double *residual2)
{
  Grid *grid = gridmoveGrid(gm);
  int color, ncolor;
  int row, entry, col, i;
  double *b;
  double residual;

  b = malloc(3*gridMaxNode(grid)*sizeof(double));

  residual = 0;

  ncolor = 2;
  for(color = 0; color < ncolor ; color++) {

    for(row=color;row<gridMaxNode(grid);row+=ncolor) 
      for(i=0;i<3;i++) b[i+3*row] = 0.0;
 
    for(row=color;row<gridMaxNode(grid);row+=ncolor) {
      if ( gridValidNode(grid,row) &&
	   !gridmoveSpecified(gm,row) && 
	   gridNodeLocal(grid,row) ) {
	for ( entry = gridmoveRowStart(gm, row) ;
	      entry < gridmoveRowStart(gm, row+1) ;
	      entry++ ) {
	  col = gridmoveRowNode(gm,entry);
	  if (row != col) {
	    for(i=0;i<3;i++) {
	      b[i+3*row] -= 
		(   gm->a[i+0*3+9*entry]*gm->xyz[0+3*col]
	          + gm->a[i+1*3+9*entry]*gm->xyz[1+3*col]
	          + gm->a[i+2*3+9*entry]*gm->xyz[2+3*col] ) ;
	    }
	  }
	}
      }
    }

    for(row=color;row<gridMaxNode(grid);row+=ncolor) {
      if ( gridValidNode(grid,row) &&
	   !gridmoveSpecified(gm,row) && 
	   gridNodeLocal(grid,row) ) {
	gridBackSolve3x3(&gm->lu[9*row], &b[3*row]);
      }
    }

    for(row=color;row<gridMaxNode(grid);row+=ncolor) {
      if ( gridValidNode(grid,row) &&
	   !gridmoveSpecified(gm,row) && 
	   gridNodeLocal(grid,row) ) {
	for(i=0;i<3;i++) 
	  residual += (   (gm->xyz[i+3*row] - b[i+3*row])
	                * (gm->xyz[i+3*row] - b[i+3*row]) );
	for(i=0;i<3;i++) gm->xyz[i+3*row] = b[i+3*row];
      }
    }
  }

  *residual2 = residual;

  free(b);
  return gm;
}

GridMove *gridmoveElasticRelaxationShutDown(GridMove *gm)
{
  Grid *grid = gridmoveGrid(gm);
  int node, i;
  double xyz0[3];

  for(node=0;node<3*gridMaxNode(grid);node++)
    gm->lastxyz[node] += gm->xyz[node];

  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridValidNode(grid,node) && !gridmoveSpecified(gm,node)) {
      gridNodeXYZ(grid,node,xyz0);
      for(i=0;i<3;i++) 
	gm->displacement[i+3*node] = gm->lastxyz[i+3*node] - xyz0[i];
    }
  }

  free(gm->rowStart); gm->rowStart=NULL;
  free(gm->compRow);  gm->compRow=NULL;

  free(gm->a);        gm->a=NULL;
  free(gm->xyz);      gm->xyz=NULL;
  free(gm->lastxyz);     gm->lastxyz=NULL;

  return gm;
}

GridMove *gridmoveElasticRelaxation(GridMove *gm, int nsteps, int subIterations)
{
  int step, iteration;
  double position;
  double rmsResidual;

  if (gm != gridmoveElasticRelaxationStartUp(gm)) return NULL;
  for(step=0;step<nsteps;step++) {
    position = (double)(step+1)/(double)nsteps;
    gridmoveElasticRelaxationStartStep(gm, position);    
    for(iteration=0;iteration<subIterations;iteration++) {
      gridmoveElasticRelaxationSubIteration(gm, &rmsResidual);
      /* printf("Iteration %4d Residual %23.15e\n",iteration,rmsResidual); */
    }
  }
  gridmoveElasticRelaxationShutDown(gm);

  return gm;
}

GridMove *gridmoveRelaxationStartUp(GridMove *gm, int relaxationScheme )
{
  switch (relaxationScheme) {
  case gridmoveELASTIC_SCHEME:
    gm->relaxationScheme = relaxationScheme;
    return gridmoveElasticRelaxationStartUp(gm);
    break;
  case gridmoveSPRING_SCHEME:
    gm->relaxationScheme = relaxationScheme;
    return gridmoveSpringRelaxationStartUp(gm);
    break;
  default:
    printf("%s: %d: ERROR relaxationScheme %d unknown\n",
	   __FILE__,__LINE__,relaxationScheme);
  }
  return NULL;
}

GridMove *gridmoveRelaxationStartStep(GridMove *gm, double position)
{
  switch (gm->relaxationScheme) {
  case gridmoveELASTIC_SCHEME:
    return gridmoveElasticRelaxationStartStep(gm, position);
    break;
  case gridmoveSPRING_SCHEME:
    return gridmoveSpringRelaxationStartStep(gm, position);
    break;
  default:
    printf("%s: %d: ERROR relaxationScheme %d unknown.\n",
	   __FILE__,__LINE__,gm->relaxationScheme);
    printf("call gridmoveRelaxationStartUp first.\n");
  }
  return NULL;
}

GridMove *gridmoveRelaxationSubIteration(GridMove *gm, double *residual2)
{
  switch (gm->relaxationScheme) {
  case gridmoveELASTIC_SCHEME:
    return gridmoveElasticRelaxationSubIteration(gm, residual2);
    break;
  case gridmoveSPRING_SCHEME:
    return gridmoveSpringRelaxationSubIteration(gm, residual2);
    break;
  default:
    printf("%s: %d: ERROR relaxationScheme %d unknown.\n",
	   __FILE__,__LINE__,gm->relaxationScheme);
    printf("call gridmoveRelaxationStartUp first.\n");
  }
  return NULL;
}

GridMove *gridmoveRelaxationShutDown(GridMove *gm)
{
  switch (gm->relaxationScheme) {
  case gridmoveELASTIC_SCHEME:
    gm->relaxationScheme = gridmoveEMPTY_SCHEME;
    return gridmoveElasticRelaxationShutDown(gm);
    break;
  case gridmoveSPRING_SCHEME:
    gm->relaxationScheme = gridmoveEMPTY_SCHEME;
    return gridmoveSpringRelaxationShutDown(gm);
    break;
  default:
    printf("%s: %d: ERROR relaxationScheme %d unknown.\n",
	   __FILE__,__LINE__,gm->relaxationScheme);
    printf("call gridmoveRelaxationStartUp first.\n");
  }
  return NULL;
}
