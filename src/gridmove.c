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
#include "gridmove.h"
#include "gridmath.h"

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

  return gm;
}

Grid *gridmoveGrid(GridMove *gm)
{
  return gm->grid;
}

void gridmoveFree(GridMove *gm)
{
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

GridMove *gridmoveSpringRelaxation(GridMove *gm, int nsteps, int subIterations)
{
  Grid *grid = gridmoveGrid(gm);
  int nsprings, *springs;
  double xyz0[3], dxyz[3], res[3];
  int s, node;
  int i, n0, n1;
  double *k;
  double *xyz, *source, *ksum, *kxyz;
  int step, iteration;
  double stepSize;
  double residual, count;

  int *c2e;

  c2e = malloc(6*gridMaxCell(grid)*sizeof(int));
  if (gm != gridmoveComputeC2E(gm, &nsprings, c2e) ) { free(c2e); return NULL; }
  springs = malloc(2*nsprings*sizeof(int));
  if (gm != gridmoveComputeSpringsWithC2E(gm, c2e, springs)) {
    free(c2e); 
    free(springs); 
    return NULL;
  }
  
  k = malloc(nsprings*sizeof(double));

  xyz = malloc(3*gridMaxNode(grid)*sizeof(double));
  source = malloc(3*gridMaxNode(grid)*sizeof(double));
  ksum = malloc(gridMaxNode(grid)*sizeof(double));
  kxyz = malloc(3*gridMaxNode(grid)*sizeof(double));
 
  for(node=0;node<gridMaxNode(grid);node++) gridNodeXYZ(grid,node,&xyz[3*node]);

  stepSize = 1.0 / (double)nsteps;
  for(step=0;step<nsteps;step++) {

    for(s=0;s<nsprings;s++) {
      n0 = springs[0+2*s];
      n1 = springs[1+2*s];
      gridSubtractVector(&xyz[3*n1],&xyz[3*n0],dxyz);
      k[s]= 1.0/gridDotProduct(dxyz,dxyz);
    }

    for(node=0;node<3*gridMaxNode(grid);node++) source[node]=0.0;
    for(s=0;s<nsprings;s++) {
      n0 = springs[0+2*s];
      n1 = springs[1+2*s];
      for(i=0;i<3;i++) {
	res[i] = k[s] * ( xyz[i+3*n0] - xyz[i+3*n1] );
	source[i+3*n0] += res[i];
	source[i+3*n1] -= res[i];
      }
    }

    for(node=0;node<gridMaxNode(grid);node++)
      if (gridValidNode(grid,node) && gridmoveSpecified(gm,node))
	for(i=0;i<3;i++) xyz[i+3*node] += stepSize*gm->displacement[i+3*node];

    for(iteration=0;iteration<subIterations;iteration++) {

      for(node=0;node<gridMaxNode(grid);node++) ksum[node]=0.0;
      for(node=0;node<3*gridMaxNode(grid);node++) kxyz[node]=0.0;
      for(s=0;s<nsprings;s++) {
	n0 = springs[0+2*s];
	n1 = springs[1+2*s];
	ksum[n0] += k[s];
	ksum[n1] += k[s];
	for(i=0;i<3;i++) {
	  kxyz[i+3*n0] += k[s] * xyz[i+3*n1];
	  kxyz[i+3*n1] += k[s] * xyz[i+3*n0];
	}
      }

      residual = 0.0; count =0.0;
      for(node=0;node<gridMaxNode(grid);node++)
	if (gridValidNode(grid,node) && !gridmoveSpecified(gm,node)) {
	  for(i=0;i<3;i++) {
	    res[i] = xyz[i+3*node] - 
	      ( kxyz[i+3*node] + source[i+3*node] ) / ksum[node];
	    xyz[i+3*node] = (  kxyz[i+3*node] + source[i+3*node] ) / ksum[node];
	  }
	  residual += gridDotProduct(res,res);
	  count += 1.0;
	}
      /*
      printf("Iteration %4d Residual %23.15e\n",iteration,sqrt(residual/count));
      */ 
   }
  }
  
  for(node=0;node<gridMaxNode(grid);node++)
    if (gridValidNode(grid,node) && !gridmoveSpecified(gm,node)) {
      gridNodeXYZ(grid,node,xyz0);
      for(i=0;i<3;i++) gm->displacement[i+3*node] = xyz[i+3*node] - xyz0[i];
    }

  free(k);

  free(xyz);
  free(source);
  free(ksum);
  free(kxyz);

  free(springs);
  free(c2e);

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
    if (grid == gridNodeXYZ(grid,node,xyz)) {
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
    if (gridGeometryFace( grid, node )) {
      gridNodeProjectionDisplacement(grid,node,displacement);
      gridmoveDisplace(gm,node,displacement);
    }
  }

  return gm;
}
