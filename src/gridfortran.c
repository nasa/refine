
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridfortran.h"
#include "grid.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "queue.h"
#include "gridmpi.h"

static Grid *grid;
static Queue *queue;

int gridcreate_( int *partId, int *nnode, double *x, double *y, double *z ,
		 int *ncell, int *maxcell, int *c2n )
{
  int node, cell;
  int nodes[4];
  double xyz[3];
  grid = gridCreate( *nnode, *ncell, 5000, 0);
  gridSetPartId(grid, *partId );
  queue = queueCreate( 9 ); /* 3:xyz + 6:m */
  for ( node=0; node<*nnode; node++) gridAddNode(grid,x[node],y[node],z[node]);
  for ( cell=0; cell<*ncell; cell++) gridAddCell( grid,
						  c2n[cell+0*(*maxcell)] - 1,
						  c2n[cell+1*(*maxcell)] - 1,
						  c2n[cell+2*(*maxcell)] - 1,
						  c2n[cell+3*(*maxcell)] - 1 );
  printf(" %6d populated                nnode%9d ncell%9d AR%14.10f\n",
	 gridPartId(grid),gridNNode(grid),gridNCell(grid),gridMinAR(grid));
  fflush(stdout);
}

int gridfree_( )
{
  queueFree(queue);
  gridFree(grid);
}

int gridinsertboundary_( int *faceId, int *nnode, int *nodedim, int *inode, 
			 int *nface, int *dim1, int *dim2, int *f2n )
{
  int face;
  int node0, node1, node2;
  for(face=0;face<*nface;face++){
    node0 = f2n[face+0*(*dim1)] - 1;
    node1 = f2n[face+1*(*dim1)] - 1;
    node2 = f2n[face+2*(*dim1)] - 1;
    node0 = inode[node0] - 1;
    node1 = inode[node1] - 1;
    node2 = inode[node2] - 1;
    if ( gridNodeGhost(grid,node0) && 
	 gridNodeGhost(grid,node1) && 
	 gridNodeGhost(grid,node2) ) {
    }else{
      gridAddFace(grid, node0, node1, node2, *faceId);
    }
  }
}

int gridsetmap_( int *nnode, double* map )
{
  int node;
  for ( node=0; node<*nnode; node++) 
    gridSetMap( grid, node,
		map[0+6*node], map[1+6*node], map[2+6*node],
		map[3+6*node], map[4+6*node], map[5+6*node] );
  printf(" %6d applied metric                                         AR%14.10f\n",
	 gridPartId(grid),gridMinAR(grid));
}

int gridsetnodepart_( int *nnode, int *part )
{
  int node;
  for ( node=0; node<*nnode; node++) gridSetNodePart(grid, node, part[node]-1);
}

int gridsetnodelocal2global_( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global )
{
  int node;
  gridSetPartId(grid, *partId );
  gridSetGlobalNNode(grid, *nnodeg );
  for ( node=0; node<*nnode; node++){ 
    gridSetNodeGlobal(grid, node, local2global[node]-1);
    if ( node < *nnode0 ) {
      gridSetNodePart(grid, node, *partId );
    }else{
      gridSetNodePart(grid, node, EMPTY );
    }
    /*printf("%d l2g node %d global %d part %d\n", gridPartId(grid),node,
      gridNodeGlobal(grid,node),gridNodePart(grid, node)); */
  }
}

int gridsetcelllocal2global_( int *ncellg, int *ncell, int *local2global )
{
  int cell;
  gridSetGlobalNCell(grid, *ncellg);
  for ( cell=0; cell<*ncell; cell++){ 
    gridSetCellGlobal(grid, cell, local2global[cell]-1);
  }
}

int gridprojectallfaces_( )
{
  int face, nodes[3], faceId;
  for(face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      gridProjectNodeToFace(grid, nodes[0], faceId );
      gridProjectNodeToFace(grid, nodes[1], faceId );
      gridProjectNodeToFace(grid, nodes[2], faceId );
    }
  }
}

int gridswap_( )
{
  gridSwap(grid);
  printf(" post swap min AR %17.15f\n",gridMinAR(grid));
}

int gridsmoothvolume_( )
{
  gridSmoothVolume(grid);
  printf( " %6d smooth volume                    %s    AR%14.10f\n",
	  gridPartId(grid),"                  ",gridMinAR(grid) );
  fflush(stdout);
}

int gridsmoothfaceinterior_( int *processor )
{
  bool localOnly;
  localOnly = (-1 == (*processor));
  gridSmoothFaceInterior(grid, localOnly );
  if (localOnly) {
    printf( " %6d smooth volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"local only        ",gridMinAR(grid) );
  } else {
    printf( " %6d smooth volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"near ghost only   ",gridMinAR(grid) );
  }
  fflush(stdout);
}

int gridadaptwithoutcad_( double *minLength, double *maxLength )
{
  gridAdaptWithOutCAD(grid,*minLength, *maxLength);
}

int gridwritetecplotsurfacezone_( )
{
  char filename[256];
  sprintf(filename, "grid%03d.t", gridPartId(grid)+1 );
  gridWriteTecplotSurfaceZone(grid,filename);
}

int gridexportfast_( )
{
  char filename[256];
  sprintf(filename, "grid%03d.fgrid", gridPartId(grid)+1 );
  gridExportFAST(grid,filename);
  
}

int gridparalleladaptwithoutcad_( int *processor, 
				  double *minLength, double *maxLength )
{
  printf(" %6d adapt processor %2d ",gridPartId(grid),*processor);
  if (*processor == -1) {
    gridParallelAdaptWithOutCAD(grid,NULL,*minLength, *maxLength);
  } else {
    gridParallelAdaptWithOutCAD(grid,queue,*minLength, *maxLength);
  } 
}

int gridparallelswap_( int *processor )
{
  printf(" %6d swap  processor %2d ",gridPartId(grid),*processor);
  if (*processor == -1) {
    gridParallelSwap(grid,NULL);
  } else {
    gridParallelSwap(grid,queue);
  } 
}

int queuedumpsize_( int *nInt, int *nDouble )
{
  queueDumpSize(queue, nInt, nDouble);
}

int queuedump_( int *nInt, int *nDouble, int *ints, double *doubles )
{
  queueDump(queue, ints, doubles);
}

int gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles )
{
  queueLoad(queue, ints, doubles);
  gridApplyQueue(grid,queue);
  queueReset(queue);
}

int gridsize_( int *nnodeg, int *ncellg )
{
  *nnodeg = gridGlobalNNode(grid);
  *ncellg = gridGlobalNCell(grid);
}

int gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset,
		      int *oldncellg, int *newncellg, int *celloffset )
{
  gridGlobalShiftNode( grid, *oldnnodeg, *newnnodeg, *nodeoffset);
  gridGlobalShiftCell( grid, *oldncellg, *newncellg, *celloffset);
}

int gridnunusedcellglobal_( int *nunused )
{
  *nunused = gridNUnusedCellGlobal( grid );
}

int gridgetunusedcellglobal_( int *nunused, int *unused )
{
  gridGetUnusedCellGlobal( grid, unused );
}

int gridjoinunusedcellglobal_( int *nunused, int *unused )
{
  int i;
  for (i=0;i<(*nunused);i++) gridJoinUnusedCellGlobal( grid, unused[i] );
}

int grideliminateunusedcellglobal_(  )
{
  gridEliminateUnusedCellGlobal( grid );
}

int gridsortfun3d_( int *nnodes0, int *nnodes01, int *nnodesg, 
		    int *ncell, int *ncellg )
{
  gridSortNodeFUN3D( grid, nnodes0 );
  *nnodes01 = gridNNode(grid);
  *nnodesg = gridGlobalNNode(grid);
  *ncell = gridNCell(grid);
  *ncellg = gridGlobalNCell(grid);
}

int gridgetnodes_( int *nnode, int *l2g, double *x, double *y, double *z)
{
  int node;
  double xyz[3];
  for (node=0;node<gridNNode(grid);node++) {
    l2g[node] = gridNodeGlobal(grid,node)+1;
    gridNodeXYZ(grid,node,xyz);
    x[node] = xyz[0];
    y[node] = xyz[1];
    z[node] = xyz[2];
  }
}

int gridgetcell_( int *cell, int *nodes, int *global )
{
  gridCell(grid,(*cell)-1,nodes);
  nodes[0]++;
  nodes[1]++;
  nodes[2]++;
  nodes[3]++;
  *global = gridCellGlobal(grid,(*cell)-1)+1;
}

int gridgetbcsize_( int *ibound, int *nface )
{
  int face, nodes[3], id;
  
  *nface = 0;
  for (face=0;face<gridMaxFace(grid);face++) {
    if ( grid == gridFace(grid,face,nodes,&id) ) {
      if ( *ibound == id ) (*nface)++;
    }
  }
}

int gridgetbc_( int *ibound, int *nface, int *ndim, int *f2n )
{
  int face, n, nodes[3], id;
  
  n = 0;
  for (face=0;face<gridMaxFace(grid);face++) {
    if ( grid == gridFace(grid,face,nodes,&id) ) {
      if ( *ibound == id ) {
	f2n[n+(*nface)*0] = nodes[0]+1;
	f2n[n+(*nface)*1] = nodes[1]+1;
	f2n[n+(*nface)*2] = nodes[2]+1;
	n++;
      }
    }
  }
}

int gridsetnaux_( int *naux )
{
  gridSetNAux(grid, *naux);
  queueFree( queue );
  queue = queueCreate( 9 + gridNAux(grid) ); /* 3:xyz + 6:m + naux */
}

int gridsetauxvector_( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    gridSetAux(grid,node,(*offset),x[node]);
  }
}

int gridsetauxmatrix_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

int gridsetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

int gridgetauxvector_( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    x[node] = gridAux(grid,node,(*offset));
  }
}

int gridgetauxmatrix_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

int gridgetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

int gridghostcount_( int *nproc, int *count )
{
  int node;
  for(node=0;node<(*nproc);node++) count[node] = 0;
  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) count[gridNodePart(grid,node)]++;
  }
}

int gridloadghostnodes_( int *nproc, int *clientindex,
			 int *clientsize, int *localnode, int *globalnode )
{
  int node, part;
  int *count;

  count = malloc( (*nproc) * sizeof(int) );

  for(node=0;node<(*nproc);node++) count[node] = 0;
  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) {
      part = gridNodePart(grid,node);
      localnode[ count[part]+clientindex[part]-1] = node+1;
      globalnode[count[part]+clientindex[part]-1] = gridNodeGlobal(grid,node)+1;
      count[part]++;
    }
  }
  free(count);
}

int gridloadglobalnodedata_( int *ndim, int *nnode, int *nodes, double *data )
{
  int node, localnode;
  double xyz[3];

  for(node=0;node<(*nnode);node++) {
    localnode = gridGlobal2Local(grid, nodes[node]-1);
    if (grid != gridNodeXYZ(grid,localnode,xyz)) 
      printf("ERROR: %s: %d: get xyz from invalid node local %d global &d.\n",
	     __FILE__, __LINE__, localnode, nodes[node]-1);
    data[0+(*ndim)*node] = xyz[0];
    data[1+(*ndim)*node] = xyz[1];
    data[2+(*ndim)*node] = xyz[2];
  }
}

int gridsetlocalnodedata_( int *ndim, int *nnode, int *nodes, double *data )
{
  int node;

  for(node=0;node<(*nnode);node++) { 
    if (grid != gridSetNodeXYZ(grid,nodes[node]-1,&data[(*ndim)*node]) )
      printf("ERROR: %s: %d: set invalid node %d .\n",
	     __FILE__, __LINE__, nodes[node]-1);
  }
  printf( " %6d update xfer                      %s    AR%14.10f\n",
	  gridPartId(grid),"                  ",gridMinAR(grid) );
  fflush(stdout);
}
