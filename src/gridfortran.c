
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

void gridcreate_( int *partId, int *nnode, double *x, double *y, double *z ,
		  int *ncell, int *maxcell, int *c2n )
{
  int node, cell;
  grid = gridCreate( *nnode, *ncell, 5000, 0);
  gridSetPartId(grid, *partId );
  queue = queueCreate( 9 ); /* 3:xyz + 6:m */
  for ( node=0; node<*nnode; node++) gridAddNode(grid,x[node],y[node],z[node]);
  for ( cell=0; cell<*ncell; cell++) gridAddCell( grid,
						  c2n[cell+0*(*maxcell)] - 1,
						  c2n[cell+1*(*maxcell)] - 1,
						  c2n[cell+2*(*maxcell)] - 1,
						  c2n[cell+3*(*maxcell)] - 1 );
#ifdef PARALLEL_VERBOSE 
  printf(" %6d populated                nnode%9d ncell%9d AR%14.10f\n",
	 gridPartId(grid),gridNNode(grid),gridNCell(grid),gridMinAR(grid));
  fflush(stdout);
#endif
}

void gridfree_( )
{
  queueFree(queue);
  gridFree(grid);
}

void gridinsertboundary_( int *faceId, int *nnode, int *nodedim, int *inode, 
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

void gridsetmap_( int *nnode, double* map )
{
  int node;
  for ( node=0; node<*nnode; node++) 
    gridSetMap( grid, node,
		map[0+6*node], map[1+6*node], map[2+6*node],
		map[3+6*node], map[4+6*node], map[5+6*node] );
#ifdef PARALLEL_VERBOSE 
  printf(" %6d applied metric                                         AR%14.10f\n",
	 gridPartId(grid),gridMinAR(grid));
#endif
}

void gridsetnodepart_( int *nnode, int *part )
{
  int node;
  for ( node=0; node<*nnode; node++) gridSetNodePart(grid, node, part[node]-1);
}

void gridsetnodelocal2global_( int *partId, int *nnodeg, 
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
  }
}

void gridsetcelllocal2global_( int *ncellg, int *ncell, int *local2global )
{
  int cell;
  gridSetGlobalNCell(grid, *ncellg);
  for ( cell=0; cell<*ncell; cell++){ 
    gridSetCellGlobal(grid, cell, local2global[cell]-1);
  }
}

void gridfreezenode_( int *node )
{
  gridFreezeNode( grid, *node );
}

void gridprojectallfaces_( )
{
  int face, node, nodes[3], faceId;
  double ar0, ar1;

  ar0 = gridMinAR(grid);
  for(face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      for(node=0;node<3;node++) {
	if ( !gridNodeFrozen(grid,nodes[node])  ) {
	  if (grid != gridProjectNodeToFace(grid, nodes[node], faceId ) )
	    printf( "ERROR: %s: %d: project failed.%d\n",
		    __FILE__,__LINE__,gridPartId(grid) );

	}
      }
    }
  }
  ar1 = gridMinAR(grid);

#ifdef PARALLEL_VERBOSE 
  printf( " %6d project faces           initial AR%14.10f final AR%14.10f\n",
	  gridPartId(grid),ar0,ar1 );
  fflush(stdout);
#endif
}

void gridminar_( double *aspectratio )
{
  *aspectratio = gridMinAR( grid );
}

void gridwritetecplotsurfacezone_( )
{
  char filename[256];
  double ar0, ar1;
  ar0 = gridMinAR(grid);
  sprintf(filename, "grid%03d.t", gridPartId(grid)+1 );
  gridWriteTecplotSurfaceZone(grid,filename);
  ar1 = gridMinAR(grid);

#ifdef PARALLEL_VERBOSE 
  printf( " %6d tecplot dump            initial AR%14.10f final AR%14.10f\n",
	  gridPartId(grid),ar0,ar1 );
  fflush(stdout);
#endif
}

void gridexportfast_( )
{
  char filename[256];
  sprintf(filename, "grid%03d.fgrid", gridPartId(grid)+1 );
  gridExportFAST(grid,filename);
}

void gridparallelswap_( int *processor, double *ARlimit )
{
#ifdef PARALLEL_VERBOSE 
  printf(" %6d swap  processor %2d      initial AR%14.10f",
	 gridPartId(grid),*processor,gridMinAR(grid));
#endif
  if (*processor == -1) {
    gridParallelSwap(grid,NULL,*ARlimit);
  } else {
    gridParallelSwap(grid,queue,*ARlimit);
  } 
}

void gridparallelsmoothfaceinterior_( int *processor )
{
  GridBool localOnly;
  localOnly = (-1 == (*processor));
  gridSmoothFaceInterior(grid, localOnly );
#ifdef PARALLEL_VERBOSE 
  if (localOnly) {
    printf( " %6d smooth volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"local only        ",gridMinAR(grid) );
  } else {
    printf( " %6d smooth volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"near ghost only   ",gridMinAR(grid) );
  }
  fflush(stdout);
#endif
}

void gridparalleladaptwithoutcad_( int *processor, 
				   double *minLength, double *maxLength )
{
#ifdef PARALLEL_VERBOSE 
  printf(" %6d adapt processor %2d ",gridPartId(grid),*processor);
#endif
  if (*processor == -1) {
    gridParallelAdaptWithOutCAD(grid,NULL,*minLength, *maxLength);
  } else {
    gridParallelAdaptWithOutCAD(grid,queue,*minLength, *maxLength);
  } 
}

void queuedumpsize_( int *nInt, int *nDouble )
{
  queueDumpSize(queue, nInt, nDouble);
}

void queuedump_( int *nInt, int *nDouble, int *ints, double *doubles )
{
  queueDump(queue, ints, doubles);
}

void gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles )
{
  queueLoad(queue, ints, doubles);
  gridApplyQueue(grid,queue);
  queueReset(queue);
}

void gridsize_( int *nnodeg, int *ncellg )
{
  *nnodeg = gridGlobalNNode(grid);
  *ncellg = gridGlobalNCell(grid);
}

void gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset,
		       int *oldncellg, int *newncellg, int *celloffset )
{
  gridGlobalShiftNode( grid, *oldnnodeg, *newnnodeg, *nodeoffset);
  gridGlobalShiftCell( grid, *oldncellg, *newncellg, *celloffset);
}

void gridnunusednodeglobal_( int *nunused )
{
  *nunused = gridNUnusedNodeGlobal( grid );
}

void gridgetunusednodeglobal_( int *nunused, int *unused )
{
  gridGetUnusedNodeGlobal( grid, unused );
}

void gridjoinunusednodeglobal_( int *nunused, int *unused )
{
  int i;
  for (i=0;i<(*nunused);i++) gridJoinUnusedNodeGlobal( grid, unused[i] );
}

void grideliminateunusednodeglobal_(  )
{
  gridEliminateUnusedNodeGlobal( grid );
}

void gridnunusedcellglobal_( int *nunused )
{
  *nunused = gridNUnusedCellGlobal( grid );
}

void gridgetunusedcellglobal_( int *nunused, int *unused )
{
  gridGetUnusedCellGlobal( grid, unused );
}

void gridjoinunusedcellglobal_( int *nunused, int *unused )
{
  int i;
  for (i=0;i<(*nunused);i++) gridJoinUnusedCellGlobal( grid, unused[i] );
}

void grideliminateunusedcellglobal_(  )
{
  gridEliminateUnusedCellGlobal( grid );
}

void gridsortfun3d_( int *nnodes0, int *nnodes01, int *nnodesg, 
		    int *ncell, int *ncellg )
{
  gridSortNodeFUN3D( grid, nnodes0 );
  *nnodes01 = gridNNode(grid);
  *nnodesg = gridGlobalNNode(grid);
  *ncell = gridNCell(grid);
  *ncellg = gridGlobalNCell(grid);
}

void gridgetnodes_( int *nnode, int *l2g, double *x, double *y, double *z)
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

void gridgetcell_( int *cell, int *nodes, int *global )
{
  gridCell(grid,(*cell)-1,nodes);
  nodes[0]++;
  nodes[1]++;
  nodes[2]++;
  nodes[3]++;
  *global = gridCellGlobal(grid,(*cell)-1)+1;
}

void gridgetbcsize_( int *ibound, int *nface )
{
  int face, nodes[3], id;
  
  *nface = 0;
  for (face=0;face<gridMaxFace(grid);face++) {
    if ( grid == gridFace(grid,face,nodes,&id) ) {
      if ( *ibound == id ) (*nface)++;
    }
  }
}

void gridgetbc_( int *ibound, int *nface, int *ndim, int *f2n )
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

void gridsetnaux_( int *naux )
{
  gridSetNAux(grid, *naux);
  queueFree( queue );
  queue = queueCreate( 9 + gridNAux(grid) ); /* 3:xyz + 6:m + naux */
}

void gridsetauxvector_( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    gridSetAux(grid,node,(*offset),x[node]);
  }
}

void gridsetauxmatrix_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

void gridsetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

void gridgetauxvector_( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    x[node] = gridAux(grid,node,(*offset));
  }
}

void gridgetauxmatrix_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

void gridgetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

void gridghostcount_( int *nproc, int *count )
{
  int node, faces;
  for(node=0;node<(*nproc);node++) count[node] = 0;
  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) { 
      count[gridNodePart(grid,node)]++;
      faces = gridNodeFaceIdDegree(grid,node);
      if (faces>0) {
	count[gridNodePart(grid,node)] += (faces+1);
      }
    }
  }
}

void gridloadghostnodes_( int *nproc, int *clientindex,
			  int *clientsize, int *localnode, int *globalnode )
{
  int node, part;
  int *count;
  int face, ids, id[MAXFACEIDDEG];

  count = malloc( (*nproc) * sizeof(int) );

  for(node=0;node<(*nproc);node++) count[node] = 0;
  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) {
      part = gridNodePart(grid,node);
      localnode[ count[part]+clientindex[part]-1] = node+1;
      globalnode[count[part]+clientindex[part]-1] = gridNodeGlobal(grid,node)+1;
      count[part]++;
      gridNodeFaceId(grid, node, MAXFACEIDDEG, &ids, id );
      if (ids>0) {
	localnode[ count[part]+clientindex[part]-1] = -ids;
	globalnode[count[part]+clientindex[part]-1] = -ids;
	count[part]++;
	for (face=0;face<ids;face++) {
	  localnode[ count[part]+clientindex[part]-1] = id[face];
	  globalnode[count[part]+clientindex[part]-1] = id[face];
	  count[part]++;
	}
      }
    }
  }
  free(count);
}

void gridloadglobalnodedata_( int *ndim, int *nnode, int *nodes, double *data )
{
  int node, localnode, face, ids, faceId;
  double uv[2], xyz[3];

  localnode=0;
  node=0;
  while (node<(*nnode)) {
    if (nodes[node] > 0) {
      localnode = gridGlobal2Local(grid, nodes[node]-1);
      if (grid != gridNodeXYZ(grid,localnode,xyz)) 
	printf("%d: ERROR: %s: %d: get xyz from invalid node local %d global %d.\n",
	       gridPartId(grid),__FILE__, __LINE__, localnode, nodes[node]-1);
      data[0+(*ndim)*node] = xyz[0];
      data[1+(*ndim)*node] = xyz[1];
      data[2+(*ndim)*node] = xyz[2];
      node++;
    } else {
      ids = -nodes[node];
      node++;
      for(face=0;face<ids;face++) {
	faceId = nodes[node];
	gridNodeUV(grid,localnode,faceId,uv);
	data[0+(*ndim)*node] = uv[0];
	data[1+(*ndim)*node] = uv[1];
	node++;
      }
    } 
  }
}

void gridsetlocalnodedata_( int *ndim, int *nnode, int *nodes, double *data )
{
  int node, localnode;
  int face, ids, faceId;

  localnode=0;
  node=0;
  while (node<(*nnode)) {
    if (nodes[node] > 0) {
      localnode = nodes[node]-1;
      if (grid != gridSetNodeXYZ(grid,localnode,&data[(*ndim)*node]) )
	printf("ERROR: %s: %d: set invalid node %d .\n",
	       __FILE__, __LINE__, nodes[node]-1);
      node++;
    } else {
      ids = -nodes[node];
      node++;
      for(face=0;face<ids;face++) {
	faceId = nodes[node];
	gridSetNodeUV(grid,localnode,faceId,
		      data[0+(*ndim)*node],data[1+(*ndim)*node]);
	node++;
      }
    } 
  }

#ifdef PARALLEL_VERBOSE 
  printf( " %6d update xfer                      %s    AR%14.10f\n",
	  gridPartId(grid),"                  ",gridMinAR(grid) );
  fflush(stdout);
#endif
}
void gridcopyabouty0_( int *symmetryFaceId, int *mirrorAux )
{
  gridCopyAboutY0(grid, *symmetryFaceId, *mirrorAux-1 );
}
