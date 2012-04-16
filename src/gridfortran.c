
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gridfortran.h"
#include "plan.h"
#include "grid.h"
#include "gridmove.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "queue.h"
#include "gridmpi.h"
#include "gridgeom.h"

/* #define TRAPFPE 1 */
#ifdef TRAPFPE
#define _GNU_SOURCE 1
#include <fenv.h>
#endif

static Grid *grid = NULL;
static GridMove *gm = NULL;
static Queue *queue = NULL;
static Plan *plan = NULL;

void FC_FUNC(gridapiversion,GRIDAPIVERSION)( int *refine_api_version )
{
  *refine_api_version = 100800000;
}

void FC_FUNC(gridcreate,GRIDCREATE)( int *partId, int *nnode, double *x, double *y, double *z )
{
  int node;

#ifdef TRAPFPE
  /* Enable some exceptions.  At startup all exceptions are masked.  */
  fprintf(stderr,"\ngridcreate_: Previous Exceptions %x\n",fetestexcept(FE_ALL_EXCEPT));
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  feclearexcept(FE_ALL_EXCEPT);
  fprintf(stderr,"gridcreate_: Floating Point Exception Handling On!\n");
  fprintf(stderr,"gridcreate_: Present  Exceptions %x\n\n",fetestexcept(FE_ALL_EXCEPT));
#endif

  grid = gridCreate( *nnode, 50000, 5000, 0);
  gridSetPartId(grid, *partId );
  gridSetCostConstraint(grid, gridCOST_CNST_VOLUME);
  queue = queueCreate( 9 ); /* 3:xyz + 6:m */
  for ( node=0; node<*nnode; node++) gridAddNode(grid,x[node],y[node],z[node]);
#ifdef PARALLEL_VERBOSE 
  printf(" %6d populated                nnode%9d ncell%9d AR%14.10f\n",
	 gridPartId(grid),gridNNode(grid),gridNCell(grid),gridMinAR(grid));
  fflush(stdout);
#endif
}

void FC_FUNC(gridfree,GRIDFREE)( void )
{
  queueFree(queue); queue = NULL;
  gridFree(grid); grid = NULL;
}

void FC_FUNC(gridinsertcells,GRIDINSERTCELLS)( int *nodes_per_cell, int *ncell, int *c2n )
{
  int cell;
 
  switch (*nodes_per_cell) {
  case 4:
    for ( cell=0; cell<*ncell; cell++) gridAddCell( grid,
						    c2n[0+4*cell] - 1,
						    c2n[1+4*cell] - 1,
						    c2n[2+4*cell] - 1,
						    c2n[3+4*cell] - 1 );
    break;
  case 5:
    for ( cell=0; cell<*ncell; cell++) gridAddPyramid( grid,
						       c2n[0+5*cell] - 1,
						       c2n[1+5*cell] - 1,
						       c2n[2+5*cell] - 1,
						       c2n[3+5*cell] - 1,
						       c2n[4+5*cell] - 1 );
    break;
  case 6:
    for ( cell=0; cell<*ncell; cell++) gridAddPrism( grid,
						     c2n[0+6*cell] - 1,
						     c2n[1+6*cell] - 1,
						     c2n[2+6*cell] - 1,
						     c2n[3+6*cell] - 1,
						     c2n[4+6*cell] - 1,
						     c2n[5+6*cell] - 1 );
    break;
  default:
    printf( "ERROR: %s: %d: Cannot handle %d node elements\n",
	    __FILE__, __LINE__, (*nodes_per_cell) );
    break;
  }
}

void FC_FUNC(gridinsertbc,GRIDINSERTBC)(int *faceId, int *nodes_per_face, int *nface, int *f2n)
{
  int face;
  int node0, node1, node2, node3;

  switch (*nodes_per_face) {
  case 3:
    for(face=0;face<*nface;face++){
      node0 = f2n[0+(*nodes_per_face)*face] - 1;
      node1 = f2n[1+(*nodes_per_face)*face] - 1;
      node2 = f2n[2+(*nodes_per_face)*face] - 1;
      if ( gridNodeGhost(grid,node0) && 
	   gridNodeGhost(grid,node1) && 
	   gridNodeGhost(grid,node2) ) {
      }else{
	gridAddFace(grid, node0, node1, node2, *faceId);
      }
    }
    break;
  case 4:
    for(face=0;face<*nface;face++){
      node0 = f2n[0+(*nodes_per_face)*face] - 1;
      node1 = f2n[1+(*nodes_per_face)*face] - 1;
      node2 = f2n[2+(*nodes_per_face)*face] - 1;
      node3 = f2n[3+(*nodes_per_face)*face] - 1;
      if ( gridNodeGhost(grid,node0) && 
	   gridNodeGhost(grid,node1) && 
	   gridNodeGhost(grid,node2) && 
	   gridNodeGhost(grid,node3) ) {
      }else{
	gridAddQuad(grid, node0, node1, node2, node3, *faceId);
      }
    }
    break;
  default:
    printf( "ERROR: %s: %d: Cannot handle %d node faces\n",
	    __FILE__, __LINE__, (*nodes_per_face) );
    break;
  }

}

void FC_FUNC(gridsetmap,GRIDSETMAP)( int *nnode, double* map )
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

void FC_FUNC(gridsetimesh,GRIDSETIMESH)( int *nnode, int *imesh )
{
  int node;
  for ( node=0; node<*nnode; node++)
    gridSetIMesh( grid, node, imesh[node] );

  if (NULL != queue) queueFree( queue );
  /* 3:xyz + 6:m + naux + 1:imesh */
  queue = queueCreate( 9 + 1 + gridNAux(grid) ); 
}

void FC_FUNC(gridsetnodepart,GRIDSETNODEPART)( int *nnode, int *part )
{
  int node;
  for ( node=0; node<*nnode; node++) gridSetNodePart(grid, node, part[node]-1);
}

void FC_FUNC(gridsetnodelocal2global,GRIDSETNODELOCAL2GLOBAL)( int *partId, int *nnodeg, 
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

void FC_FUNC(gridfreezenode,GRIDFREEZENODE)( int *nodeFortran )
{
  int nodeC;
  nodeC = (*nodeFortran)-1;
  gridFreezeNode( grid, nodeC );
}

void FC_FUNC(gridparallelloadcapri,GRIDPARALLELLOADCAPRI)( char *url, char *modeler, char *capriProject,
                             int *status )
{
  gridSetCostConstraint(grid, gridCostConstraint(grid)|gridCOST_CNST_AREAUV);
  if( grid != gridParallelGeomLoad( grid, url, modeler, capriProject ) ) {
    printf( "ERROR: %s: %d: failed to load part %s, partition %d, Modeler %s\n",
	    __FILE__,__LINE__,capriProject,gridPartId(grid),modeler );
    *status = 0;
  } else {
    *status = 1;
  }
}

void FC_FUNC(gridparallelsavecapri,GRIDPARALLELSAVECAPRI)( char *capriProject )
{
  gridParallelGeomSave( grid, capriProject );
}

void FC_FUNC(gridprojectallfaces,GRIDPROJECTALLFACES)( void )
{
  int face, node, nodes[3], faceId;
#ifdef PARALLEL_VERBOSE 
  double ar0, ar1;
#endif

#ifdef PARALLEL_VERBOSE 
  ar0 = gridMinThawedAR(grid);
#endif
  for(face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      for(node=0;node<3;node++) {
	if ( !gridNodeFrozen(grid,nodes[node]) ) {
	  if (grid != gridProjectNodeToFace(grid, nodes[node], faceId ) )
	    printf( "ERROR: %s: %d: project failed on part %d.\n",
		    __FILE__,__LINE__,gridPartId(grid) );

	}
      }
    }
  }
#ifdef PARALLEL_VERBOSE 
  ar1 = gridMinThawedAR(grid);

  printf( " %6d project faces           initial AR%14.10f final AR%14.10f\n",
	  gridPartId(grid),ar0,ar1 );
  fflush(stdout);
#endif
}


void FC_FUNC(gridtestcadparameters,GRIDTESTCADPARAMETERS)( void )
{
  int global, local; 
  int nodes[3], edge, edgeId, face, faceId;
  double oldXYZ[3], oldT, oldUV[2];
  double newXYZ[3], newT, newUV[2];
  AdjIterator it;

  for(global=0;global<gridGlobalNNode(grid);global++) {
    local = gridGlobal2Local(grid,global);
    if (EMPTY!=global) {
      for ( it = adjFirst(gridFaceAdj(grid),local); 
	    adjValid(it); 
	    it = adjNext(it) ){
	face = adjItem(it);
	gridFace(grid, face, nodes, &faceId);
	gridNodeXYZ(grid,local,oldXYZ);
	gridNodeUV(grid,local,faceId,oldUV);
	gridProjectNodeToFace(grid, local, faceId );
	gridNodeXYZ(grid,local,newXYZ);
	gridNodeUV(grid,local,faceId,newUV);
      }
      for ( it = adjFirst(gridEdgeAdj(grid),local); 
	    adjValid(it); 
	    it = adjNext(it) ){
	edge = adjItem(it);
	gridEdge(grid, edge, nodes, &edgeId);
	gridNodeXYZ(grid,local,oldXYZ);
	gridNodeT(grid,local,edgeId,&oldT);
	gridProjectNodeToEdge(grid, local, edgeId );
	gridNodeXYZ(grid,local,newXYZ);
	gridNodeT(grid,local,edgeId,&newT);
      }
    }
  }

}

void FC_FUNC(gridminar,GRIDMINAR)( double *aspectratio )
{
  *aspectratio = gridMinThawedAR( grid );
#ifdef PARALLEL_VERBOSE 
  printf( " %6d Minimum Aspect Ratio %g\n", gridPartId(grid), *aspectratio );
  fflush(stdout);
#endif
}

void FC_FUNC(gridwritetecplotsurfacezone,GRIDWRITETECPLOTSURFACEZONE)( void )
{
  char filename[256];

  sprintf(filename, "grid%04d.t", gridPartId(grid)+1 );
  gridWriteTecplotInvalid(grid,filename);
}

void FC_FUNC(gridexportfast,GRIDEXPORTFAST)( void )
{
  char filename[256];
  sprintf(filename, "grid%04d.fgrid", gridPartId(grid)+1 );
  gridExportFAST(grid,filename);
}

void FC_FUNC(gridsetcostconstraint,GRIDSETCOSTCONSTRAINT)( int *cost_constraint )
{
  gridSetCostConstraint(grid, *cost_constraint);
}

void FC_FUNC(gridconstrainsurfacenode,GRIDCONSTRAINSURFACENODE)( void )
{
  gridConstrainSurfaceNode(grid);
}

void FC_FUNC(gridparallelswap,GRIDPARALLELSWAP)( int *processor, double *ARlimit )
{
#ifdef PARALLEL_VERBOSE 
  printf(" %6d swap  processor %2d      initial AR%14.10f",
	 gridPartId(grid),*processor,gridMinAR(grid));
#endif
  if (*processor == -1) {
    int plan_size_guess, plan_chunk_size;
    int cell, nodes[4];
    double ar;
    int nodes_on_surface;
    gridParallelSwap(grid,NULL,*ARlimit);
    plan_size_guess = gridNCell(grid)/10;
    plan_chunk_size = 5000;
    plan = planCreate(plan_size_guess, plan_chunk_size);
    for (cell=0;cell<gridMaxCell(grid);cell++){
      if (grid == gridCell( grid, cell, nodes) ) {
	if ( gridCellHasGhostNode(grid,nodes)  ||
	     gridNodeNearGhost(grid, nodes[0]) ||
	     gridNodeNearGhost(grid, nodes[1]) ||
	     gridNodeNearGhost(grid, nodes[2]) ||
	     gridNodeNearGhost(grid, nodes[3]) ) {
	  /* if there are four nodes on surface it may be a two face cell */
	  nodes_on_surface = 0;
	  if ( gridGeometryFace(grid, nodes[0]) ) nodes_on_surface++;
	  if ( gridGeometryFace(grid, nodes[1]) ) nodes_on_surface++;
	  if ( gridGeometryFace(grid, nodes[2]) ) nodes_on_surface++;
	  if ( gridGeometryFace(grid, nodes[3]) ) nodes_on_surface++;
	  if ( 4 == nodes_on_surface) {
	    planAddItemWithPriority(plan,cell,0.0); /* highest priority */
	  } else {
	    /* add poor quality cells */
	    ar = gridAR(grid, nodes);
	    if ( ar < *ARlimit ) planAddItemWithPriority(plan,cell,ar);
	  }
	}
      } 
    }
    planDeriveRankingsFromPriorities(plan);
  } else {
    int ranking;
    int cell, nodes[4];
    for (ranking = 0 ; ranking < planSize(plan) ; ranking++) {
      cell = planItemWithThisRanking(plan, ranking);
      if ( grid==gridCell( grid, cell, nodes) ) {
	/* if ( grid == gridRemoveTwoFaceCell(grid, queue, cell) ) continue;*/
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[1] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[2] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[3] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[1], nodes[2] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[1], nodes[3] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[2], nodes[3] ) )
	  continue;
      }
    }
    planFree(plan); plan = NULL;
  } 

}

void FC_FUNC(gridparallelsmooth,GRIDPARALLELSMOOTH)( int *processor,
			  double *optimizationLimit, double *laplacianLimit,
                          int *geometryAllowed )
{
  GridBool localOnly, smoothOnSurface;
  localOnly = (-1 == (*processor));
  smoothOnSurface = (0 != (*geometryAllowed));
  gridParallelSmooth(grid, localOnly, *optimizationLimit, *laplacianLimit,
                     smoothOnSurface);
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

void FC_FUNC(gridparallelrelaxneg,GRIDPARALLELRELAXNEG)( int *processor, int *geometryAllowed )
{
  GridBool localOnly, smoothOnSurface;
  localOnly = (-1 == (*processor));
  smoothOnSurface = (0 != (*geometryAllowed));
  gridParallelRelaxNegativeCells(grid, localOnly, smoothOnSurface);
#ifdef PARALLEL_VERBOSE 
  if (localOnly) {
    printf( " %6d relaxN volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"local only        ",gridMinVolume(grid) );
  } else {
    printf( " %6d relaxN volume and face interior  %s    AR%14.10f\n",
	    gridPartId(grid),"near ghost only   ",gridMinVolume(grid) );
  }
  fflush(stdout);
#endif
}

void FC_FUNC(gridparallelrelaxsurf,GRIDPARALLELRELAXSURF)( int *processor )
{
  GridBool localOnly;
  localOnly = (-1 == (*processor));
  gridParallelRelaxNegativeFaceAreaUV(grid, localOnly);
}

void FC_FUNC(gridparalleladapt,GRIDPARALLELADAPT)( int *processor, 
			 double *minLength, double *maxLength )
{
#ifdef PARALLEL_VERBOSE 
  printf(" %6d adapt processor %2d ",gridPartId(grid),*processor);
#endif
  if (*processor == -1) {
    gridParallelAdapt(grid,NULL,*minLength, *maxLength);
  } else {
    gridParallelAdapt(grid,queue,*minLength, *maxLength);
  } 
}

void FC_FUNC(gridparallelpreproject,GRIDPARALLELPREPROJECT)( int *processor )
{
#ifdef PARALLEL_VERBOSE 
  printf(" %6d prepj processor %2d ",gridPartId(grid),*processor);
#endif
  if (*processor == -1) {
    gridParallelPreProject(grid,NULL);
  } else {
    gridParallelPreProject(grid,queue);
  } 
}

void FC_FUNC(queuedumpsize,QUEUEDUMPSIZE)( int *nInt, int *nDouble )
{
  queueDumpSize(queue, nInt, nDouble);
}

void FC_FUNC(queuedump,QUEUEDUMP)( int *nInt, int *nDouble, int *ints, double *doubles )
{
  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nInt);
  SUPRESS_UNUSED_COMPILER_WARNING(nDouble);

  queueDump(queue, ints, doubles);
  queueReset(queue);
}

void FC_FUNC(gridapplyqueue,GRIDAPPLYQUEUE)( int *nInt, int *nDouble, int *ints, double *doubles )
{
  Queue *appliedQueue;

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nInt);
  SUPRESS_UNUSED_COMPILER_WARNING(nDouble);

  appliedQueue = queueCreate( queueNodeSize( queue ) );
  queueLoad(appliedQueue, ints, doubles);
  gridApplyQueue(grid,appliedQueue);
  queueFree(appliedQueue);
}

void FC_FUNC(gridglobalnnode,GRIDGLOBALNNODE)( int *nnodeg )
{
  *nnodeg = gridGlobalNNode(grid);
}

void FC_FUNC(gridglobalshift,GRIDGLOBALSHIFT)( int *oldnnodeg, int *newnnodeg, int *nodeoffset )
{
  gridGlobalShiftNode( grid, *oldnnodeg, *newnnodeg, *nodeoffset);
  queueGlobalShiftNode( queue, *oldnnodeg, *nodeoffset);
}

void FC_FUNC(gridrenumberglobalnodes,GRIDRENUMBERGLOBALNODES)( int *nnode, int *new2old )
{
  int i;
  for (i=0;i<*nnode;i++) new2old[i]--;
  gridRenumberGlobalNodes( grid, *nnode, new2old );
  for (i=0;i<*nnode;i++) new2old[i]++;
}

void FC_FUNC(gridnunusednodeglobal,GRIDNUNUSEDNODEGLOBAL)( int *nunused )
{
  *nunused = gridNUnusedNodeGlobal( grid );
}

void FC_FUNC(gridgetunusednodeglobal,GRIDGETUNUSEDNODEGLOBAL)( int *nunused, int *unused )
{
  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nunused);

  gridGetUnusedNodeGlobal( grid, unused );
}

void FC_FUNC(gridjoinunusednodeglobal,GRIDJOINUNUSEDNODEGLOBAL)( int *nunused, int *unused )
{
  int i;
  for (i=0;i<(*nunused);i++) gridJoinUnusedNodeGlobal( grid, unused[i] );
}

void FC_FUNC(gridcopyunusednodeglobal,GRIDCOPYUNUSEDNODEGLOBAL)( int *nunused, int *unused )
{
  int i;

  if (NULL != grid->unusedNodeGlobal)
    free(grid->unusedNodeGlobal);

  grid->nUnusedNodeGlobal = grid->maxUnusedNodeGlobal = *nunused;
  grid->unusedNodeGlobal =
    (int *)malloc(grid->maxUnusedNodeGlobal * sizeof(int));

  for (i=0;i<(*nunused);i++)
    grid->unusedNodeGlobal[i] = unused[i];
}

void FC_FUNC(grideliminateunusednodeglobal,GRIDELIMINATEUNUSEDNODEGLOBAL)(  )
{
  gridEliminateUnusedNodeGlobal( grid );
}

void FC_FUNC(gridsortfun3d,GRIDSORTFUN3D)( int *nnodes0, int *nnodes01, int *nnodesg )
{
  gridSortNodeFUN3D( grid, nnodes0 );
  *nnodes01 = gridNNode(grid);
  *nnodesg = gridGlobalNNode(grid);
}

void FC_FUNC(gridgetnodes,GRIDGETNODES)( int *nnode, int *l2g, double *x, double *y, double *z)
{
  int node;
  double xyz[3];

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nnode);

  for (node=0;node<gridNNode(grid);node++) {
    l2g[node] = gridNodeGlobal(grid,node)+1;
    gridNodeXYZ(grid,node,xyz);
    x[node] = xyz[0];
    y[node] = xyz[1];
    z[node] = xyz[2];
  }
}

void FC_FUNC(gridgetmap,GRIDGETMAP)( int *nnode, double *map)
{
  int node;

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nnode);

  for (node=0;node<gridNNode(grid);node++) {
    gridMap(grid,node,&map[6*node]);
  }
}

void FC_FUNC(gridgetfreezestate,GRIDGETFREEZESTATE)( int *nnode, int *state)
{
  int node;

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nnode);

  for (node=0;node<gridNNode(grid);node++) 
    {
      if ( gridNodeFrozen( grid, node ) )
	{
	  state[node] = 1;
	} 
      else 
	{
	  state[node] = 0;
	}
    }
}

void FC_FUNC(gridgetimesh,GRIDGETIMESH)( int *nnode, int *imesh)
{
  int node;

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(nnode);

  for (node=0;node<gridNNode(grid);node++) {
    imesh[node] = gridIMesh(grid,node);
  }
}

void FC_FUNC(gridgetncell,GRIDGETNCELL)( int *nodes_per_cell, int *ncell )
{

  switch (*nodes_per_cell) {
  case 4:
    *ncell = gridNCell(grid);
    break;
  case 5:
    *ncell = gridNPyramid(grid);
    break;
  case 6:
    *ncell = gridNPrism(grid);
    break;
  default:
    *ncell = 0;
    printf( "ERROR: %s: %d: Cannot handle %d node elements\n",
	    __FILE__, __LINE__, (*nodes_per_cell) );
    break;
  }

}

void FC_FUNC(gridgetcell,GRIDGETCELL)( int *nodes_per_cell, int *ncell, int *c2n )
{
  int cell, total;
  int node, nodes[6];

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(ncell);

  switch (*nodes_per_cell) {
  case 4:
    total = 0;
    for ( cell = 0 ; cell < gridMaxCell(grid) ; cell++ )
      {
	if ( grid == gridCell(grid,cell,nodes) )
	  {
	    for ( node = 0 ; node < (*nodes_per_cell) ; node++ )
	      c2n[node+(*nodes_per_cell)*total] = nodes[node] + 1;
	    total++;
	  }
      }
    break;
  case 5:
    total = 0;
    for ( cell = 0 ; cell < gridNPyramid(grid) ; cell++ )
      {
	if ( grid == gridPyramid(grid,cell,nodes) )
	  {
	    for ( node = 0 ; node < (*nodes_per_cell) ; node++ )
	      c2n[node+(*nodes_per_cell)*total] = nodes[node] + 1;
	    total++;
	  }
      }
    break;
  case 6:
    total = 0;
    for ( cell = 0 ; cell < gridNPrism(grid) ; cell++ )
      {
	if ( grid == gridPrism(grid,cell,nodes) )
	  {
	    for ( node = 0 ; node < (*nodes_per_cell) ; node++ )
	      c2n[node+(*nodes_per_cell)*total] = nodes[node] + 1;
	    total++;
	  }
      }
    break;
  default:
    printf( "ERROR: %s: %d: Cannot handle %d node elements\n",
	    __FILE__, __LINE__, (*nodes_per_cell) );
    break;
  }

}

void FC_FUNC(gridgetbcsize,GRIDGETBCSIZE)( int *ibound, int *nodes_per_face, int *nface )
{
  int face, nodes[4], id;
  
  switch (*nodes_per_face) {
  case 3:
    *nface = 0;
    for (face=0;face<gridMaxFace(grid);face++) {
      if ( grid == gridFace(grid,face,nodes,&id) ) {
	if ( *ibound == id ) (*nface)++;
      }
    }
    break;
  case 4:
    *nface = 0;
    for (face=0;face<gridNQuad(grid);face++) {
      if ( grid == gridQuad(grid,face,nodes,&id) ) {
	if ( *ibound == id ) (*nface)++;
      }
    }
    break;
  default:
    *nface = 0;
    printf( "ERROR: %s: %d: Cannot handle %d node faces\n",
	    __FILE__, __LINE__, (*nodes_per_face) );
    break;
  }
}

void FC_FUNC(gridgetbc,GRIDGETBC)( int *ibound, int *nodes_per_face, int *nface, int *f2n )
{
  int face, n, nodes[4], id;

  switch (*nodes_per_face) {
  case 3:
    n = 0;
    for (face=0;face<gridMaxFace(grid);face++) {
      if ( grid == gridFace(grid,face,nodes,&id) ) {
	if ( *ibound == id ) {
	  f2n[0+(*nodes_per_face)*n] = nodes[0]+1;
	  f2n[1+(*nodes_per_face)*n] = nodes[1]+1;
	  f2n[2+(*nodes_per_face)*n] = nodes[2]+1;
	  n++;
	}
      }
    }
    break;
  case 4:
    n = 0;
    for (face=0;face<gridNQuad(grid);face++) {
      if ( grid == gridQuad(grid,face,nodes,&id) ) {
	if ( *ibound == id ) {
	  f2n[0+(*nodes_per_face)*n] = nodes[0]+1;
	  f2n[1+(*nodes_per_face)*n] = nodes[1]+1;
	  f2n[2+(*nodes_per_face)*n] = nodes[2]+1;
	  f2n[3+(*nodes_per_face)*n] = nodes[3]+1;
	  n++;
	}
      }
    }
    break;
  default:
    *nface = 0;
    printf( "ERROR: %s: %d: Cannot handle %d node faces\n",
	    __FILE__, __LINE__, (*nodes_per_face) );
    break;
  }

}

void FC_FUNC(gridsetnaux,GRIDSETNAUX)( int *naux )
{
  int imesh_index;

  gridSetNAux(grid, *naux);

  imesh_index = (gridHaveIMesh(grid)?1:0);
  if (NULL != queue) queueFree( queue );
  /* 3:xyz + 6:m + naux + 1:imesh */
  queue = queueCreate( 9 + imesh_index + gridNAux(grid) ); 
}

void FC_FUNC(gridgetnaux,GRIDGETNAUX)( int *naux )
{
  *naux = gridNAux(grid);
}

void FC_FUNC(gridsetauxvector,GRIDSETAUXVECTOR)( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    gridSetAux(grid,node,(*offset),x[node]);
  }
}

void FC_FUNC(gridsetauxmatrix,GRIDSETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

void FC_FUNC(gridsetauxmatrix3,GRIDSETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      gridSetAux(grid,node,(*offset)+dim,x[dim+(*ndim)*node]);
    }
  }
}

void FC_FUNC(gridgetauxvector,GRIDGETAUXVECTOR)( int *nnode, int *offset, double *x )
{
  int node;
  for (node=0;node<(*nnode);node++) {
    x[node] = gridAux(grid,node,(*offset));
  }
}

void FC_FUNC(gridgetauxmatrix,GRIDGETAUXMATRIX)( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

void FC_FUNC(gridgetauxmatrix3,GRIDGETAUXMATRIX3)( int *ndim, int *nnode, int *offset, double *x )
{
  int node, dim;
  for (node=0;node<(*nnode);node++) {
    for (dim=0;dim<(*ndim);dim++){
      x[dim+(*ndim)*node] = gridAux(grid,node,(*offset)+dim);
    }
  }
}

void FC_FUNC(gridghostcount,GRIDGHOSTCOUNT)( int *nproc, int *count )
{
  if (grid!=gridGhostDataCountByPartition(grid, (*nproc), count))
    printf("%s: %d: gridghostcount_ %s returned error\n", __FILE__, __LINE__,
	   "gridGhostDataCountByPartition");
}

void FC_FUNC(gridloadghostnodes,GRIDLOADGHOSTNODES)( int *nproc, int *clientindex,
			  int *clientsize, int *localnode, int *globalnode )
{
  int node, part;
  int *count;
  int face, faceids, faceid[MAXFACEIDDEG];
  int edge, edgeids, edgeid[MAXEDGEIDDEG];

  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(clientsize);

  count = malloc( (*nproc) * sizeof(int) );

  for(node=0;node<(*nproc);node++) count[node] = 0;
  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) {
      part = gridNodePart(grid,node);
      localnode[ count[part]+clientindex[part]-1] = node+1;
      globalnode[count[part]+clientindex[part]-1] = gridNodeGlobal(grid,node)+1;
      count[part]++;
      gridNodeFaceId(grid, node, MAXFACEIDDEG, &faceids, faceid );
      if (faceids>0) {
	localnode[ count[part]+clientindex[part]-1] = -faceids;
	globalnode[count[part]+clientindex[part]-1] = -faceids;
	count[part]++;
	for (face=0;face<faceids;face++) {
	  localnode[ count[part]+clientindex[part]-1] = -faceid[face];
	  globalnode[count[part]+clientindex[part]-1] = -faceid[face];
	  count[part]++;
	}
      }
      gridNodeEdgeId(grid, node, MAXEDGEIDDEG, &edgeids, edgeid );
      if (edgeids>0) {
	localnode[ count[part]+clientindex[part]-1] = -edgeids;
	globalnode[count[part]+clientindex[part]-1] = -edgeids;
	count[part]++;
	for (edge=0;edge<edgeids;edge++) {
	  localnode[ count[part]+clientindex[part]-1] = -edgeid[edge];
	  globalnode[count[part]+clientindex[part]-1] = -edgeid[edge];
	  count[part]++;
	}
      }
    }
  }
  free(count);
}

void FC_FUNC(gridloadlocalnodes,GRIDLOADLOCALNODES)( int *nnode, int *global, int *local )
{
  int node, globalnode, localnode;

  localnode=0;
  node=0;
  while (node<(*nnode)) {
    if (global[node] > 0) {
      globalnode = global[node]-1;
      localnode = gridGlobal2Local(grid, globalnode);
      if (!gridValidNode(grid,localnode)) 
	printf("%d: ERROR: %s: %d: invalid node local %d global %d.\n",
	       gridPartId(grid),__FILE__, __LINE__, localnode, globalnode);
      local[node] = 1+localnode;
      node++;
    } else {
      local[node] = global[node];
      node++;
    } 
  }
}

void FC_FUNC(gridloadglobalnodedata,GRIDLOADGLOBALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data )
{
  int node, localnode;
  int face, faceids, faceId;
  int edge, edgeids, edgeId;
  double t, uv[2], xyz[3];

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
      if ( (*ndim) > 3 ) {
	data[3+(*ndim)*node] = 0.0;
	if (gridNodeFrozen( grid, localnode )) data[3+(*ndim)*node] = 1.0;
      }
      node++;
    } else {
      faceids = -nodes[node];
      node++;
      for(face=0;face<faceids;face++) {
	faceId = -nodes[node];
	gridNodeUV(grid,localnode,faceId,uv);
	data[0+(*ndim)*node] = uv[0];
	data[1+(*ndim)*node] = uv[1];
	node++;
      }
      if (node<(*nnode) && nodes[node] < 0) {
	edgeids = -nodes[node];
	node++;
	for(edge=0;edge<edgeids;edge++) {
	  edgeId = -nodes[node];
	  gridNodeT(grid,localnode,edgeId,&t);
	  data[0+(*ndim)*node] = t;
	  node++;
	}
      }
    } 
  }
}

void FC_FUNC(gridsetlocalnodedata,GRIDSETLOCALNODEDATA)( int *ndim, int *nnode, int *nodes, double *data )
{
  int node, localnode;
  int face, faceids, faceId;
  int edge, edgeids, edgeId;

  localnode=0;
  node=0;
  while (node<(*nnode)) {
    if (nodes[node] > 0) {
      localnode = nodes[node]-1;
      if (grid != gridSetNodeXYZ(grid,localnode,&data[(*ndim)*node]) )
	printf("ERROR: %s: %d: set invalid node %d .\n",
	       __FILE__, __LINE__, nodes[node]-1);
      if ( (*ndim) > 3 ) {
	gridThawNode( grid, localnode );
	if ( data[3+(*ndim)*node] > 0.5 ) gridFreezeNode( grid, localnode );
      }
      node++;
    } else {
      faceids = -nodes[node];
      node++;
      for(face=0;face<faceids;face++) {
	faceId = -nodes[node];
	gridSetNodeUV(grid,localnode,faceId,
		      data[0+(*ndim)*node],data[1+(*ndim)*node]);
	node++;
      }
      if (node<(*nnode) && nodes[node] < 0) {
	edgeids = -nodes[node];
	node++;
	for(edge=0;edge<edgeids;edge++) {
	  edgeId = -nodes[node];
	  gridSetNodeT(grid,localnode,edgeId,data[0+(*ndim)*node]);
	  node++;
	}
      }
    } 
  }

#ifdef PARALLEL_VERBOSE 
  printf( " %6d update xfer                      %s    AR%14.10f\n",
	  gridPartId(grid),"                  ",gridMinAR(grid) );
  fflush(stdout);
#endif
}

void FC_FUNC(gridmovesetprojectiondisp,GRIDMOVESETPROJECTIONDISP)( void )
{
  gm = gridmoveCreate( grid );
  gridmoveProjectionDisplacements( gm );
}

void FC_FUNC(gridmoverelaxstartup,GRIDMOVERELAXSTARTUP)( int *relaxationScheme )
{
  gridmoveRelaxationStartUp(gm, *relaxationScheme);
}

void FC_FUNC(gridmoverelaxstartstep,GRIDMOVERELAXSTARTSTEP)( double *position)
{
  gridmoveRelaxationStartStep( gm, *position );
}

void FC_FUNC(gridmoverelaxsubiter,GRIDMOVERELAXSUBITER)( double *residual)
{
  gridmoveRelaxationSubIteration( gm, residual );
}

void FC_FUNC(gridmoverelaxshutdown,GRIDMOVERELAXSHUTDOWN)( void )
{
  gridmoveRelaxationShutDown(gm);
}

void FC_FUNC(gridmoveapplydisplacements,GRIDMOVEAPPLYDISPLACEMENTS)( void )
{
  gridmoveApplyDisplacements(gm);
}

void FC_FUNC(gridmovedataleadingdim,GRIDMOVEDATALEADINGDIM)( int *ndim )
{
  gridmoveDataLeadingDimension( gm, ndim );
}

void FC_FUNC(gridmoveinitializempitest,GRIDMOVEINITIALIZEMPITEST)( void )
{
  gridmoveInitializeMPITest(gm);
}

void FC_FUNC(gridmovecompletempitest,GRIDMOVECOMPLETEMPITEST)( void )
{
  gridmoveCompleteMPITest(gm);
}

void FC_FUNC(gridmoveloadlocalnodedata,GRIDMOVELOADLOCALNODEDATA)( int *ndim, int *nnode, 
				 int *nodes, double *data )
{
  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(ndim);

  gridmoveLoadFortranNodeData( gm, *nnode, nodes, data);
}

void FC_FUNC(gridmovesetlocalnodedata,GRIDMOVESETLOCALNODEDATA)( int *ndim, int *nnode, 
				 int *nodes, double *data )
{
  /* this is for the fortran interface */
  SUPRESS_UNUSED_COMPILER_WARNING(ndim);

  gridmoveSetFortranNodeData( gm, *nnode, nodes, data);
}

void FC_FUNC(gridmovefree,GRIDMOVEFREE)( void )
{
  gridmoveFree( gm ); gm = NULL;
}

void FC_FUNC(gridgeomsize,GRIDGEOMSIZE)( int *nGeomNode, int *nGeomEdge, int *nGeomFace )
{
  *nGeomNode = gridNGeomNode(grid);
  *nGeomEdge = gridNGeomEdge(grid);
  *nGeomFace = gridNGeomFace(grid);
}

void FC_FUNC(gridlocalboundnode,GRIDLOCALBOUNDNODE)( int *nBoundNode )
{
  int node, nnode;
  nnode = 0;
  for(node=0;node<gridMaxNode(grid);node++) 
    if ( gridGeometryFace(grid,node) && gridNodeLocal(grid,node) ) nnode++;

  *nBoundNode = nnode;
}

void FC_FUNC(gridgeomedgeendpoints,GRIDGEOMEDGEENDPOINTS)( int *edgeId, int *endPoints )
{
  endPoints[0] = 1 + gridGeomEdgeStart(grid,*edgeId);
  endPoints[1] = 1 + gridGeomEdgeEnd(grid,*edgeId);
}

void FC_FUNC(gridmaxedge,GRIDMAXEDGE)( int *maxedge )
{
  *maxedge = gridMaxEdge( grid );
}

void FC_FUNC(gridedge,GRIDEDGE)( int *edge, int *edgeId,
		int *globalnodes, int *nodeparts,
		double *t, double *xyz)
{
  int localnodes[2];
  if (grid==gridEdge(grid, (*edge)-1, localnodes, edgeId)) {
    globalnodes[0] = 1 + gridNodeGlobal(grid,localnodes[0]);
    globalnodes[1] = 1 + gridNodeGlobal(grid,localnodes[1]);
    nodeparts[0] = gridNodePart(grid,localnodes[0]);
    nodeparts[1] = gridNodePart(grid,localnodes[1]);
    gridNodeT(grid,localnodes[0],*edgeId,&t[0]);
    gridNodeT(grid,localnodes[1],*edgeId,&t[1]);
    gridNodeXYZ(grid,localnodes[0],&xyz[0*3]);
    gridNodeXYZ(grid,localnodes[1],&xyz[1*3]);
  }else{
    *edgeId = EMPTY;
  }
}

void FC_FUNC(gridupdateedgegrid,GRIDUPDATEEDGEGRID)(int *edgeId, int *nCurveNode, double *xyz, double *t)
{
  gridUpdateEdgeGrid( grid, *edgeId, *nCurveNode, xyz, t);
}

void FC_FUNC(gridmaxface,GRIDMAXFACE)( int *maxface )
{
  *maxface = gridMaxFace( grid );
}

void FC_FUNC(gridface,GRIDFACE)( int *face, int *faceId, 
		int *globalnodes, int *nodeparts,
		double *uv, double *xyz)
{
  int localnodes[3];
  if (grid==gridFace(grid, (*face)-1, localnodes, faceId)) {
    globalnodes[0] = 1 + gridNodeGlobal(grid,localnodes[0]);
    globalnodes[1] = 1 + gridNodeGlobal(grid,localnodes[1]);
    globalnodes[2] = 1 + gridNodeGlobal(grid,localnodes[2]);
    nodeparts[0] = gridNodePart(grid,localnodes[0]);
    nodeparts[1] = gridNodePart(grid,localnodes[1]);
    nodeparts[2] = gridNodePart(grid,localnodes[2]);
    gridNodeUV(grid,localnodes[0],*faceId,&uv[0*2]);
    gridNodeUV(grid,localnodes[1],*faceId,&uv[1*2]);
    gridNodeUV(grid,localnodes[2],*faceId,&uv[2*2]);
    gridNodeXYZ(grid,localnodes[0],&xyz[0*3]);
    gridNodeXYZ(grid,localnodes[1],&xyz[1*3]);
    gridNodeXYZ(grid,localnodes[2],&xyz[2*3]);
  }else{
    *faceId = EMPTY;
  }
}

void FC_FUNC(gridfaceedgecount,GRIDFACEEDGECOUNT)( int *faceId, int *faceEdgeCount )
{
  *faceEdgeCount = gridFaceEdgeCount( *faceId );
}

void FC_FUNC(gridfaceedgel2g,GRIDFACEEDGEL2G)( int *faceId, int *faceEdgeCount, int *local2global )
{
  int i;
  gridFaceEdgeLocal2Global( grid, *faceId, *faceEdgeCount, local2global );
  for (i=0;i<*faceEdgeCount;i++) local2global[i]++;
}

void FC_FUNC(gridupdategeometryface,GRIDUPDATEGEOMETRYFACE)( int *faceId, int *nnode, double *xyz, double *uv,
			      int *nface, int *f2n )
{
  int i;
  for( i=0 ; i<3*(*nface) ; i++) f2n[i]--;
  gridUpdateGeometryFace( grid, *faceId, *nnode, xyz, uv,
			  *nface, f2n );
  for( i=0 ; i<3*(*nface) ; i++) f2n[i]++; 
}

void FC_FUNC(gridcreateshellfromfaces,GRIDCREATESHELLFROMFACES)( void )
{
  gridCreateShellFromFaces( grid );
}

