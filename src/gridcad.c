
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
#include "CADGeom/CADGeom.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "adj.h"
#include "gridStruct.h"

Grid *gridForceNodeToEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  t = DBL_MAX;

  if (!CADGeom_NearestOnEdge( vol, edgeId, xyz, &t, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;

  return grid;
}

Grid *gridForceNodeToFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  uv[0] = DBL_MAX;
  uv[1] = DBL_MAX;

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;

  return grid;
}

Grid *gridProjectNodeToEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL;

  if (!CADGeom_NearestOnEdge( vol, edgeId, xyz, &t, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;
  if ( grid != gridSetNodeT( grid, node, edgeId, t ) ) return NULL;

  return grid;
}

Grid *gridProjectNodeToFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;
  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  return grid;
}

Grid *gridSafeProjectNode(Grid *grid, int node, double ratio )
{
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(grid->edgeAdj, node));
    edgeId = grid->edgeId[edge];
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId, ratio ) ) 
      return NULL;
    for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
      face = adjItem(it);
      faceId = grid->faceId[face];
      if ( grid != gridSafeProjectNodeToFace( grid, node, faceId, 1.0 ) ) 
	return NULL;
    }
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId, 1.0 ) ) return NULL;
    return grid;
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(grid->faceAdj, node));
    faceId = grid->faceId[face];
    if ( grid != gridSafeProjectNodeToFace( grid, node, faceId, ratio ) ) 
      return NULL;
    return grid;
  }

  return grid;
}

Grid *gridSafeProjectNodeToEdge(Grid *grid, int node, int edgeId, double ratio )
{
  double origxyz[3], origt;
  double projxyz[3];
  int attempts;

  origxyz[0] = grid->xyz[0+3*node];
  origxyz[1] = grid->xyz[1+3*node];
  origxyz[2] = grid->xyz[2+3*node];
  if ( grid != gridNodeT( grid, node, edgeId, &origt ) ) return NULL; 

  if ( grid != gridProjectNodeToEdge( grid, node, edgeId ) ) return NULL;

  if (!gridNegCellAroundNode( grid, node )) return grid;

  gridSetNodeT( grid, node, edgeId, origt );

  projxyz[0] = grid->xyz[0+3*node];
  projxyz[1] = grid->xyz[1+3*node];
  projxyz[2] = grid->xyz[2+3*node];

  attempts = 0;
  while ( attempts < 10 && gridNegCellAroundNode( grid, node ) ) {
    attempts++;
    grid->xyz[0+3*node] = (1.0-ratio)*grid->xyz[0+3*node] + ratio*origxyz[0];
    grid->xyz[1+3*node] = (1.0-ratio)*grid->xyz[1+3*node] + ratio*origxyz[1];
    grid->xyz[2+3*node] = (1.0-ratio)*grid->xyz[2+3*node] + ratio*origxyz[2];
  }

  if (gridNegCellAroundNode( grid, node )) {
    grid->xyz[0+3*node] = origxyz[0];
    grid->xyz[1+3*node] = origxyz[1];
    grid->xyz[2+3*node] = origxyz[2];
  }

  return NULL;
}

Grid *gridSafeProjectNodeToFace(Grid *grid, int node, int faceId, double ratio )
{
  double origxyz[3], origuv[2];
  double projxyz[3];
  int attempts;

  origxyz[0] = grid->xyz[0+3*node];
  origxyz[1] = grid->xyz[1+3*node];
  origxyz[2] = grid->xyz[2+3*node];
  if ( grid != gridNodeUV( grid, node, faceId, origuv ) ) return NULL; 

  if ( grid != gridProjectNodeToFace( grid, node, faceId ) ) return NULL;

  if (!gridNegCellAroundNode( grid, node )) return grid;

  gridSetNodeUV( grid, node, faceId, origuv[0], origuv[1] );

  projxyz[0] = grid->xyz[0+3*node];
  projxyz[1] = grid->xyz[1+3*node];
  projxyz[2] = grid->xyz[2+3*node];

  attempts = 0;
  while ( attempts < 10 && gridNegCellAroundNode( grid, node ) ) {
    attempts++;
    grid->xyz[0+3*node] = (1.0-ratio)*grid->xyz[0+3*node] + ratio*origxyz[0];
    grid->xyz[1+3*node] = (1.0-ratio)*grid->xyz[1+3*node] + ratio*origxyz[1];
    grid->xyz[2+3*node] = (1.0-ratio)*grid->xyz[2+3*node] + ratio*origxyz[2];
  }

  if (gridNegCellAroundNode( grid, node )) {
    grid->xyz[0+3*node] = origxyz[0];
    grid->xyz[1+3*node] = origxyz[1];
    grid->xyz[2+3*node] = origxyz[2];
  }
 

  return NULL;
}

Grid *gridProject(Grid *grid)
{
  int node;
  int notProjected;
  notProjected = 0;

  for (node=0;node<grid->maxnode;node++)
    if ( gridValidNode( grid, node ) )
      if ( gridSafeProjectNode(grid,node,1.0) != grid ) notProjected++;

  if (notProjected > 0){
    printf("gridProject: %d of %d nodes not projected.\n",
	   notProjected,grid->nnode);
    return NULL;
  }

  return grid;
}
Grid *gridRobustProject(Grid *grid)
{
  int  node;
  int notProjected;
  notProjected = 0;
  for (node=0;node<grid->maxnode;node++)
    if ( gridValidNode( grid, node ) && !gridNodeFrozen(grid, node) ) 
      if (gridRobustProjectNode( grid, node)!= grid ) 
	notProjected++;

  if (notProjected > 0){
    printf("gridRobustProject: %d of %d nodes not projected.\n",
	   notProjected,grid->nnode);
    return NULL;
  }

  return grid;
}

Grid *gridRobustProjectNode(Grid *grid, int node)
{
  int i, j, goodnode, nodes[4], level2nodes[4];
  AdjIterator it, level2;
  double xyz[3];

  if ( !gridValidNode( grid, node ) ) return NULL;
  
  if ( gridSafeProjectNode(grid,node,0.95) != grid ) {
    for ( it = adjFirst(grid->cellAdj,node); 
	  adjValid(it); 
	  it = adjNext(it) ){
      gridCell(grid, adjItem(it), nodes);
      for (i=0;i<4;i++)
	if (!gridGeometryFace( grid, nodes[i])) 
	  gridSmoothNode( grid, nodes[i]);
    }      
    gridSwapNearNodeExceptBoundary( grid, node);
    for ( it = adjFirst(grid->cellAdj,node); 
	  adjValid(it); 
	  it = adjNext(it) ){
      gridCell(grid, adjItem(it), nodes);
      for (i=0;i<4;i++)
	if (!gridGeometryFace( grid, nodes[i])) 
	  gridSmoothNode( grid, nodes[i]);
    }      
    if ( gridSafeProjectNode(grid,node,0.9) != grid ) {
      for ( it = adjFirst(grid->cellAdj,node); 
	    adjValid(it); 
	    it = adjNext(it) ){
	gridCell(grid, adjItem(it), nodes);
	for (i=0;i<4;i++) {
	  for ( level2 = adjFirst(grid->cellAdj,nodes[i]); 
		adjValid(level2); 
		level2 = adjNext(level2) ){
	    gridCell(grid, adjItem(level2), level2nodes);
	    for (j=0;j<4;j++) {
	      if (!gridGeometryFace( grid, level2nodes[j])) 
		gridSmoothNode( grid, level2nodes[j]);
	    }
	  }
	}
      }
      if ( gridSafeProjectNode(grid,node,1.0) != grid ){
	gridNodeXYZ(grid,node,xyz);
	printf(" try to collapse-project %d X %10.5f Y %10.5f Z %10.5f\n",
	       node,xyz[0],xyz[1],xyz[2]);
	for ( it = adjFirst(grid->cellAdj,node); 
	      adjValid(it); 
	      it = adjNext(it) ){
	  gridCell(grid, adjItem(it), nodes);
	  for (i=0;i<4;i++) {
	    goodnode = nodes[i];
	    if ( node != goodnode && 
		 gridGeometryFace( grid, goodnode) &&
		 grid == gridSafeProjectNode( grid, goodnode, 1.0 ) ) { 
	      if ( grid == gridCollapseEdge(grid, goodnode, node, 0.0 ) ){
		printf(" got it ! %d\n",goodnode);
		return grid;
	      }
	    }
	  }
	}      
	return NULL;
      }
    }
  }

  return grid;
}

Grid *gridSmoothNode(Grid *grid, int node )
{
  double xyzProj[3], uv[2], t;
  double ar, dARdx[3];
  double mr, dMRdx[3];
  double du[3], dv[3], dt[3];
  double dARdu[2], dARdt;
  int vol =1;
  int face, faceId;
  int edge, edgeId;
  int maxsmooth;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(grid->edgeAdj, node));
    edgeId = grid->edgeId[edge];
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    gridNodeFaceMRDerivative ( grid, node, &mr, dMRdx);
    gridNodeT( grid, node, edgeId, &t);
    if ( !CADGeom_PointOnEdge( vol, edgeId,   
			       t, xyzProj, 1, dt, NULL) )
      printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
    if (ar<mr) {
      dARdt = dARdx[0]*dt[0] + dARdx[1]*dt[1] + dARdx[2]*dt[2];
    }else{
      dARdt = dMRdx[0]*dt[0] + dMRdx[1]*dt[1] + dMRdx[2]*dt[2];
    }
    return gridOptimizeT( grid, node, dARdt );
  }
  if ( gridGeometryFace( grid, node ) ) {
    for (maxsmooth=0;maxsmooth<3;maxsmooth++) {
      face = adjItem(adjFirst(grid->faceAdj, node));
      faceId = grid->faceId[face];
      gridNodeARDerivative ( grid, node, &ar, dARdx);
      gridNodeFaceMRDerivative ( grid, node, &mr, dMRdx);
      gridNodeUV( grid, node, faceId, uv);
      if ( !CADGeom_PointOnFace( vol, faceId,   
				 uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
	printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
      
      if (ar<mr) {
	dARdu[0] = dARdx[0]*du[0] + dARdx[1]*du[1] + dARdx[2]*du[2] ; 
	dARdu[1] = dARdx[0]*dv[0] + dARdx[1]*dv[1] + dARdx[2]*dv[2] ; 
      }else{
	dARdu[0] = dMRdx[0]*du[0] + dMRdx[1]*du[1] + dMRdx[2]*du[2] ; 
	dARdu[1] = dMRdx[0]*dv[0] + dMRdx[1]*dv[1] + dMRdx[2]*dv[2] ; 
      }
      if (grid != gridOptimizeUV( grid, node, dARdu ) ) return NULL;
    }
    return grid;
  }
  if (FALSE) {
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    return gridOptimizeXYZ( grid, node, dARdx );
  }else{
    maxsmooth = 40;
    while ( grid == gridSmoothNodeQP(grid,node) && maxsmooth >0) maxsmooth--;
    return grid;
  }
}


Grid *gridSmoothNodeFaceMR(Grid *grid, int node )
{
  double xyzProj[3], uv[2], t;
  double mr, dMRdx[3];
  double du[3], dv[3], dt[3];
  double dMRdu[2], dMRdt;
  int vol =1;
  int face, faceId;
  int edge, edgeId;

  if ( gridGeometryEdge( grid, node ) ) return grid;
  if ( !gridGeometryFace( grid, node ) ) return grid;


  face = adjItem(adjFirst(grid->faceAdj, node));
  faceId = grid->faceId[face];
  gridNodeFaceMRDerivative ( grid, node, &mr, dMRdx);
  gridNodeUV( grid, node, faceId, uv);
  if ( !CADGeom_PointOnFace( vol, faceId,   
			     uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );

  dMRdu[0] = dMRdx[0]*du[0] + dMRdx[1]*du[1] + dMRdx[2]*du[2] ; 
  dMRdu[1] = dMRdx[0]*dv[0] + dMRdx[1]*dv[1] + dMRdx[2]*dv[2] ; 

  return gridOptimizeFaceUV( grid, node, dMRdu );
}

Grid *gridOptimizeT(Grid *grid, int node, double dt )
{
  int vol =1;
  double tOrig, t;
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;
  double gold;
  double alpha[2], ar[2], mr;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  edge = adjItem(adjFirst(grid->edgeAdj, node));
  edgeId = grid->edgeId[edge];

  gridNodeT( grid, node, edgeId, &tOrig);

  alpha[0] = 0.0;
  t = tOrig + alpha[0]*dt;
  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[0] );
  gridNodeFaceMR( grid, node, &mr );
  ar[0] = MIN(mr, ar[0]);

  alpha[1] = 1.0e-10;
  t = tOrig + alpha[1]*dt;
  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[1] );
  gridNodeFaceMR( grid, node, &mr );
  ar[1] = MIN(mr, ar[1]);

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    t = tOrig + alpha[1]*dt;
    if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			       0, NULL, NULL) )
      printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
    gridNodeAR( grid, node, &ar[1] );
    gridNodeFaceMR( grid, node, &mr );
    ar[1] = MIN(mr, ar[1]);
  }

  t = tOrig + alpha[0]*dt;
  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  gridSetNodeT(grid, node, edgeId, t);

  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    faceId = grid->faceId[face];
    if ( grid != gridSafeProjectNodeToFace( grid, node, faceId, 1.0 ) ) 
      return NULL;
  }

  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  
  return grid;
}

Grid *gridOptimizeUV(Grid *grid, int node, double *dudv )
{
  int vol =1;
  double uvOrig[2], uv[2];
  int face, faceId;
  double gold;
  double alpha[2], ar[2], mr;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(grid->faceAdj, node));
  faceId = grid->faceId[face];

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[0] );
  gridNodeFaceMR( grid, node, &mr );
  ar[0] = MIN(mr, ar[0]);

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[1] );
  gridNodeFaceMR( grid, node, &mr );
  ar[1] = MIN(mr, ar[1]);

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];

    if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			       0, NULL, NULL, NULL, NULL, NULL) )
      printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
    gridNodeAR( grid, node, &ar[1] );
  gridNodeFaceMR( grid, node, &mr );
  ar[1] = MIN(mr, ar[1]);
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridSetNodeUV(grid, node, faceId, uv[0], uv[1]);

//printf("node %d alpha %e ar %f uv %e %e\n",node,alpha[0],ar[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeFaceUV(Grid *grid, int node, double *dudv )
{
  int vol =1;
  double uvOrig[2], uv[2];
  int face, faceId;
  double gold;
  double alpha[2], mr[2], ar;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(grid->faceAdj, node));
  faceId = grid->faceId[face];

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridNodeFaceMR( grid, node, &mr[0] );

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridNodeFaceMR( grid, node, &mr[1] );

  iter = 0;
  gridNodeAR( grid, node, &ar ); 
  while ( mr[1] > mr[0] && mr[1] > 0.0 && iter < 100 && ar > 0.1){
    iter++;
    alpha[0] = alpha[1]; mr[0] = mr[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];

    if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			       0, NULL, NULL, NULL, NULL, NULL) )
      printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
    gridNodeFaceMR( grid, node, &mr[1] );
    gridNodeAR( grid, node, &ar );
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridSetNodeUV(grid, node, faceId, uv[0], uv[1]);

//printf("node %d alpha %e mr %f uv %e %e\n",node,alpha[0],mr[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeXYZ(Grid *grid, int node, double *dxdydz )
{
  double xyzOrig[2];
  double gold;
  double alpha[2], ar[2];
  int i, iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  for(i=0;i<3;i++) xyzOrig[i] = grid->xyz[i+3*node];

  alpha[0] = 0.0;
  for(i=0;i<3;i++) grid->xyz[i+3*node] = xyzOrig[i] + alpha[0]*dxdydz[i];
  gridNodeAR( grid, node, &ar[0] );

  alpha[1] = 1.0e-10;
  for(i=0;i<3;i++) grid->xyz[i+3*node] = xyzOrig[i] + alpha[1]*dxdydz[i];
  gridNodeAR( grid, node, &ar[1] );

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    for(i=0;i<3;i++) grid->xyz[i+3*node] = xyzOrig[i] + alpha[1]*dxdydz[i];
    gridNodeAR( grid, node, &ar[1] );
  }

  for(i=0;i<3;i++) grid->xyz[i+3*node] = xyzOrig[i] + alpha[0]*dxdydz[i];
  
  return grid;
}

Grid *gridSmooth( Grid *grid )
{
  int node;
  double ar, optimizationLimit, laplacianLimit;
  optimizationLimit =0.40;
  laplacianLimit =0.60;
  for (node=0;node<grid->maxnode;node++) {
    if ( gridValidNode( grid, node ) && !gridNodeFrozen( grid, node ) ) {
      gridNodeAR(grid,node,&ar);
      if (ar < optimizationLimit) {
	gridSmoothNode( grid, node );
      }else{
	if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	  gridSmartLaplacian( grid, node ); 
	}
      }
    }
  }
  return grid;
}

Grid *gridSmoothFaceMR( Grid *grid, double optimizationLimit )
{
  int node;
  double mr;
 for (node=0;node<grid->maxnode;node++) {
    if ( gridValidNode( grid, node ) && gridGeometryFace( grid, node )) {
      gridNodeFaceMR(grid,node,&mr);
      if (mr < optimizationLimit) {
	gridSmoothNodeFaceMR( grid, node );
	gridSmoothNodeFaceMR( grid, node );
      }
    }
  }
  return grid;
}

Grid *gridSmoothVolume( Grid *grid )
{
  int node;
  double ar, optimizationLimit, laplacianLimit;
  optimizationLimit =0.30;
  laplacianLimit =0.60;
  for (node=0;node<grid->maxnode;node++) {
    if ( gridValidNode( grid, node ) && !gridGeometryFace( grid, node ) ) {
      gridNodeAR(grid,node,&ar);
      if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	gridSmartLaplacian( grid, node ); 
	gridNodeAR(grid,node,&ar);
      }
      if (ar < optimizationLimit) {
	gridSmoothNode( grid, node );
      }
    }
  }
  return grid;
}

Grid *gridSmartLaplacian(Grid *grid, int node )
{
  double origAR, newAR, origXYZ[3], xyz[3], oneOverNCell;
  AdjIterator it;
  int cell, ncell, inode, ixyz, n;
  
  gridNodeAR(grid, node, &origAR);
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
  ncell =0;

  for ( it = adjFirst(grid->cellAdj,node); adjValid(it) ; it = adjNext(it) ){
    ncell++;
    cell = adjItem(it);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      n = grid->c2n[inode+4*cell];
      for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) xyz[ixyz] += grid->xyz[ixyz+3*n];
    }

  }
  oneOverNCell = 1.0/(double)(ncell*3);
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){  
    xyz[ixyz] -= grid->xyz[ixyz+3*node] * (double)ncell ;
    grid->xyz[ixyz+3*node] = xyz[ixyz] * oneOverNCell;
  }

  gridNodeAR(grid, node, &newAR);
  
  if ( origAR > newAR ) {
    for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) grid->xyz[ixyz+3*node] = origXYZ[ixyz];
    return NULL;
  }

  return grid;
}

Grid *gridSmoothNodeQP(Grid *grid, int node )
{
  int i, minCell, nearestCell;
  double minAR, nearestAR, newAR, searchDirection[3];
  double g00, g01, g11, minRatio, nearestRatio;
  double length, projection;
  double deltaAR, currentAlpha, alpha, lastAlpha;
  double predictedImprovement, actualImprovement, lastImprovement;
  double origXYZ[3], xyz[3];
  bool searchFlag, goodStep;
  int iteration;

  if ( grid != gridNodeXYZ(grid, node, origXYZ)) return NULL;
  if ( grid != gridStoreARDerivative(grid, node ) ) return NULL;

  minAR =2.1;
  minCell = EMPTY;
  for (i=0;i<gridStoreARDegree(grid);i++){
    if (grid->AR[i]<minAR){
      minAR = grid->AR[i];
      minCell = i;
    }
  }

  searchFlag = FALSE;
  if (searchFlag) {
    for (i=0;i<3;i++) searchDirection[i] = grid->dARdX[i+minCell*3];
  }else{
    nearestCell=EMPTY;
    nearestAR = 2.1;
    for (i=0;i<gridStoreARDegree(grid);i++){
      if ( i != minCell){
	if (ABS(grid->AR[i]-minAR)<nearestAR) {
	  nearestCell=i;
	  nearestAR = ABS(grid->AR[i]-minAR);
	}
      }
    }
    if (nearestCell == EMPTY || nearestAR > 0.001 ){
      for (i=0;i<3;i++) searchDirection[i] = grid->dARdX[i+minCell*3];
    }else{
      g00 
	= grid->dARdX[0+minCell*3]*grid->dARdX[0+minCell*3]
	+ grid->dARdX[1+minCell*3]*grid->dARdX[1+minCell*3]
	+ grid->dARdX[2+minCell*3]*grid->dARdX[2+minCell*3];
      g11 
	= grid->dARdX[0+nearestCell*3]*grid->dARdX[0+nearestCell*3]
	+ grid->dARdX[1+nearestCell*3]*grid->dARdX[1+nearestCell*3]
	+ grid->dARdX[2+nearestCell*3]*grid->dARdX[2+nearestCell*3];
      g01 
	= grid->dARdX[0+minCell*3]*grid->dARdX[0+nearestCell*3]
	+ grid->dARdX[1+minCell*3]*grid->dARdX[1+nearestCell*3]
	+ grid->dARdX[2+minCell*3]*grid->dARdX[2+nearestCell*3];
      nearestRatio = (g00-g01)/(g00 + g11 - 2*g01);
      if (nearestRatio > 1.0 || nearestRatio < 0.0 ) nearestRatio = 0.0;
      minRatio = 1.0 - nearestRatio;
      for (i=0;i<3;i++) searchDirection[i] 
			  = minRatio*grid->dARdX[i+minCell*3]
			  + nearestRatio*grid->dARdX[i+nearestCell*3];
      /* reset length to the projection of min cell to search dir*/
      length 
	= searchDirection[0]*searchDirection[0]
	+ searchDirection[1]*searchDirection[1]
	+ searchDirection[2]*searchDirection[2];
      length = sqrt(length);
      for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;
      projection
	= searchDirection[0]*grid->dARdX[0+minCell*3]
	+ searchDirection[1]*grid->dARdX[1+minCell*3]
	+ searchDirection[2]*grid->dARdX[2+minCell*3];
      for (i=0;i<3;i++) searchDirection[i] = projection*searchDirection[i];
      //printf("node %5d min %10.7f near %10.7f\n",node,minRatio,nearestRatio);
    }
  }

  length 
    = searchDirection[0]*searchDirection[0]
    + searchDirection[1]*searchDirection[1]
    + searchDirection[2]*searchDirection[2];
  length = sqrt(length);
  for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;

  alpha = 1.0;
  for (i=0;i<gridStoreARDegree(grid);i++){
    if (i != minCell ) {
      projection
	= searchDirection[0]*grid->dARdX[0+i*3]
	+ searchDirection[1]*grid->dARdX[1+i*3]
	+ searchDirection[2]*grid->dARdX[2+i*3];
      deltaAR = grid->AR[i] - minAR;
      if (ABS(length-projection) < 1e-12){
	currentAlpha=0.0; /* no intersection */
      }else{
	currentAlpha = deltaAR / ( length + projection);
      }
      if (currentAlpha > 0 && currentAlpha < alpha ) alpha = currentAlpha;
    }
  }

  //printf( "node %5d deg %3d active %3d old %12.9f\n",
  //	  node, gridStoreARDegree(grid), minCell, minAR );
  
  goodStep = FALSE;
  actualImprovement = 0.0;
  lastImprovement = -10.0;
  iteration = 0;
  while (alpha > 10.e-10 && !goodStep && iteration < 30 ) {
    iteration++;

    predictedImprovement = length*alpha;
  
    for (i=0;i<3;i++) xyz[i] = origXYZ[i] + alpha*searchDirection[i];
    gridSetNodeXYZ(grid,node,xyz);
    gridNodeAR(grid,node,&newAR);
    actualImprovement = newAR-minAR;
    //printf(" alpha %12.5e predicted %12.9f actual %12.9f new %12.9f\n",
    //	   alpha, predictedImprovement, actualImprovement, newAR);

    if ( actualImprovement > 0.0 && actualImprovement < lastImprovement) {
      for (i=0;i<3;i++) xyz[i] = origXYZ[i] + lastAlpha*searchDirection[i];
      gridSetNodeXYZ(grid,node,xyz);
      gridNodeAR(grid,node,&newAR);
      actualImprovement = newAR-minAR;
      goodStep = TRUE;
    }
    
    if ( actualImprovement > 0.9*predictedImprovement  ){
      goodStep = TRUE;
    }else{
      lastImprovement = actualImprovement;
      lastAlpha = alpha;
      alpha =alpha*0.5;
    }
  }

  //printf( "node %5d deg %3d active %3d old %8.5f new %8.5f\n",
  //  node, gridStoreARDegree(grid), minCell, minAR, newAR );

  if ( actualImprovement <= 0.0  ){
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  if ( newAR > 0.6) return NULL;
  if ( actualImprovement <= 0.000001 ) return NULL;

  return grid;
}
