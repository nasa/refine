
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
#ifdef HAVE_SDK
#include "CADGeom/CADGeom.h"
#else
#include "FAKEGeom.h"
#endif
#include "gridmetric.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridcad.h"

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

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: CADGeom_NearestOnFace failed.",__FILE__,__LINE__);
    return NULL;  
  }

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

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: CADGeom_NearestOnFace failed.",__FILE__,__LINE__);
    return NULL;  
  }

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;
  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  return grid;
}

Grid *gridSetUVofFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: CADGeom_NearestOnFace failed.",__FILE__,__LINE__);
    return NULL;
  }

  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  return grid;
}

Grid *gridEvaluateEdgeAtT(Grid *grid, int node, double t )
{
  int vol =1;
  int nodes[2];
  int edge, edgeId;
  double xyz[3];

  if (gridGeometryNode(grid,node)) return NULL;
  if ( grid != gridNodeXYZ(grid, node, xyz) ) return NULL;

  edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
  if ( grid != gridEdge(grid,edge,nodes,&edgeId)) return NULL;

  if ( !CADGeom_PointOnEdge( vol, edgeId, t, xyz, 
			     0, NULL, NULL) ) {
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  gridSetNodeT(grid, node, edgeId, t);
  gridSetNodeXYZ(grid, node, xyz);
  return grid;
}

Grid *gridEvaluateFaceAtUV(Grid *grid, int node, double *uv )
{
  int vol =1;
  int nodes[3];
  int face, faceId;
  double xyz[3];

  if ( gridGeometryEdge(grid, node) ) return NULL;
  if ( grid != gridNodeXYZ(grid, node, xyz) ) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if (EMPTY == face) return NULL;

  gridFace(grid,face,nodes,&faceId);

  if ( !CADGeom_PointOnFace( vol, faceId, uv, xyz, 
			     0, NULL, NULL, NULL, NULL, NULL) ){
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  gridSetNodeUV(grid, node, faceId, uv[0], uv[1]);
  gridSetNodeXYZ(grid, node, xyz);
  return grid;
}

Grid *gridUpdateFaceParameter(Grid *grid, int node ){

  AdjIterator it;
  int nodes[3];
  int face, faceId;

  for ( it = adjFirst(gridFaceAdj(grid),node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    if ( grid != gridFace(grid, face, nodes, &faceId) )  return NULL;
    if ( grid != gridSafeProjectNodeToFace( grid, node, faceId, 1.0 ) ) 
      return NULL;
  }

  return grid;
}

Grid *gridProjectToFace(Grid *grid, int faceId, 
			double *xyz, double *uv, double *newxyz )
{
  int vol = 1;

  if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, newxyz) ) {
    printf("%s: %d: CADGeom_NearestOnFace failed.",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridFaceNormalAtUV(Grid *grid, int faceId,
			 double *uv, double *xyz, double *normal )

{
  int vol =1;

  if ( !CADGeom_NormalToFace( vol, faceId, uv, xyz, normal) ){
    printf ( "ERROR: CADGeom_NormalToFace, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  return grid;
}

Grid *gridFaceNormalAtXYZ(Grid *grid, int faceId, double *xyz, double *normal )
{
  double uv[2], newxyz[3];

  uv[0] = DBL_MAX;
  uv[1] = DBL_MAX;
  if (grid != gridProjectToFace(grid,faceId,xyz,uv,newxyz)) return NULL;
  if (grid != gridFaceNormalAtUV(grid, faceId, uv, newxyz, normal)) return NULL;

  return grid;
}

Grid *gridSafeProjectNode(Grid *grid, int node, double ratio )
{
  int nodes[3];
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
    gridEdge(grid, edge, nodes, &edgeId );
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId, ratio ) ) 
      return NULL;
    for ( it = adjFirst(gridFaceAdj(grid),node); adjValid(it); it = adjNext(it) ){
      face = adjItem(it);
      gridFace(grid, face, nodes, &faceId );
      if ( grid != gridSetUVofFace( grid, node, faceId ) ) 
	return NULL;
    }
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId, 1.0 ) ) return NULL;
    return grid;
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(gridFaceAdj(grid), node));
    gridFace(grid, face, nodes, &faceId );
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
  double xyz[3];
  int attempts;

  if ( grid != gridNodeXYZ( grid, node, origxyz ) ) return NULL; 
  if ( grid != gridNodeT( grid, node, edgeId, &origt ) ) return NULL; 

  if ( grid != gridProjectNodeToEdge( grid, node, edgeId ) ) return NULL;

  if (!gridNegCellAroundNode( grid, node )) return grid;

  gridSetNodeT( grid, node, edgeId, origt );

  if ( grid != gridNodeXYZ( grid, node, projxyz ) ) return NULL; 
  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;

  attempts = 0;
  while ( attempts < 10 && gridNegCellAroundNode( grid, node ) ) {
    attempts++;
    xyz[0] = (1.0-ratio)*xyz[0] + ratio*origxyz[0];
    xyz[1] = (1.0-ratio)*xyz[1] + ratio*origxyz[1];
    xyz[2] = (1.0-ratio)*xyz[2] + ratio*origxyz[2];
    gridSetNodeXYZ( grid, node, xyz );
  }

  if (gridNegCellAroundNode( grid, node )) 
    gridSetNodeXYZ( grid, node, origxyz );

  return NULL;
}

Grid *gridSafeProjectNodeToFace(Grid *grid, int node, int faceId, double ratio )
{
  double origxyz[3], origuv[2];
  double projxyz[3];
  double xyz[3];
  int attempts;

  if ( grid != gridNodeXYZ( grid, node, origxyz ) ) return NULL; 
  if ( grid != gridNodeUV( grid, node, faceId, origuv ) ) return NULL; 

  if ( grid != gridProjectNodeToFace( grid, node, faceId ) ) return NULL;

  if (!gridNegCellAroundNode( grid, node )) return grid;

  gridSetNodeUV( grid, node, faceId, origuv[0], origuv[1] );

  if ( grid != gridNodeXYZ( grid, node, projxyz ) ) return NULL; 
  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;

  attempts = 0;
  while ( attempts < 10 && gridNegCellAroundNode( grid, node ) ) {
    attempts++;
    xyz[0] = (1.0-ratio)*xyz[0] + ratio*origxyz[0];
    xyz[1] = (1.0-ratio)*xyz[1] + ratio*origxyz[1];
    xyz[2] = (1.0-ratio)*xyz[2] + ratio*origxyz[2];
    gridSetNodeXYZ( grid, node, xyz );
  }

  if (gridNegCellAroundNode( grid, node )) 
    gridSetNodeXYZ( grid, node, origxyz );

  return NULL;
}

Grid *gridProject(Grid *grid)
{
  int node;
  int notProjected;
  notProjected = 0;

  for (node=0;node<gridMaxNode(grid);node++)
    if ( gridValidNode( grid, node ) )
      if ( gridSafeProjectNode(grid,node,1.0) != grid ) notProjected++;

  if (notProjected > 0){
    // printf("gridProject: %d of %d nodes not projected.\n", notProjected,gridNNode(grid));
    return NULL;
  }

  return grid;
}

Grid *gridNodeProjectionDisplacement(Grid *grid, int node, double *displacement )
{
  int nodes[3];
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;
  int vol = 1;
  double t, uv[2], xyz[3], xyznew[3];

  displacement[0] = displacement[1] = displacement[2] = 0.0;

  if ( gridGeometryNode( grid, node ) ) return grid;

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;

  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
    gridEdge(grid, edge, nodes, &edgeId );
    if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL;
    if (!CADGeom_NearestOnEdge( vol, edgeId, xyz, &t, xyznew) ) return NULL;  
    if ( grid != gridSetNodeT( grid, node, edgeId, t ) ) return NULL;
    for ( it = adjFirst(gridFaceAdj(grid),node); 
	  adjValid(it); 
	  it = adjNext(it) ){
      face = adjItem(it);
      gridFace(grid, face, nodes, &faceId );
      if ( grid != gridSetUVofFace( grid, node, faceId ) ) return NULL;
    }
    gridSubtractVector(xyznew,xyz,displacement);
    return grid;
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(gridFaceAdj(grid), node));
    gridFace(grid, face, nodes, &faceId );
    if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;
    if (!CADGeom_NearestOnFace( vol, faceId, xyz, uv, xyznew) )  {
      printf("%s: %d: CADGeom_NearestOnFace failed.",__FILE__,__LINE__);
      return NULL;  
    }
    gridSetNodeUV(grid, node, faceId, uv[0], uv[1]);
    gridSubtractVector(xyznew,xyz,displacement);
    return grid;
  }

  return grid;
}

Grid *gridRobustProject(Grid *grid)
{
  int  node;
  int notProjected;
  notProjected = 0;
  for (node=0;node<gridMaxNode(grid);node++)
    if ( gridValidNode( grid, node ) && !gridNodeFrozen(grid, node) ) 
      if (gridRobustProjectNode( grid, node)!= grid ) 
	notProjected++;

  if (notProjected > 0){
    printf("gridRobustProject: %d of %d nodes not projected.\n",
	   notProjected,gridNNode(grid));
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
  
  if ( gridSafeProjectNode(grid,node,0.95) == grid ) return grid;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++)
      if (!gridGeometryFace( grid, nodes[i])) 
	gridSmoothNode( grid, nodes[i]);
  }      
  gridSwapNearNodeExceptBoundary( grid, node);
  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++)
      if (!gridGeometryFace( grid, nodes[i])) 
	gridSmoothNode( grid, nodes[i]);
  }
      
  if ( gridSafeProjectNode(grid,node,0.95) == grid ) return grid;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++) {
      for ( level2 = adjFirst(gridCellAdj(grid),nodes[i]); 
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

  if ( gridSafeProjectNode(grid,node,1.0) == grid ) return grid;
  
  gridNodeXYZ(grid,node,xyz);
  printf(" try to c-p %d X %10.5f Y %10.5f Z %10.5f ...",
	 node,xyz[0],xyz[1],xyz[2]);
  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++) {
      goodnode = nodes[i];
      if ( node != goodnode && 
	   gridGeometryFace( grid, goodnode) &&
	   grid == gridSafeProjectNode( grid, goodnode, 1.0 ) ) { 
	if ( grid == gridCollapseEdge(grid, NULL, goodnode, node, 0.0 ) ){
	  printf(" got it ! %d\n",goodnode);
	  return grid;
	}
      }
    }
  }
  printf("\n");

  return NULL;
}

Grid *gridSmoothNearNode1(Grid *grid, int node )
{
#define SMOOTHDEG (500)
  int i, nodes[4];
  int nlist, smooth, look, nodelist[SMOOTHDEG];
  GridBool looking;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return NULL;
  
  nlist =0;
  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++) {
      if (!gridNodeFrozen( grid, nodes[i])) {
	looking = (nlist<=SMOOTHDEG);
	look = 0;
	for (look=0;look<nlist && looking ; look++){
	  looking = (nodelist[look] != nodes[i]);
	}
	if (looking && nlist<=SMOOTHDEG){
	  nodelist[nlist] = nodes[i];
	  nlist++;
	}
      }
    }
  }      

  for (smooth=0;smooth<5;smooth++)
    for (i=0;i<nlist;i++) 
      gridSmoothNode( grid, nodelist[nlist]);

  return grid;
}

Grid *gridSmoothNearNode(Grid *grid, int node )
{
  int i, i0, nodes[4], nodes0[4];
  int nlist, smooth, look, nodelist[SMOOTHDEG];
  GridBool looking;
  AdjIterator it, it0;

  if (!gridValidNode(grid,node)) return NULL;

  nlist =0;
  for ( it0 = adjFirst(gridCellAdj(grid),node); 
	adjValid(it0); 
	it0 = adjNext(it0) ){
    gridCell(grid, adjItem(it0), nodes0);
    for (i0=0;i0<4;i0++) {
      for ( it = adjFirst(gridCellAdj(grid),nodes0[i0]); 
	    adjValid(it); 
	    it = adjNext(it) ){
	gridCell(grid, adjItem(it), nodes);
	for (i=0;i<4;i++) {
	  if (!gridNodeFrozen( grid, nodes[i])) {
	    looking = (nlist<=SMOOTHDEG);
	    look = 0;
	    for (look=0;look<nlist && looking ; look++){
	      looking = (nodelist[look] != nodes[i]);
	    }
	    if (looking && nlist<=SMOOTHDEG){
	      nodelist[nlist] = nodes[i];
	      nlist++;
	    }
	  }
	}
      }
    }
  }      

  for (smooth=0;smooth<2;smooth++)
    for (i=0;i<nlist;i++) 
      gridSmoothNode( grid, nodelist[nlist]);

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
  int nodes[3];
  int face, faceId;
  int edge, edgeId;
  int maxsmooth;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryBetweenFace( grid, node ) &&
       !gridGeometryEdge( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
    gridEdge(grid,edge,nodes,&edgeId);
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
      face = adjItem(adjFirst(gridFaceAdj(grid), node));
      gridFace(grid,face,nodes,&faceId);
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
  double xyzProj[3], uv[2];
  double mr, dMRdx[3];
  double du[3], dv[3];
  double dMRdu[2];
  int vol =1;
  int nodes[3];
  int face, faceId;

  if ( gridGeometryEdge( grid, node ) ) return grid;
  if ( !gridGeometryFace( grid, node ) ) return grid;


  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  gridFace(grid,face,nodes,&faceId);
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
  double tOrig, t;
  double gold;
  double alpha[2], ar[2], mr;
  int iter;
  int nodes[2];
  int edge, edgeId;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
  if (EMPTY == edge) return NULL;
  if ( grid != gridEdge(grid,edge,nodes,&edgeId) ) return NULL;
  gridNodeT( grid, node, edgeId, &tOrig);

  alpha[0] = 0.0;
  t = tOrig + alpha[0]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
  gridNodeAR( grid, node, &ar[0] );
  gridNodeFaceMR( grid, node, &mr );
  ar[0] = MIN(mr, ar[0]);

  alpha[1] = 1.0e-10;
  t = tOrig + alpha[1]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
  gridNodeAR( grid, node, &ar[1] );
  gridNodeFaceMR( grid, node, &mr );
  ar[1] = MIN(mr, ar[1]);

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    t = tOrig + alpha[1]*dt;
    if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
    gridNodeAR( grid, node, &ar[1] );
    gridNodeFaceMR( grid, node, &mr );
    ar[1] = MIN(mr, ar[1]);
  }

  t = tOrig + alpha[0]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
  if (grid != gridUpdateFaceParameter(grid, node )) return NULL;

  return grid;
}

Grid *gridOptimizeUV(Grid *grid, int node, double *dudv )
{
  double uvOrig[2], uv[2];
  int nodes[3];
  int face, faceId;
  double gold;
  double alpha[2], ar[2], mr;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeAR( grid, node, &ar[0] );
  gridNodeFaceMR( grid, node, &mr );
  ar[0] = MIN(mr, ar[0]);

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
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
    if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
    gridNodeAR( grid, node, &ar[1] );
    gridNodeFaceMR( grid, node, &mr );
    ar[1] = MIN(mr, ar[1]);
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;

//printf("node %d alpha %e ar %f uv %e %e\n",node,alpha[0],ar[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeFaceUV(Grid *grid, int node, double *dudv )
{
  double uvOrig[2], uv[2];
  int nodes[3];
  int face, faceId;
  double gold;
  double alpha[2], mr[2], ar;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  
  gridNodeFaceMR( grid, node, &mr[0] );

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeFaceMR( grid, node, &mr[1] );

  iter = 0;
  gridNodeAR( grid, node, &ar ); 
  while ( mr[1] > mr[0] && mr[1] > 0.0 && iter < 100 && ar > 0.1){
    iter++;
    alpha[0] = alpha[1]; mr[0] = mr[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];
    if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
    gridNodeFaceMR( grid, node, &mr[1] );
    gridNodeAR( grid, node, &ar );
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;

//printf("node %d alpha %e mr %f uv %e %e\n",node,alpha[0],mr[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeXYZ(Grid *grid, int node, double *dxdydz )
{
  double xyzOrig[3];
  double xyz[3];
  double gold;
  double alpha[2], ar[2];
  int i, iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  if (grid != gridNodeXYZ(grid,node,xyzOrig)) return NULL;

  alpha[0] = 0.0;
  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[0]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR( grid, node, &ar[0] );

  alpha[1] = 1.0e-10;
  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[1]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR( grid, node, &ar[1] );

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[1]*dxdydz[i];
    gridSetNodeXYZ(grid,node,xyz);
    gridNodeAR( grid, node, &ar[1] );
  }

  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[0]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  
  return grid;
}

Grid *gridSmooth( Grid *grid )
{
  int node;
  double ar, optimizationLimit, laplacianLimit;
  optimizationLimit =0.40;
  laplacianLimit =0.60;
  for (node=0;node<gridMaxNode(grid);node++) {
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
 for (node=0;node<gridMaxNode(grid);node++) {
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
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridGeometryFace( grid, node ) &&
	 gridNodeLocal(grid,node) ) {
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
  double origAR, newAR;
  double origXYZ[3], xyz[3], nodeXYZ[3];
  double oneOverNCell;
  AdjIterator it;
  int nodes[4], ncell, inode, ixyz;
  
  gridNodeAR(grid, node, &origAR);
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
  ncell =0;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it) ; 
	it = adjNext(it) ){
    ncell++;
    gridCell(grid,adjItem(it),nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      gridNodeXYZ(grid,nodes[inode],nodeXYZ);
      for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) xyz[ixyz] += nodeXYZ[ixyz];
    }
  }
  oneOverNCell = 1.0/(double)(ncell*3);
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){  
    xyz[ixyz] -= origXYZ[ixyz] * (double)ncell ;
    xyz[ixyz] = xyz[ixyz] * oneOverNCell;
  }
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR(grid, node, &newAR);
  
  if ( origAR > newAR ) {
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  return grid;
}

Grid *gridSmartVolumeLaplacian(Grid *grid, int node )
{
  double origVol, newVol;
  double origXYZ[3], xyz[3], nodeXYZ[3];
  double oneOverNCell;
  AdjIterator it;
  int nodes[4], ncell, inode, ixyz;
  
  gridNodeVolume(grid, node, &origVol);
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
  ncell =0;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it) ; 
	it = adjNext(it) ){
    ncell++;
    gridCell(grid,adjItem(it),nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      gridNodeXYZ(grid,nodes[inode],nodeXYZ);
      for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) xyz[ixyz] += nodeXYZ[ixyz];
    }
  }
  oneOverNCell = 1.0/(double)(ncell*3);
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){  
    xyz[ixyz] -= origXYZ[ixyz] * (double)ncell ;
    xyz[ixyz] = xyz[ixyz] * oneOverNCell;
  }
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR(grid, node, &newVol);
  
  if ( origVol > newVol ) {
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  return grid;
}

Grid *gridSmoothNodeQP(Grid *grid, int node )
{
  int i, minCell, nearestCell;
  double minAR, nearestAR, nearestDifference, newAR, searchDirection[3];
  double g00, g01, g11, minRatio, nearestRatio;
  double length, projection;
  double deltaAR, currentAlpha, alpha, lastAlpha;
  double predictedImprovement, actualImprovement, lastImprovement;
  double minDirection[3], nearestDirection[3], dARdX[3];
  double origXYZ[3], xyz[3];
  GridBool searchFlag, goodStep;
  int iteration;

  if ( grid != gridNodeXYZ(grid, node, origXYZ)) return NULL;
  if ( grid != gridStoreARDerivative(grid, node ) ) return NULL;

  minAR =2.1;
  minCell = EMPTY;
  for (i=0;i<gridStoredARDegree(grid);i++){
    if (gridStoredAR(grid,i)<minAR){
      minAR = gridStoredAR(grid,i);
      minCell = i;
    }
  }

  searchFlag = FALSE;
  if (searchFlag) {
    gridStoredARDerivative(grid, minCell, searchDirection);
  }else{
    nearestCell=EMPTY;
    nearestAR = 2.1;
    for (i=0;i<gridStoredARDegree(grid);i++){
      if ( i != minCell){
	nearestDifference = ABS(gridStoredAR(grid,i)-minAR);
	if (nearestDifference<nearestAR) {
	  nearestCell=i;
	  nearestAR = nearestDifference;
	}
      }
    }
    if (nearestCell == EMPTY || nearestAR > 0.001 ){
      gridStoredARDerivative(grid, minCell, searchDirection);
      gridStoredARDerivative(grid, minCell, minDirection);
    }else{
      gridStoredARDerivative(grid, minCell, minDirection);
      gridStoredARDerivative(grid, nearestCell, nearestDirection);
      g00 = gridDotProduct(minDirection,minDirection);
      g11 = gridDotProduct(nearestDirection,nearestDirection);
      g01 = gridDotProduct(minDirection,nearestDirection);
      nearestRatio = (g00-g01)/(g00 + g11 - 2*g01);
      if (nearestRatio > 1.0 || nearestRatio < 0.0 ) nearestRatio = 0.0;
      minRatio = 1.0 - nearestRatio;
      for (i=0;i<3;i++) searchDirection[i] 
			  = minRatio*minDirection[i]
			  + nearestRatio*nearestDirection[i];
      /* reset length to the projection of min cell to search dir*/
      length = sqrt(gridDotProduct(searchDirection,searchDirection));
      for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;
      projection = gridDotProduct(searchDirection,minDirection);
      for (i=0;i<3;i++) searchDirection[i] = projection*searchDirection[i];
      //printf("node %5d min %10.7f near %10.7f\n",node,minRatio,nearestRatio);
    }
  }

  length = sqrt(gridDotProduct(searchDirection,searchDirection));
  for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;

  alpha = 1.0;
  for (i=0;i<gridStoredARDegree(grid);i++){
    if (i != minCell ) {
      gridStoredARDerivative(grid,i,dARdX);
      projection = gridDotProduct(searchDirection,dARdX);
      deltaAR = gridStoredAR(grid,i) - minAR;
      if (ABS(length-projection) < 1e-12){
	currentAlpha=0.0; /* no intersection */
      }else{
	currentAlpha = deltaAR / ( length + projection);
      }
      if (currentAlpha > 0 && currentAlpha < alpha ) alpha = currentAlpha;
    }
  }

  //printf( "node %5d deg %3d active %3d old %12.9f\n",
  //	  node, gridStoredARDegree(grid), minCell, minAR );
  
  goodStep = FALSE;
  actualImprovement = 0.0;
  lastImprovement = -10.0;
  lastAlpha = alpha;
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
  //  node, gridStoredARDegree(grid), minCell, minAR, newAR );

  if ( actualImprovement <= 0.0  ){
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  if ( newAR > 0.6) return NULL;
  if ( actualImprovement <= 0.000001 ) return NULL;

  return grid;
}

static double reflect( Grid *grid,
		       double simplex[4][3], double volume[4], double avgXYZ[3],
		       int node, int worst, double factor );

static Grid *gridMakeFacesFromSimplex(Grid *grid, 
				      double simplex[4][3], int faceId);

Grid *gridSmoothNodeVolume( Grid *grid, int node )
{
  int evaluations;
  int s, i;
  double origXYZ[3], avgXYZ[3];
  double simplex[4][3];
  double volume[4];
  double lengthScale;
  int best, worst, secondworst; 
  double newVolume, savedVolume;
  GridBool makefaces = FALSE;
  int faceId = 1;

  gridSmartVolumeLaplacian( grid, node );

  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  lengthScale = gridAverageEdgeLength(grid, node );

  for(s=0;s<4;s++)
    for(i=0;i<3;i++)
      simplex[s][i] = origXYZ[i];

  simplex[1][0] += lengthScale;
  simplex[2][1] += lengthScale;
  simplex[3][2] += lengthScale;

  for(s=0;s<4;s++) {
    gridSetNodeXYZ(grid, node, simplex[s]);
    gridNodeVolume(grid,node,&volume[s]);
  }

  for(i=0;i<3;i++) avgXYZ[i] = 0.0;
  for(s=0;s<4;s++)
    for(i=0;i<3;i++) avgXYZ[i] += simplex[s][i];

  evaluations = 4;
  while (evaluations < 200 ) {


    best = 0;
    if ( volume[0] > volume[1] ) {
      secondworst = 0;
      worst = 1;
    }else{
      secondworst = 1;
      worst = 0;
    }
    
    for(s=0;s<4;s++) {
      if (volume[s]>=volume[best]) best = s;
      if (volume[s]<volume[worst]) {
	secondworst = worst;
	worst = s;
      }else{
	if ( s!=worst && volume[s]<volume[secondworst]) secondworst = s;
      }
    }

    printf("evaluations%6d best%20.15f worst%20.15f\n", 
	   evaluations, volume[best], volume[worst]);
    if (makefaces) gridMakeFacesFromSimplex(grid, simplex, ++faceId);

    if (volume[best]-volume[worst] < 1.0e-5*volume[best]) break;

    evaluations++;
    newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, -1.0 );
    if ( newVolume >= volume[best] ) {
      evaluations++;
      newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, 2.0 );
    } else {
      if (newVolume <= volume[secondworst]) {
	savedVolume = volume[worst];
	evaluations++;
	newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, 0.5 );
	if (newVolume <= savedVolume) {
	  for(s=0;s<4;s++) {
	    if (s != best) {
	      for(i=0;i<3;i++) 
		simplex[s][i]=0.5*(simplex[s][i]+simplex[best][i]);
	      gridSetNodeXYZ(grid, node, simplex[s]);
	      gridNodeVolume(grid,node,&volume[s]);
	    }
	  }
	}      
      }
    }
  }    

  best = 0;
  for(s=1;s<4;s++) if (volume[s]>=volume[best]) best = s;

  gridSetNodeXYZ(grid, node, simplex[best]);

  return grid;
}

static double reflect( Grid *grid,
		       double simplex[4][3], double volume[4], double avgXYZ[3],
		       int node, int worst, double factor)
{
  int i;
  double factor1, factor2;
  double reflectedXYZ[3];
  double reflectedVolume;

  factor1 = (1.0-factor) / 3.0;
  factor2 = factor1 - factor;

  for(i=0;i<3;i++) 
    reflectedXYZ[i] = factor1*avgXYZ[i] - factor2*simplex[worst][i];

  gridSetNodeXYZ(grid,node,reflectedXYZ );
  gridNodeVolume(grid,node,&reflectedVolume);

  if ( reflectedVolume > volume[worst] ) {
    volume[worst] = reflectedVolume;
    for(i=0;i<3;i++) avgXYZ[i] += ( reflectedXYZ[i] - simplex[worst][i] );
    for(i=0;i<3;i++) simplex[worst][i] = reflectedXYZ[i];
  }

  return reflectedVolume;
}

static Grid *gridMakeFacesFromSimplex(Grid *grid, 
				      double simplex[4][3], int faceId)
{
  int node;
  int nodes[4];
  for(node=0;node<4;node++)
    nodes[node]=gridAddNode(grid,
			    simplex[node][0],
			    simplex[node][1],
			    simplex[node][2]);
  gridAddFace(grid,nodes[0],nodes[1],nodes[2],faceId);
  gridAddFace(grid,nodes[0],nodes[1],nodes[3],faceId);
  gridAddFace(grid,nodes[1],nodes[2],nodes[3],faceId);
  gridAddFace(grid,nodes[0],nodes[3],nodes[2],faceId);
  return grid;
}
