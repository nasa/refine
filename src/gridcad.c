
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
#include "gridcad.h"
#include "adj.h"
#include "gridStruct.h"

Grid *gridProjectNodeToEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyznew[3];

  if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL;

  if (!CADGeom_NearestOnEdge( vol, edgeId, &grid->xyz[3*node], &t, xyznew) ) 
    return NULL;  

  if ( grid != gridSetNodeT( grid, node, edgeId, t ) ) return NULL;

  grid->xyz[0+3*node] = xyznew[0];
  grid->xyz[1+3*node] = xyznew[1];
  grid->xyz[2+3*node] = xyznew[2];

  return grid;
}

Grid *gridProjectNodeToFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyznew[3];

  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;

  if (!CADGeom_NearestOnFace( vol, faceId, &grid->xyz[3*node], uv, xyznew) ) 
    return NULL;  

  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  grid->xyz[0+3*node] = xyznew[0];
  grid->xyz[1+3*node] = xyznew[1];
  grid->xyz[2+3*node] = xyznew[2];

  return grid;
}

Grid *gridSafeProjectNode(Grid *grid, int node )
{
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(grid->edgeAdj, node));
    edgeId = grid->edgeId[edge];
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId ) ) return NULL;
    for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
      face = adjItem(it);
      faceId = grid->faceId[face];
      if ( grid != gridSafeProjectNodeToFace( grid, node, faceId ) ) 
	return NULL;
    }
    if ( grid != gridSafeProjectNodeToEdge( grid, node, edgeId ) ) return NULL;
    return grid;
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(grid->faceAdj, node));
    faceId = grid->faceId[face];
    if ( grid != gridSafeProjectNodeToFace( grid, node, faceId ) ) return NULL;
    return grid;
  }

  return grid;
}

Grid *gridSafeProjectNodeToEdge(Grid *grid, int node, int edgeId )
{
  double xyz[3], t;

  xyz[0] = grid->xyz[0+3*node];
  xyz[1] = grid->xyz[1+3*node];
  xyz[2] = grid->xyz[2+3*node];
  if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL; 
  if ( grid != gridProjectNodeToEdge( grid, node, edgeId ) ) return NULL;
  if ( gridNegCellAroundNode( grid, node ) ) {
    grid->xyz[0+3*node] = xyz[0];
    grid->xyz[1+3*node] = xyz[1];
    grid->xyz[2+3*node] = xyz[2];
    gridSetNodeT( grid, node, edgeId, t );
    return NULL;
  }
  return grid;
}

Grid *gridSafeProjectNodeToFace(Grid *grid, int node, int faceId )
{
  double xyz[3], uv[2];

  xyz[0] = grid->xyz[0+3*node];
  xyz[1] = grid->xyz[1+3*node];
  xyz[2] = grid->xyz[2+3*node];
  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL; 
  if ( grid != gridProjectNodeToFace( grid, node, faceId ) ) return NULL;
  if ( gridNegCellAroundNode( grid, node ) ) {
    grid->xyz[0+3*node] = xyz[0];
    grid->xyz[1+3*node] = xyz[1];
    grid->xyz[2+3*node] = xyz[2];
    gridSetNodeUV( grid, node, faceId, uv[0], uv[1] );
    return NULL;
  }
  return grid;
}

Grid *gridProject(Grid *grid)
{
  int node;
  int notProjected;
  notProjected = 0;

  for (node=0;node<grid->maxnode;node++)
    if ( gridValidNode( grid, node ) )
      if ( gridSafeProjectNode(grid,node) != grid ) notProjected++;

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
    if ( gridValidNode( grid, node ) ) 
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
  int i, j, nodes[4], level2nodes[4];
  AdjIterator it, level2;

  if ( !gridValidNode( grid, node ) ) return NULL;
  
  if ( gridSafeProjectNode(grid,node) != grid ) {
    for ( it = adjFirst(grid->cellAdj,node); 
	  adjValid(it); 
	  it = adjNext(it) ){
      gridCell(grid, adjItem(it), nodes);
      for (i=0;i<4;i++)
	if (!gridGeometryFace( grid, nodes[i])) 
	  gridSmoothNode( grid, nodes[i]);
    }      
    gridSwapNearNode( grid, node);
    for ( it = adjFirst(grid->cellAdj,node); 
	  adjValid(it); 
	  it = adjNext(it) ){
      gridCell(grid, adjItem(it), nodes);
      for (i=0;i<4;i++)
	if (!gridGeometryFace( grid, nodes[i])) 
	  gridSmoothNode( grid, nodes[i]);
    }      
    if ( gridSafeProjectNode(grid,node) != grid ) {
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
      if ( gridSafeProjectNode(grid,node) != grid )return NULL;
    }
  }

  return grid;
}

Grid *gridSmoothNode(Grid *grid, int node )
{
  double xyzProj[3], uv[2], t;
  double ar, dARdx[3];
  double du[3], dv[3], dt[3];
  double dARdu[2], dARdt;
  int vol =1;
  int face, faceId;
  int edge, edgeId;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(grid->edgeAdj, node));
    edgeId = grid->edgeId[edge];
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    gridNodeT( grid, node, edgeId, &t);
    if ( !CADGeom_PointOnEdge( vol, edgeId,   
			       t, xyzProj, 1, dt, NULL) )
      printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
    dARdt = dARdx[0]*dt[0] + dARdx[1]*dt[1] + dARdx[2]*dt[2];
    return gridOptimizeT( grid, node, dARdt );
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(grid->faceAdj, node));
    faceId = grid->faceId[face];
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    gridNodeUV( grid, node, faceId, uv);
    if ( !CADGeom_PointOnFace( vol, faceId,   
			       uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
      printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );

    dARdu[0] = dARdx[0]*du[0] + dARdx[1]*du[1] + dARdx[2]*du[2] ; 
    dARdu[1] = dARdx[0]*dv[0] + dARdx[1]*dv[1] + dARdx[2]*dv[2] ; 

    return gridOptimizeUV( grid, node, dARdu );
  }    
  gridNodeARDerivative ( grid, node, &ar, dARdx);
  return gridOptimizeXYZ( grid, node, dARdx );
}

Grid *gridOptimizeT(Grid *grid, int node, double dt )
{
  int vol =1;
  double tOrig, t;
  int edge, edgeId;
  int face, faceId;
  AdjIterator it;
  double gold;
  double alpha[2], ar[2];
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

  alpha[1] = 1.0e-10;
  t = tOrig + alpha[1]*dt;
  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[1] );

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
  }

  t = tOrig + alpha[0]*dt;
  if ( !CADGeom_PointOnEdge( vol, edgeId, t, &grid->xyz[3*node], 
			     0, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
  gridSetNodeT(grid, node, edgeId, t);

  for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
    face = adjItem(it);
    faceId = grid->faceId[face];
    if ( grid != gridSafeProjectNodeToFace( grid, node, faceId ) ) 
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
  double alpha[2], ar[2];
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

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if ( !CADGeom_PointOnFace( vol, faceId, uv, &grid->xyz[3*node], 
			     0, NULL, NULL, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  gridNodeAR( grid, node, &ar[1] );

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
  optimizationLimit =0.30;
  laplacianLimit =0.60;
  for (node=0;node<grid->maxnode;node++) {
    if ( gridValidNode( grid, node ) ) {
      gridNodeAR(grid,node,&ar);
      if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	gridSmartLaplacian( grid, node ); 
	gridNodeAR(grid,node,&ar);
      }
      if (ar < optimizationLimit) {
	gridSmoothNode( grid, node );
	gridSmoothNode( grid, node );
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
