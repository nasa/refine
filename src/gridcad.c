
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

  for (node=0;node<grid->nnode;node++) 
    if ( gridSafeProjectNode(grid,node) != grid ) notProjected++;

  if (notProjected > 0){
    printf("gridProject: %d of %d nodes not projected.\n",
	   notProjected,grid->nnode);
    return NULL;
  }

  return grid;
}

Grid *gridSmoothNode(Grid *grid, int node )
{
  double xyz[3], xyzProj[3], uv[2];
  double ar, dARdx[4];
  double du[3], dv[3];
  double dARdu[2];
  int vol =1;
  int face, faceId;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryEdge( grid, node ) ) return grid;
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(grid->faceAdj, node));
    faceId = grid->faceId[face];
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    gridNodeXYZ( grid, node, xyz);
    gridNodeUV( grid, node, faceId, uv);
    if ( !CADGeom_PointOnFace( vol, faceId,   
			       uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
      printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );

    dARdu[0] = dARdx[0]*du[0] + dARdx[1]*du[1] + dARdx[2]*du[2] ; 
    dARdu[1] = dARdx[0]*dv[0] + dARdx[1]*dv[1] + dARdx[2]*dv[2] ; 

    return gridOptimizeUV( grid, node, dARdu );
  }    
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

Grid *gridSmooth( Grid *grid )
{
  int node;
  double ar,limit;
  limit =0.99;
  for (node=0;node<grid->nnode;node++) {
    gridNodeAR(grid,node,&ar);
    if (ar < 0.99) gridSmoothNode( grid, node );
  }
  for (node=0;node<grid->nnode;node++) {
    gridNodeAR(grid,node,&ar);
    if (ar < 0.99) gridSmoothNode( grid, node );
  }
  return grid;
}

