
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
    if ( grid != gridProjectNodeToEdge( grid, node, edgeId ) ) return NULL;
    for ( it = adjFirst(grid->faceAdj,node); adjValid(it); it = adjNext(it) ){
      face = adjItem(it);
      faceId = grid->faceId[face];
      if ( grid != gridSafeProjectNodeToFace( grid, node, faceId ) ) 
	return NULL;
    }
    if ( grid != gridProjectNodeToEdge( grid, node, edgeId ) ) return NULL;    
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

  for (node=0;node<grid->nnode;node++) gridSafeProjectNode(grid,node);  

  return grid;
}
