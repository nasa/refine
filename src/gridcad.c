
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

Grid *gridSafeProject(Grid *grid, int node )
{
  return NULL;
}

