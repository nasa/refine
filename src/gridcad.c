
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

Grid *gridProjectNodeToEdge(Grid *grid, int node, int edge )
{
  int vol = 1;
  double t, xyznew[3];

  t = DBL_MAX;
  if (!CADGeom_NearestOnEdge( vol, edge, &grid->xyz[3*node], &t, xyznew) ) 
    return NULL;  

  grid->xyz[0+3*node] = xyznew[0];
  grid->xyz[1+3*node] = xyznew[1];
  grid->xyz[2+3*node] = xyznew[2];

  return grid;
}

Grid *gridSafeProject(Grid *grid, int node )
{
  return NULL;
}

