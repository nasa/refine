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

#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridedger.h"

GridEdger *gridedgerCreate( Grid *grid, int edgeId )
{
  GridEdger *ge;
  ge = malloc(sizeof(GridEdger));
  ge->grid = grid;

  ge->edgeId = edgeId;

  gridAttachPacker( grid, gridedgerPack, (void *)ge );
  gridAttachNodeSorter( grid, gridedgerSortNode, (void *)ge );
  gridAttachReallocator( grid, gridedgerReallocator, (void *)ge );
  gridAttachFreeNotifier( grid, gridedgerGridHasBeenFreed, (void *)ge );

  return ge;
}

Grid *gridedgerGrid(GridEdger *ge)
{
  return ge->grid;
}

void gridedgerFree(GridEdger *ge)
{
  if (NULL != ge->grid) { 
    gridDetachPacker( ge->grid );
    gridDetachNodeSorter( ge->grid );
    gridDetachReallocator( ge->grid );
    gridDetachFreeNotifier( ge->grid );
  }
  free(ge);
}

void gridedgerPack(void *voidGridEdger, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n)
{
  //GridEdger *ge = (GridEdger *)voidGridEdger;
}

void gridedgerSortNode(void *voidGridEdger, int maxnode, int *o2n)
{
  //GridEdger *ge = (GridEdger *)voidGridEdger;
}

void gridedgerReallocator(void *voidGridEdger, int reallocType, 
			 int lastSize, int newSize)
{
  //GridEdger *ge = (GridEdger *)voidGridEdger;
}

void gridedgerGridHasBeenFreed(void *voidGridEdger )
{
  GridEdger *ge = (GridEdger *)voidGridEdger;
  ge->grid = NULL;
}

int gridedgerEdgeId(GridEdger *ge)
{
  return ge->edgeId;
}

GridEdger *gridedgerDiscreteSegmentAndRatio(GridEdger *ge, double segment, 
					    int *discrete_segment, 
					    double *segment_ratio )
{
  int size;
  int segment_index;
  double ratio;
  Grid *grid = gridedgerGrid( ge );
  *discrete_segment = EMPTY;
  *segment_ratio = DBL_MAX;

  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );

  segment_index = (int)segment;
  ratio = segment - (double)segment_index;

  /* allow a little slop at the end points */
  if ( segment_index == -1 && ratio > (1.0-1.0e-12) ) {
    segment_index = 0;
    ratio = 0.0;
  }
  if ( segment_index == size-1 && ratio < 1.0e-12 ) {
    segment_index = size-2;
    ratio = 1.0;
  }

  if ( segment_index < 0 || segment_index > size-2 ) return NULL;

  *discrete_segment = segment_index;
  *segment_ratio = ratio;

  return ge;
}

GridEdger *gridedgerSegmentT(GridEdger *ge, double segment, double *t )
{
  int size;
  int segment_index;
  double ratio;
  int *curve;
  double t0, t1;
  Grid *grid = gridedgerGrid( ge );

  *t = DBL_MAX;

  if ( ge != gridedgerDiscreteSegmentAndRatio(ge, segment, 
					      &segment_index, 
					      &ratio ) ) return NULL;

  /* collect the edge segments into a curve */
  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  curve = malloc( size * sizeof(int) );
  if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
    free(curve);
    return NULL;
  }

  /* get the end points in t of the curve segment that we are in */
  gridNodeT(grid, curve[segment_index],   gridedgerEdgeId( ge ), &t0 );
  gridNodeT(grid, curve[segment_index+1], gridedgerEdgeId( ge ), &t1 );

  /* allowing a curve memory leak would suck */
  free(curve);

  /* linearally interpolate t in segment */
  *t = ratio*t1 + (1.0-ratio)*t0;

  /* make sure that this t value is supported by a sucessful cad evaluation */
  if ( !gridNewGeometryEdgeSiteAllowedAt( grid, 
					  curve[segment_index], 
					  curve[segment_index+1],
					  *t ) ) return NULL;
  return ge;
}

GridEdger *gridedgerSegmentMap( GridEdger *ge, double segment, double *map )
{
  int size;
  int segment_index;
  double ratio;
  int *curve;
  double map0[6], map1[6];
  Grid *grid = gridedgerGrid( ge );
  int i;

  for(i=0;i<6;i++) map[i] = DBL_MAX;

  if ( ge != gridedgerDiscreteSegmentAndRatio(ge, segment, 
					      &segment_index, 
					      &ratio ) ) return NULL;

  /* collect the edge segments into a curve */
  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  curve = malloc( size * sizeof(int) );
  if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
    free(curve);
    return NULL;
  }

  /* get the end points in t of the curve segment that we are in */
  gridMap(grid, curve[segment_index],   map0 );
  gridMap(grid, curve[segment_index+1], map1 );

  /* allowing a curve memory leak would suck */
  free(curve);

  /* linearally interpolate map in segment */
  for(i=0;i<6;i++) map[i] = ratio*map1[i] + (1.0-ratio)*map0[i];

  return ge;
}

GridEdger *gridedgerLengthToS(GridEdger *ge, double segment, double length,
			      double *next_s )
{
  *next_s = 1.0;
  
  /* bracket search to a single discrete edge segment
     by finding first segment end point that is too long */

  /* if the last segment end point for CAD curve is too short return it */

  /* do n-r to find the desired s */

  return ge;
}
