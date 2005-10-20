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
  int i;
  GridEdger *gm;
  gm = malloc(sizeof(GridEdger));
  gm->grid = grid;

  gm->edgeId = edgeId;

  gridAttachPacker( grid, gridedgerPack, (void *)gm );
  gridAttachNodeSorter( grid, gridedgerSortNode, (void *)gm );
  gridAttachReallocator( grid, gridedgerReallocator, (void *)gm );
  gridAttachFreeNotifier( grid, gridedgerGridHasBeenFreed, (void *)gm );

  return gm;
}

Grid *gridedgerGrid(GridEdger *gm)
{
  return gm->grid;
}

void gridedgerFree(GridEdger *gm)
{
  if (NULL != gm->grid) { 
    gridDetachPacker( gm->grid );
    gridDetachNodeSorter( gm->grid );
    gridDetachReallocator( gm->grid );
    gridDetachFreeNotifier( gm->grid );
  }
  free(gm);
}

void gridedgerPack(void *voidGridEdger, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerSortNode(void *voidGridEdger, int maxnode, int *o2n)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerReallocator(void *voidGridEdger, int reallocType, 
			 int lastSize, int newSize)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerGridHasBeenFreed(void *voidGridEdger )
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
  gm->grid = NULL;
}

int gridedgerEdgeId(GridEdger *gm)
{
  return gm->edgeId;
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

  /* collect the edge segments into a curve */
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
