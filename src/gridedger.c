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
  ge->nodes = 0;
  ge->s = NULL;

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
  if ( NULL != ge->s ) free( ge->s );
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

int gridedgerNodes(GridEdger *ge)
{
  return ge->nodes;
}

GridEdger *gridedgerNodeS(GridEdger *ge, int node, double *segment )
{
  if (NULL == ge->s) return NULL;
  if ( 0 > node || gridedgerNodes(ge) < node ) return NULL;
  *segment = ge->s[node];
  return ge;
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
  int node0, node1;
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

  /* to get segment end points */
  node0 = curve[segment_index];
  node1 = curve[segment_index+1];

  /* allowing a curve memory leak would suck */
  free(curve);

  /* get the end points in t of the curve segment that we are in */
  gridNodeT(grid, node0, gridedgerEdgeId( ge ), &t0 );
  gridNodeT(grid, node1, gridedgerEdgeId( ge ), &t1 );

  /* linearally interpolate t in segment */
  *t = ratio*t1 + (1.0-ratio)*t0;

  /* make sure that this t value is supported by a sucessful cad evaluation */
  if ( !gridNewGeometryEdgeSiteAllowedAt( grid, 
					  node0, 
					  node1,
					  *t ) ) return NULL;
  return ge;
}

GridEdger *gridedgerSegmentMap( GridEdger *ge, double segment, double *map )
{
  int size;
  int segment_index;
  double ratio;
  int *curve;
  int node0, node1;
  double map0[6], map1[6];
  int i;

  Grid *grid = gridedgerGrid( ge );

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

  /* to get segment end points */
  node0 = curve[segment_index];
  node1 = curve[segment_index+1];

  /* allowing a curve memory leak would suck */
  free(curve);

  /* get the end points in t of the curve segment that we are in */
  gridMap(grid, node0, map0 );
  gridMap(grid, node1, map1 );

  /* linearally interpolate map in segment */
  for(i=0;i<6;i++) map[i] = ratio*map1[i] + (1.0-ratio)*map0[i];

  return ge;
}

GridEdger *gridedgerLengthBetween(GridEdger *ge, 
				  double segment0, double segment1, 
				  double *length )
{
  double map0[6], map1[6], map[6];
  double t0, t1;
  double xyz0[3], xyz1[3];
  double dx, dy, dz;

  Grid *grid = gridedgerGrid( ge );

  *length = -1.0;

  if ( ge != gridedgerSegmentT( ge, segment0, &t0 ) ) return NULL;
  if ( ge != gridedgerSegmentT( ge, segment1, &t1 ) ) return NULL;

  if ( grid != gridEvaluateOnEdge( grid, gridedgerEdgeId( ge ), t0, xyz0 ) ) 
    return NULL;
  if ( grid != gridEvaluateOnEdge( grid, gridedgerEdgeId( ge ), t1, xyz1 ) ) 
    return NULL;

  dx = xyz1[0] - xyz0[0];
  dy = xyz1[1] - xyz0[1];
  dz = xyz1[2] - xyz0[2];

  if ( ge != gridedgerSegmentMap( ge, segment0, map0 ) ) return NULL;
  if ( ge != gridedgerSegmentMap( ge, segment1, map1 ) ) return NULL;
  
  map[0] = 0.5*(map0[0]+map1[0]);
  map[1] = 0.5*(map0[1]+map1[1]);
  map[2] = 0.5*(map0[2]+map1[2]);
  map[3] = 0.5*(map0[3]+map1[3]);
  map[4] = 0.5*(map0[4]+map1[4]);
  map[5] = 0.5*(map0[5]+map1[5]);

  *length =  sqrt ( dx * ( map[0]*dx + map[1]*dy + map[2]*dz ) +
		    dy * ( map[1]*dx + map[3]*dy + map[4]*dz ) +
		    dz * ( map[2]*dx + map[4]*dy + map[5]*dz ) );
  return ge;
}

GridEdger *gridedgerLengthToS(GridEdger *ge, double segment, double length,
			      double *next_s )
{
  int size;
  int segment_index;
  double ratio;
  int in_this_segment;
  int seg;
  double length_to_end_of_segment;
  double max, min, mid;
  int iteration;
  double mid_length;

  Grid *grid = gridedgerGrid( ge );

  *next_s = -1.0;
  
  if ( ge != gridedgerDiscreteSegmentAndRatio(ge, segment, 
					      &segment_index, 
					      &ratio ) ) return NULL;

  /* bracket search to a single discrete edge segment
     by finding first segment end point that is too long */

  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );

  in_this_segment = EMPTY;
  for ( seg = segment_index; seg < size-1; seg++ ) {
    if ( ge != gridedgerLengthBetween( ge, segment, 1.0 + (double)seg,
				       &length_to_end_of_segment ) ) 
      return NULL;
    if ( length_to_end_of_segment >= length ) {
      in_this_segment = seg; 
      break;
    }
  }

  /* if the last segment end point for CAD curve is too short
     because we could not find a segment the next point is in return it */
  if ( EMPTY == in_this_segment ) {
    *next_s = (double)(size-1);
    return ge;
  }

  /* do binary search to find the desired s */
  /* n-r would be better (quadratic convergence */

  /* 30 iterations gives 0.5^30 = 1e-10 convergence */
  max = (double)(in_this_segment+1);
  min = segment;
  mid = 0.5*(min+max);
  for ( iteration = 0 ; iteration < 40 ; iteration++ ) {
    if ( ge != gridedgerLengthBetween( ge, segment, mid, 
				       &mid_length ) ) return NULL;
    if ( mid_length > length ) {
      max = mid;
    }else{
      min = mid;
    }
    mid = 0.5*(min+max);
  }
  *next_s = mid;
  return ge;
}

GridEdger *gridedgerDiscretize(GridEdger *ge, double length )
{
  int max_size, chunk_size;
  int node;
  double endpoint;

  int size;

  Grid *grid = gridedgerGrid( ge );

  ge->nodes = 0;
  if ( NULL != ge->s ) {
    free( ge->s );
    ge->s = NULL;
  }
  chunk_size = 1000;
  max_size = 5*chunk_size;

  ge->s = (double *)malloc( max_size * sizeof(double) );

  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  endpoint = (double)(size-1);

  ge->s[0] = 0.0;
  node = 0;
  while ( ge->s[node] < (endpoint-1.0e-12) ) {
    node++;
    if (node >= max_size) {
      max_size += chunk_size;
      ge->s = (double *)realloc( ge->s, max_size * sizeof(double) );
    }
    if (ge != gridedgerLengthToS(ge, ge->s[node-1], length, &(ge->s[node]) )){
      free( ge->s );
      ge->s = NULL;
    }
  }

  ge->nodes = node+1;

  ge->s = (double *)realloc( ge->s, ge->nodes * sizeof(double) );  

  return ge;
}

GridEdger *gridedgerDiscretizeEvenly(GridEdger *ge )
{
  double length;
  int iteration;

  double s0, s1;
  double last_length;

  int size;

  length = 1.0;

  for (iteration = 1; iteration <= 5 ; iteration++) {
    if (ge != gridedgerDiscretize( ge, length ) ) return NULL;

    size = gridedgerNodes( ge );

    gridedgerNodeS( ge, size-2, &s0 );
    gridedgerNodeS( ge, size-1, &s1 );
    gridedgerLengthBetween( ge, s0, s1, &last_length );
    
    length = length * (double)(size-1) / 
      ((double)(size-1)+(length-last_length));

  }

  printf("edge%4d nodes%4d len %f %f\n",
	 gridedgerEdgeId( ge ), size, length, last_length);

  return ge;
}
