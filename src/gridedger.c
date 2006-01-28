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

#include "queue.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "gridedger.h"

GridEdger *gridedgerCreate( Grid *grid, int edgeId )
{
  GridEdger *ge;
  ge = malloc(sizeof(GridEdger));
  ge->grid = grid;

  ge->edgeId = edgeId;
  ge->ideal_nodes = 0;
  ge->t = NULL;

  ge->total_unused = 0;
  ge->unused = NULL;

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
  if ( NULL != ge->t ) free( ge->t );
  if ( NULL != ge->unused ) free( ge->unused );
  free(ge);
}

void gridedgerPack(void *voidGridEdger, 
		   int nnode, int maxnode, int *nodeo2n,
		   int ncell, int maxcell, int *cello2n,
		   int nface, int maxface, int *faceo2n,
		   int nedge, int maxedge, int *edgeo2n)
{
  int i;
  GridEdger *ge = (GridEdger *)voidGridEdger;
  for ( i = 0 ; i < gridedgerUnusedNodes(ge) ; i++ ) {
    ge->unused[i] = nodeo2n[ge->unused[i]];
  }
}

void gridedgerSortNode(void *voidGridEdger, int maxnode, int *o2n)
{
  int i;
  GridEdger *ge = (GridEdger *)voidGridEdger;
  for ( i = 0 ; i < gridedgerUnusedNodes(ge) ; i++ ) {
    ge->unused[i] = o2n[ge->unused[i]];
  }
}

void gridedgerReallocator(void *voidGridEdger, int reallocType, 
			  int lastSize, int newSize)
{

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

int gridedgerIdealNodes(GridEdger *ge)
{
  return ge->ideal_nodes;
}

GridEdger *gridedgerIdealNodeT(GridEdger *ge, int node, double *t )
{
  if (NULL == ge->t) return NULL;
  if ( 0 > node || gridedgerIdealNodes(ge) <= node ) return NULL;
  *t = ge->t[node];
  return ge;
}

int gridedgerUnusedNodes(GridEdger *ge)
{
  return ge->total_unused;
}

int gridedgerUnusedNode(GridEdger *ge, int index )
{
  if (NULL == ge->unused) return EMPTY;
  if ( 0 > index || gridedgerUnusedNodes(ge) <= index ) return EMPTY;
  return ge->unused[index];
}

GridEdger *gridedgerSupportingSegment(GridEdger *ge, 
				      double t, double *segment )
{
  int size;
  int *curve;
  int node0, node1;
  double t0, t1;
  int in_this_segment, segment_index;
  double ratio;

  Grid *grid = gridedgerGrid( ge );

  *segment = DBL_MAX;

  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  curve = malloc( size * sizeof(int) );
  if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
    free(curve);
    return NULL;
  }

  node0 = curve[0];
  gridNodeT(grid, node0, gridedgerEdgeId( ge ), &t0 );
  if ( t - t0 < -1.0e-12 ) {
    free(curve);
    return NULL;
  }
  
  node1 = curve[size-1];
  gridNodeT(grid, node1, gridedgerEdgeId( ge ), &t1 );
  if ( t - t1 > 1.0e-12 ) {
    free(curve);
    return NULL;
  }

  in_this_segment = EMPTY;
  for ( segment_index = 1; segment_index < size; segment_index++ ) {
    node0 = curve[segment_index];
    gridNodeT(grid, node0, gridedgerEdgeId( ge ), &t0 );

    if ( t - t0 < 1.0e-12 ) {
      in_this_segment = segment_index-1; 
      break;
    }
  }
  if (EMPTY == in_this_segment) {
    free(curve);
    return NULL;
  }

  node0 = curve[in_this_segment];
  gridNodeT(grid, node0, gridedgerEdgeId( ge ), &t0 );
  node1 = curve[in_this_segment+1];
  gridNodeT(grid, node1, gridedgerEdgeId( ge ), &t1 );
  /* t = ratio*t1 + (1.0-ratio)*t0 */
  /* t = ratio*t1 + t0 - ratio*t0 */
  /* t-t0 = ratio*(t1-t0) */
  ratio = (t-t0)/(t1-t0);

  *segment = (double)in_this_segment + ratio;

  free(curve);
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
  /* n-r would be better (quadratic convergence instead of linear) */

  /* 30 iterations gives 0.5^30 = 1e-10 convergence */
  max = (double)(in_this_segment+1);
  min = MAX(segment,((double)in_this_segment));
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
  double s0, s1;
  int size;

  Grid *grid = gridedgerGrid( ge );

  ge->ideal_nodes = 0;
  if ( NULL != ge->t ) {
    free( ge->t );
    ge->t = NULL;
  }
  chunk_size = 1000;
  max_size = 5*chunk_size;

  ge->t = (double *)malloc( max_size * sizeof(double) );

  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  endpoint = (double)(size-1);

  gridedgerSegmentT( ge, 0.0, &(ge->t[0]));
  s0 = 0.0;
  node = 0;
  while ( s0 < (endpoint-1.0e-12) ) {
    node++;
    if (node >= max_size) {
      max_size += chunk_size;
      ge->t = (double *)realloc( ge->t, max_size * sizeof(double) );
    }
    if (ge != gridedgerLengthToS(ge, s0, length, &s1 )){
      free( ge->t );
      ge->t = NULL;
      printf( "%s: %d: gridedgerDiscretize: gridedgerLengthToS NULL.\n",
	      __FILE__, __LINE__ );
      return NULL;
    }
    gridedgerSegmentT( ge, s1, &(ge->t[node]));
    s0 = s1;
  }

  ge->ideal_nodes = node+1;

  ge->t = (double *)realloc( ge->t, ge->ideal_nodes * sizeof(double) );  

  return ge;
}

static GridEdger *gridedgerReport(GridEdger *ge, double ratio )
{
  double longest_ratio, shortest_ratio;
  int node;
  double length;
  double t0, s0, t1, s1;

  longest_ratio = -1.0;
  shortest_ratio = 999.0;
  for ( node = 1 ; node < gridedgerIdealNodes( ge ) ; node++ ) {
    gridedgerIdealNodeT( ge, node-1, &t0 );
    gridedgerIdealNodeT( ge, node, &t1 );
    gridedgerSupportingSegment(ge, t0, &s0 );
    gridedgerSupportingSegment(ge, t1, &s1 );
    gridedgerLengthBetween( ge, s0, s1, &length );
    longest_ratio = MAX( longest_ratio, length );
    shortest_ratio = MIN( shortest_ratio, length );
  }
  printf("edge%4d longest %f shortest %f desired %f\n",
	 gridedgerEdgeId( ge ), longest_ratio, shortest_ratio, ratio);
  return ge;
}

GridEdger *gridedgerDiscretizeEvenly(GridEdger *ge )
{
  double length;
  int iteration;

  double t0, s0, s1;
  double last_length;
  double w, next_length;

  int size;
  int last_size;

  Grid *grid = gridedgerGrid( ge );

  length = 1.0;

  s1 = (double)(gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) )-1);

  w = 0.5;
  next_length = 0.9;
  for (iteration = 1; iteration <= 20 ; iteration++) {
    if (ge != gridedgerDiscretize( ge, length ) ) {
      printf( "%s: %d: %s: on iter %d for len %f try another length\n",
	      __FILE__, __LINE__, "gridedgerDiscretizeEvenly",
	      iteration, length );
      length = w*next_length + (1.0-w)*length;
      continue;
    }

    if (1==iteration) size = gridedgerIdealNodes( ge );
    last_size = gridedgerIdealNodes( ge );

    gridedgerIdealNodeT( ge, size-2, &t0 );
    gridedgerSupportingSegment(ge, t0, &s0 );

    gridedgerLengthBetween( ge, s0, s1, &last_length );
    
    /*
    printf("edge%4d nodes%4d len %f %f\n",
	   gridedgerEdgeId( ge ), last_size, length, last_length);
    */

    if (ABS(last_length-length) < 0.01) break;

    next_length = (((double)(last_size-2))*length+last_length) / 
                   ((double)(size-1));

    length = w*next_length + (1.0-w)*length;

  }

  gridedgerReport(ge, length);

  return ge;
}

GridEdger *gridedgerDiscretizeOnce(GridEdger *ge )
{
  double length;
  int size;
  double t0, s0, t1, s1;
  double last_length;
  double ratio, shrink;

  double *orig;
  int node;

  Grid *grid = gridedgerGrid( ge );

  length = 1.0;
  if (ge != gridedgerDiscretize( ge, length ) ) return NULL;

  size = gridedgerIdealNodes( ge );

  gridedgerIdealNodeT( ge, size-2, &t0 );
  gridedgerSupportingSegment(ge, t0, &s0 );
  s1 = (double)(gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) )-1);
  gridedgerLengthBetween( ge, s0, s1, &last_length );
    
  ratio = (((double)(size-2))*length+last_length) / ((double)(size-1));

  orig = (double *)malloc( gridedgerIdealNodes( ge ) * sizeof(double) );  
  for ( node = 0 ; node < gridedgerIdealNodes( ge ) ; node++ )
    gridedgerIdealNodeT( ge, node, &(orig[node]) );
  
  shrink = 1.0-ratio;
  for ( node = 1 ; node < (gridedgerIdealNodes( ge )-1) ; node++ ) {
    t0 = orig[node-1];
    t1 = orig[node];
    ratio = 1.0-(shrink*((double)node));
    ge->t[node] = (1.0-ratio)*t0 + ratio*t1;
  }
  free(orig);

  gridedgerReport(ge, 1.0-shrink);

  return ge;
}

GridEdger *gridedgerDiscretizeSupport(GridEdger *ge, int subintervals )
{
  int segments, subsegments;
  double *s, *l;
  int parent;
  int node;
  double s0, s1, delta_length;
  int target_size;
  double target_delta_length;
  double target_length;
  double target_s;

  Grid *grid = gridedgerGrid( ge );

  ge->ideal_nodes = 0;
  if ( NULL != ge->t ) {
    free( ge->t );
    ge->t = NULL;
  }

  segments = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) ) - 1;

  subsegments = segments * subintervals;

  s = (double *)malloc( (subsegments+1) * sizeof(double) );  
  l = (double *)malloc( (subsegments+1) * sizeof(double) );

  for ( node = 0 ; node <= subsegments ; node++ ) {
    parent = node/subintervals;
    s[node] = ((double)parent) + 
      ((double)(node-parent*subintervals)) / ((double)subintervals);
  }
  l[0] = 0.0;
  for ( node = 1 ; node <= subsegments ; node++ ) {
    s0 = s[node-1];
    s1 = s[node];
    gridedgerLengthBetween( ge, s0, s1, &delta_length );
    l[node] = l[node-1] + delta_length;
  }

  target_size = ((int)l[subsegments])+1;
  target_delta_length = l[subsegments]/((double)target_size) ;

  ge->t = (double *)malloc( (target_size+1) * sizeof(double) );

  gridedgerSegmentT( ge,              0.0, &(ge->t[0]) );
  gridedgerSegmentT( ge, (double)segments, &(ge->t[target_size]) );

  s0 = 0.0;
  for ( node = 1 ; node < target_size ; node++ ) {
    target_length = target_delta_length * ((double)node);
    target_s = 0.0;
    // target_s = f(last_s,target_length);
    gridedgerSegmentT( ge, target_s, &(ge->t[node]) );
  }


  free(s);
  free(l);

  return ge;
}

GridEdger *gridedgerInsert(GridEdger *ge )
{
  int node;
  double t;
  double segment;
  int segment_index;
  double segment_ratio;
  int size;
  int *curve;
  int node0, node1;
  int newnode;
  int unused;
  GridBool found;

  Queue *queue = NULL;
  Grid *grid = gridedgerGrid( ge );

  /* save off existing edge interior nodes as unused to mark for removal */
  size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
  curve = malloc( size * sizeof(int) );
  if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
    free(curve);
    return NULL;
  }
  ge->total_unused = size-2;
  if (NULL==ge->unused) free(ge->unused);
  ge->unused = (int *)malloc( ge->total_unused * sizeof(int) );
  for ( node = 1 ; node < (size-1) ; node++ ) ge->unused[node-1] = curve[node];
  free(curve);  

  if ( 1 > gridedgerIdealNodes( ge ) ) return NULL;
  
  for ( node = 1 ; node < (gridedgerIdealNodes( ge )-1) ; node++ ) {

    /* the ideal location in t for this next node */
    if (ge != gridedgerIdealNodeT(ge, node, &t )) {
      ge->total_unused = 0; free(ge->unused); ge->unused=NULL;
      return NULL;
    }
    /* find the segment that this t value lies inside  */
    if (ge != gridedgerSupportingSegment(ge, t, &segment )) {
      ge->total_unused = 0; free(ge->unused); ge->unused=NULL;
      return NULL;
    }
    if (ge != gridedgerDiscreteSegmentAndRatio(ge, segment, 
					       &segment_index, 
					       &segment_ratio ) ) {
      ge->total_unused = 0; free(ge->unused); ge->unused=NULL;
      return NULL;
    }

    /* collect the edge segments into a curve */
    size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
    curve = malloc( size * sizeof(int) );
    if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
      ge->total_unused = 0; free(ge->unused); ge->unused=NULL;
      free(curve);
      return NULL;
    }

    /* to get segment end points */
    node0 = curve[segment_index];
    node1 = curve[segment_index+1];

    if ( segment_ratio < 0.000001 || segment_ratio > 0.999999 ) {
      /* reuse existing node at the ideal location */
      newnode = node0;
      if (segment_ratio > 0.5) newnode = node1;

      /* extract newnode from unused list */
      found = FALSE;
      for ( unused = 0 ; unused < ge->total_unused ; unused++ ) {
	found = found || ( newnode == ge->unused[unused]);
	if (found && unused < ge->total_unused-1) 
	  ge->unused[unused] = ge->unused[unused+1];
      }
      
      if (found) {
	ge->total_unused--;
      }else{
	printf( "%s: %d: gridedgerInsert: unable to remove used from unused.\n",
		__FILE__, __LINE__ );
	newnode = EMPTY;
      }
    }else{
      /* actually insert the new node in the ideal location */
      newnode = gridSplitEdgeRatio( grid, queue, node0, node1, segment_ratio );
    }

    /* allowing a curve memory leak would suck */
    free(curve);

    /* return NULL if the new node was not added */
    if (EMPTY == newnode) {
      ge->total_unused = 0; free(ge->unused); ge->unused=NULL;
      printf( "ERROR: gridedgerInsert: %s: %d: new node not added\n",
	      __FILE__, __LINE__ );
      return NULL;
    }

  }

  return ge;
}

GridEdger *gridedgerRemoveUnused(GridEdger *ge )
{
  int unused;
  int nodeA, nodeB, node1;
  double current_cost, cost1, costA, costB;
  int size;
  int *curve;
  int i, target;

  GridBool removed;

  Queue *queue = NULL;
  Grid *grid = gridedgerGrid( ge );

  unused = 0;
  while ( unused < ge->total_unused ) {
    node1 = ge->unused[unused];
    
    /* I should look at connected edges instead of curve for efficiency */

    /* collect the edge segments into a curve */
    size = gridGeomEdgeSize( grid, gridedgerEdgeId( ge ) );
    curve = malloc( size * sizeof(int) );
    if (grid != gridGeomEdge( grid, gridedgerEdgeId( ge ), curve )) {
      free(curve);
      return NULL;
    }
    target = EMPTY;
    for ( i = 1 ; i < (size-1) ; i++) if (node1 == curve[i]) target = i;
    if ( EMPTY == target ) {
      printf( "ERROR: gridedgerRemoveUnused: %s: %d: target not found\n",
	      __FILE__, __LINE__ );
      free(curve);
      return NULL;
    }
    nodeA = curve[target-1];
    nodeB = curve[target+1];
    free(curve);

    /* collapse to another node on curve */
    if (grid!=gridCollapseCost(grid, nodeA, node1, 
			       &current_cost, &costA, &cost1)) {
      printf( "ERROR: gridedgerRemoveUnused: %s: %d: gridCollapseCost NULL\n",
	      __FILE__, __LINE__ );
      return NULL;
    }
    if (grid!=gridCollapseCost(grid, nodeB, node1, 
			       &current_cost, &costB, &cost1)) {
      printf( "ERROR: gridedgerRemoveUnused: %s: %d: gridCollapseCost NULL\n",
	      __FILE__, __LINE__ );
      return NULL;
    }

    removed = FALSE;
    if (costA>costB) {
      removed = ( grid==gridCollapseEdge(grid, queue, nodeA, node1, 0.0 ) );
      if (!removed) { 
	removed = ( grid==gridCollapseEdge(grid, queue, nodeB, node1, 0.0 ) );
      }
    }else{
      removed = ( grid==gridCollapseEdge(grid, queue, nodeB, node1, 0.0 ) );
      if (!removed) { 
	removed = ( grid==gridCollapseEdge(grid, queue, nodeA, node1, 0.0 ) );
      }
    }

    if (removed) { 
      ge->total_unused--;
      for ( i = unused ; i < ge->total_unused ; i++) 
	ge->unused[i] = ge->unused[i+1];
    }else{
      unused++;
    }
  }

  if ( ge->total_unused == 0 ) {
    free(ge->unused); ge->unused=NULL;
  }else{
    ge->unused = (int *)realloc( ge->unused, ge->total_unused * sizeof(int) );
    printf("edge%4d unused remaining %d\n",
	   gridedgerEdgeId( ge ), ge->total_unused );
  }

  return ge;
}
