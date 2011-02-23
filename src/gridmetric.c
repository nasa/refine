
/* Computes metrics from faces and tets 
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifdef __APPLE__       /* Not needed on Mac OS X */
#include <float.h>
#else
#include <values.h>
#endif
#ifdef HAVE_SDK
#include "CADGeom/CADGeom.h"
#else
#include "FAKEGeom.h"
#endif
#include "gridshape.h"
#include "gridmetric.h"

#define GRID_VOLUME_TOL (1.0e-40)
#define GRID_AREA_TOL (1.0e-40)

Grid *gridWriteTecplotInvalid(Grid *grid, char *filename )
{
  char comment[256];
  int cell, nodes[4];
  double cost, costs[4];
  double min_cost;
  int faceid;

  gridWriteTecplotSurfaceGeom(grid,filename);

  min_cost = 999.0;
  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      cost = gridAR(grid,nodes);
      min_cost = MIN(min_cost,cost);
      if ( -0.5 > cost ) {
	sprintf(comment,"cell cost of %f detected.",cost);
	gridWriteTecplotComment(grid, comment);
	costs[0] = costs[1] = costs[2] = costs[3] = cost;
	gridWriteTecplotCellGeom(grid,nodes,costs,filename);
      }
    }
  }

  if (min_cost < -1.5) {
    for (faceid=1;faceid<=gridNGeomFace(grid);faceid++) {
      if (0 < gridTrianglesOnFaceId(grid, faceid )) {
	gridWriteTecplotGeomFaceUV(grid,filename,faceid);
      }
    }
  }
  
  return grid;
}

Grid *gridSetMapWithSpacingVectors(Grid *grid, int node,
				   double *v1, double *v2, double *v3,
                                   double s1, double s2, double s3)
{
  double m11, m12, m13, m22, m23, m33;
  double e1, e2, e3;

  e1 = 1.0/s1/s1;
  e2 = 1.0/s2/s2;
  e3 = 1.0/s3/s3;

  m11 = v1[0]*v1[0]*e1 + v2[0]*v2[0]*e2 + v3[0]*v3[0]*e3;
  m12 = v1[0]*v1[1]*e1 + v2[0]*v2[1]*e2 + v3[0]*v3[1]*e3;
  m13 = v1[0]*v1[2]*e1 + v2[0]*v2[2]*e2 + v3[0]*v3[2]*e3;
  m22 = v1[1]*v1[1]*e1 + v2[1]*v2[1]*e2 + v3[1]*v3[1]*e3;
  m23 = v1[1]*v1[2]*e1 + v2[1]*v2[2]*e2 + v3[1]*v3[2]*e3;
  m33 = v1[2]*v1[2]*e1 + v2[2]*v2[2]*e2 + v3[2]*v3[2]*e3;

  return gridSetMap(grid,node,m11,m12,m13,m22,m23,m33);
}

Grid *gridSetMapMatrixToAverageOfNodes2(Grid *grid, int avgNode,
					int n0, int n1 )
{
  return gridInterpolateMap2(grid, n0, n1, 0.5, avgNode);
}

Grid *gridSetMapMatrixToAverageOfNodes3(Grid *grid, int avgNode,
					int n0, int n1, int n2 )
{
  int i;
  double map[6];
  double *m0, *m1, *m2;
  double map0[6], map1[6], map2[6];

  if ( !gridValidNode(grid,n0) ||
       !gridValidNode(grid,n1) ||
       !gridValidNode(grid,n2) ) return NULL;

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid, n0);
    m1 = gridMapPointer(grid, n1);
    m2 = gridMapPointer(grid, n2);
    for (i=0;i<6;i++) map[i] = (1.0/3.0)*(m0[i]+m1[i]+m2[i]);
  }else{
    gridMap(grid, n0, map0);
    gridMap(grid, n1, map1);
    gridMap(grid, n2, map2);
    for (i=0;i<6;i++) map[i] = (1.0/3.0)*(map0[i]+map1[i]+map2[i]);
  }

  if (grid != gridSetMap( grid, avgNode, 
			  map[0], map[1], map[2], 
			  map[3], map[4], map[5] ) ) return NULL;

  return grid;
}

Grid *gridSetMapMatrixToAverageOfNodes4(Grid *grid, int avgNode,
					int n0, int n1, int n2, int n3 )
{
  int i;
  double map[6];
  double *m0, *m1, *m2, *m3;
  double map0[6], map1[6], map2[6], map3[6];

  if ( !gridValidNode(grid,n0) ||
       !gridValidNode(grid,n1) ||
       !gridValidNode(grid,n2) ||
       !gridValidNode(grid,n3) ) return NULL;

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid, n0);
    m1 = gridMapPointer(grid, n1);
    m2 = gridMapPointer(grid, n2);
    m3 = gridMapPointer(grid, n3);
    for (i=0;i<6;i++) map[i] = (1.0/4.0)*(m0[i]+m1[i]+m2[i]+m3[i]);
  }else{
    gridMap(grid, n0, map0);
    gridMap(grid, n1, map1);
    gridMap(grid, n2, map2);
    gridMap(grid, n3, map3);
    for (i=0;i<6;i++) map[i] = (1.0/4.0)*(map0[i]+map1[i]+map2[i]+map3[i]);
  }

  if (grid != gridSetMap( grid, avgNode, 
			  map[0], map[1], map[2], 
			  map[3], map[4], map[5] ) ) return NULL;

  return grid;
}

void gridMapXYZWithJ( double *j,
		      double *x, double *y, double *z )
{
  double mapx, mapy, mapz;
  
  mapx = j[0] * *x + j[1] * *y + j[2] * *z; 
  mapy = j[3] * *x + j[4] * *y + j[5] * *z; 
  mapz = j[6] * *x + j[7] * *y + j[8] * *z; 

  *x = mapx;
  *y = mapy;
  *z = mapz;
}

double gridEdgeLength(Grid *grid, int n0, int n1 )
{
  double xyz0[3], xyz1[3];
  double dx, dy, dz;

  if (grid != gridNodeXYZ(grid, n0, xyz0) ) return -1.0;
  if (grid != gridNodeXYZ(grid, n1, xyz1) ) return -1.0;
  
  dx = xyz1[0] - xyz0[0];
  dy = xyz1[1] - xyz0[1];
  dz = xyz1[2] - xyz0[2];

  return  sqrt(dx*dx+dy*dy+dz*dz);
}

Grid *gridEdgeRatioTolerence(Grid *grid, double longest, double shortest,
			     int *active_edges, int *out_of_tolerence_edges )
{
  int conn, nodes[2];
  double ratio;

  gridCreateConn(grid);

  *active_edges = 0;
  *out_of_tolerence_edges = 0;

  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    if ( !gridNodeFrozen(grid, nodes[0]) &&
	 !gridNodeFrozen(grid, nodes[0]) )
      {
	(*active_edges) += 1;
	ratio = gridEdgeRatio(grid, nodes[0], nodes[1]);
	if ( ratio > longest || shortest > ratio )
	  (*out_of_tolerence_edges) += 1;
      }
  }

  gridEraseConn(grid);

  return grid;
}

Grid *gridEdgeRatioRange(Grid *grid, double *longest, double *shortest )
{
  int cell, nodes[4];
  double ratio;

  *shortest = DBL_MAX;
  *longest = -DBL_MAX;

  for( cell=0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      ratio = gridEdgeRatio(grid, nodes[0], nodes[1]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      ratio = gridEdgeRatio(grid, nodes[0], nodes[2]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      ratio = gridEdgeRatio(grid, nodes[0], nodes[3]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      ratio = gridEdgeRatio(grid, nodes[1], nodes[2]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      ratio = gridEdgeRatio(grid, nodes[1], nodes[3]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      ratio = gridEdgeRatio(grid, nodes[2], nodes[3]);
      *longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
    }
  }
  return grid;
}

Grid *gridEdgeRatioRangeInVolume(Grid *grid, double *longest, double *shortest )
{
  int cell, nodes[4];
  double ratio;

  *shortest = DBL_MAX;
  *longest = -DBL_MAX;

  for( cell=0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      if ( 0 == gridParentGeometry(grid, nodes[0], nodes[1]) ) {
	ratio = gridEdgeRatio(grid, nodes[0], nodes[1]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
      if ( 0 == gridParentGeometry(grid, nodes[0], nodes[2]) ) {
	ratio = gridEdgeRatio(grid, nodes[0], nodes[2]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
      if ( 0 == gridParentGeometry(grid, nodes[0], nodes[3]) ) {
	ratio = gridEdgeRatio(grid, nodes[0], nodes[3]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
      if ( 0 == gridParentGeometry(grid, nodes[1], nodes[2]) ) {
	ratio = gridEdgeRatio(grid, nodes[1], nodes[2]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
      if ( 0 == gridParentGeometry(grid, nodes[1], nodes[3]) ) {
	ratio = gridEdgeRatio(grid, nodes[1], nodes[3]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
      if ( 0 == gridParentGeometry(grid, nodes[2], nodes[3]) ) {
	ratio = gridEdgeRatio(grid, nodes[2], nodes[3]);
	*longest = MAX( *longest, ratio);*shortest = MIN( *shortest, ratio);
      }
    }
  }
  return grid;
}

double gridEdgeRatio(Grid *grid, int n0, int n1 )
{
  int i;
  double *xyz0, *xyz1;
  double dx, dy, dz;
  double *m0, *m1;
  double map0[6], map1[6];
  double m[6];

  if (!gridValidNode(grid, n0) || !gridValidNode(grid, n1)) return -1.0;

  xyz0 = gridNodeXYZPointer(grid, n0);
  xyz1 = gridNodeXYZPointer(grid, n1);
  
  dx = xyz1[0] - xyz0[0];
  dy = xyz1[1] - xyz0[1];
  dz = xyz1[2] - xyz0[2];

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid, n0);
    m1 = gridMapPointer(grid, n1);

    for (i=0;i<6;i++) m[i] = 0.5*(m0[i]+m1[i]);
  }else{
    if ( ( grid != gridMap(grid, n0, map0) ) ||
	 ( grid != gridMap(grid, n1, map1) ) ) {
      printf("%s: %d: gridEdgeRatio: gridMap NULL\n",__FILE__,__LINE__);
      return -1.0;
    }

    for (i=0;i<6;i++) m[i] = 0.5*(map0[i]+map1[i]);
  }

  return sqrt ( dx * ( m[0]*dx + m[1]*dy + m[2]*dz ) +
		dy * ( m[1]*dx + m[3]*dy + m[4]*dz ) +
	        dz * ( m[2]*dx + m[4]*dy + m[5]*dz ) );
}

Grid *gridEdgeRatio3(Grid *grid, int n0, int n1, double *ratio )
{
  double *xyz0, *xyz1;
  double dx, dy, dz;
  double *m0, *m1;
  double map0[6], map1[6];
  double m[6];

  if (!gridValidNode(grid, n0) || !gridValidNode(grid, n1)) return NULL;

  xyz0 = gridNodeXYZPointer(grid, n0);
  xyz1 = gridNodeXYZPointer(grid, n1);
  
  dx = xyz1[0] - xyz0[0];
  dy = xyz1[1] - xyz0[1];
  dz = xyz1[2] - xyz0[2];

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid, n0);
    m1 = gridMapPointer(grid, n1);
  }else{
    gridMap(grid,n0,map0); m0 = map0;
    gridMap(grid,n1,map1); m1 = map1;
  }
  
  ratio[0] = sqrt (
      dx * ( m0[0]*dx + m0[1]*dy + m0[2]*dz )
    + dy * ( m0[1]*dx + m0[3]*dy + m0[4]*dz )
    + dz * ( m0[2]*dx + m0[4]*dy + m0[5]*dz ) );

  ratio[1] = sqrt (
      dx * ( m1[0]*dx + m1[1]*dy + m1[2]*dz )
    + dy * ( m1[1]*dx + m1[3]*dy + m1[4]*dz )
    + dz * ( m1[2]*dx + m1[4]*dy + m1[5]*dz ) );
  
  m[0] = 0.5*(m0[0]+m1[0]);
  m[1] = 0.5*(m0[1]+m1[1]);
  m[2] = 0.5*(m0[2]+m1[2]);
  m[3] = 0.5*(m0[3]+m1[3]);
  m[4] = 0.5*(m0[4]+m1[4]);
  m[5] = 0.5*(m0[5]+m1[5]);

  ratio[2] = sqrt (
      dx * ( m[0]*dx + m[1]*dy + m[2]*dz )
    + dy * ( m[1]*dx + m[3]*dy + m[4]*dz )
    + dz * ( m[2]*dx + m[4]*dy + m[5]*dz ) );
  
  return grid;
}

double gridEdgeRatioError(Grid *grid, int n0, int n1 )
{
  double ratio;
  ratio = gridEdgeRatio(grid, n0, n1 );
  if (ratio<0.0) return ratio;
  return ABS((1.0-ratio)/(1.0+ratio));
}

double gridAverageEdgeLength(Grid *grid, int node )
{
  AdjIterator it;
  int ncell, cell, nodes[4];
  double length, celllength;
  ncell = 0;
  length = 0.0;
  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    ncell++;
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    celllength 
      = gridEdgeLength( grid, nodes[0], nodes[1] )
      + gridEdgeLength( grid, nodes[0], nodes[2] )
      + gridEdgeLength( grid, nodes[0], nodes[3] )
      + gridEdgeLength( grid, nodes[1], nodes[2] )
      + gridEdgeLength( grid, nodes[1], nodes[3] )
      + gridEdgeLength( grid, nodes[2], nodes[3] ) ;
    length += celllength/6.0;      
  }

  return length/(double)ncell;
}

Grid *gridLargestRatioEdge(Grid *grid, int node, 
			   int *edgeNode, double *ratio )
{
  AdjIterator it;
  int i, cell, nodes[4];
  double currentRatio;
  
  *edgeNode = EMPTY;
  *ratio = -1.0;
  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    for (i=0;i<4;i++){
      if (node != nodes[i]) {
	currentRatio = gridEdgeRatio( grid, node, nodes[i] );
	if ( currentRatio > *ratio ){
	  *ratio = currentRatio;
	  *edgeNode = nodes[i];
	}
      }
    }
  }
  
  return grid;
}

Grid *gridSmallestRatioEdge(Grid *grid, int node, 
			    int *edgeNode, double *ratio )
{
  AdjIterator it;
  int i, cell, nodes[4];
  double currentRatio;
  
  *edgeNode = EMPTY;
  *ratio = DBL_MAX;
  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    for (i=0;i<4;i++){
      if (node != nodes[i]) {
	currentRatio = gridEdgeRatio( grid, node, nodes[i] );
	if ( currentRatio < *ratio ){
	  *ratio = currentRatio;
	  *edgeNode = nodes[i];
	}
      }
    }
  }
  
  return grid;
}

double gridSpacing(Grid *grid, int node )
{
  double map[6];
  if (grid != gridMap(grid, node, map)) return -1.0;
  return 1.0/sqrt(map[0]);
}

Grid *gridSetSpacing(Grid *grid, int node, double spacing )
{
  double spacingInverse;
  spacingInverse = 1.0/spacing;
  spacingInverse = spacingInverse * spacingInverse;
  if ( grid != 
      gridSetMap(grid,node,
		 spacingInverse,0,0,
		 spacingInverse,0,
		 spacingInverse)) return(NULL);
  return grid;
}

Grid *gridResetSpacing(Grid *grid )
{
  int node;
  double spacingInverse;
  for ( node=0; node < gridMaxNode(grid); node++) {
    if ( gridValidNode(grid,node) ) {
      spacingInverse = 1.0/gridAverageEdgeLength( grid, node );
      spacingInverse = spacingInverse * spacingInverse;
      gridSetMap(grid,node,
		 spacingInverse,0,0,
		 spacingInverse,0,
		 spacingInverse);
    }
  } 
  return grid;
}

Grid *gridScaleSpacing(Grid *grid, int node, double scale )
{
  double map[6];
  double s;
  if (grid != gridMap(grid, node, map)) return NULL;
  s = 1.0/(scale*scale);
  gridSetMap(grid, node, 
	     map[0]*s, map[1]*s, map[2]*s, 
	     map[3]*s, map[4]*s, map[5]*s);
  return grid;
}

Grid *gridScaleSpacingSphere( Grid *grid, 
			      double x, double y, double z, double r,
			      double scale )
{
  int node;
  double xyz[3];
  double dx, dy, dz;
  double distanceSquared, radiusSquared;
  radiusSquared = r*r;
  
  for ( node=0; node<gridMaxNode(grid); node++ ) {
    if (grid == gridNodeXYZ(grid,node,xyz)) {
      dx = xyz[0] - x;
      dy = xyz[1] - y;
      dz = xyz[2] - z;
      distanceSquared = dx*dx + dy*dy + dz*dz;
      if (radiusSquared >= distanceSquared) gridScaleSpacing(grid,node,scale);
    }
  }

  return grid;
}

Grid *gridScaleSpacingSphereDirection( Grid *grid, 
			      double x, double y, double z, double r,
			      double scalex, double scaley, double scalez )
{
  int node;
  double xyz[3];
  double dx, dy, dz;
  double distanceSquared, radiusSquared;
  double map[6];
  radiusSquared = r*r;
  
  for ( node=0; node<gridMaxNode(grid); node++ ) {
    if (grid == gridNodeXYZ(grid,node,xyz)) {
      dx = xyz[0] - x;
      dy = xyz[1] - y;
      dz = xyz[2] - z;
      distanceSquared = dx*dx + dy*dy + dz*dz;
      if (radiusSquared >= distanceSquared){
	gridMap(grid, node, map);
	gridSetMap(grid, node, 
		   map[0] / (scalex*scalex), map[1], map[2], 
		   map[3] / (scaley*scaley), map[4], 
		   map[5] / (scalez*scalez));
      } 
    }
  }

  return grid;
}

Grid *gridSetGlobalMap(Grid *grid,
		       double m11, double m12, double m13,
		                   double m22, double m23,
		                               double m33)
{
  int node;
  for ( node=0; node<gridMaxNode(grid); node++ )
    if ( gridValidNode(grid,node) )
      gridSetMap(grid, node,
		 m11, m12, m13,
   		      m22, m23,
		           m33);
  return grid;
}

Grid *gridCopySpacing(Grid *grid, int originalNode, int newNode)
{
  double map[6];
  if (grid != gridMap(grid, originalNode, map)) return NULL;
  if (grid != gridSetMap( grid, newNode, 
			  map[0], map[1], map[2], 
			  map[3], map[4], map[5] ) ) return NULL;

  return grid;
}

Grid *gridConvertMetricToJacobian(Grid *grid, double *m, double *j)
{
  double d[3], e[3], v0[3], v1[3], v2[3];
  double e0, e1, e2;

  gridTriDiag3x3(m, d, e, v0, v1, v2);
  if ( !gridEigTriDiag3x3(d, e, v0, v1, v2 )) {
    printf("%s: %d: gridConvertMetricToJacobian: gridEigTriDiag3x3 FAILED.\n",
	   __FILE__,__LINE__);
    return NULL;
  }

  /* the new EigTriDiag should be ortho-normal */
  gridEigOrtho3x3( v0, v1, v2 );

  e0 = sqrt(d[0]);
  e1 = sqrt(d[1]);
  e2 = sqrt(d[2]);

  /* sqrt(eigenValues) * transpose(eigenVectors) */
  j[0] = e0*v0[0];
  j[1] = e0*v0[1];
  j[2] = e0*v0[2];
  j[3] = e1*v1[0];
  j[4] = e1*v1[1];
  j[5] = e1*v1[2];
  j[6] = e2*v2[0];
  j[7] = e2*v2[1];
  j[8] = e2*v2[2];

  return grid;
}

double gridVolume(Grid *grid, int *nodes )
{
  double *xyz0, *xyz1, *xyz2, *xyz3;
  double edge1[3], edge2[3], edge3[3];
  double norm[3];
  
  if ( !gridValidNode(grid, nodes[0]) || 
       !gridValidNode(grid, nodes[1]) ||
       !gridValidNode(grid, nodes[2]) ||
       !gridValidNode(grid, nodes[3]) ) return -1.0;

  xyz0=gridNodeXYZPointer(grid,nodes[0]);
  xyz1=gridNodeXYZPointer(grid,nodes[1]);
  xyz2=gridNodeXYZPointer(grid,nodes[2]);
  xyz3=gridNodeXYZPointer(grid,nodes[3]);

  gridSubtractVector( xyz1, xyz0, edge1);
  gridSubtractVector( xyz2, xyz0, edge2);
  gridSubtractVector( xyz3, xyz0, edge3);
  gridCrossProduct( edge1, edge2, norm );

  return gridDotProduct(norm,edge3)/6.0;
}

Grid *gridNodeCostValid(Grid *grid, int node, double *valid )
{
  AdjIterator it;
  int cell, nodes[4];
  double local_valid;

  *valid = 0.0;

  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    local_valid = gridCostValid(grid, nodes);
    if ( local_valid < *valid ) *valid = local_valid;
  }

  return grid;
}

Grid *gridNodeAR(Grid *grid, int node, double *ar )
{
  AdjIterator it;
  int cell, nodes[4];
  double local_ar;

  *ar = 1.0;

  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    local_ar = gridAR(grid, nodes);
    if ( local_ar < *ar ) *ar = local_ar;
  }

  return grid;
}

Grid *gridNodeVolume(Grid *grid, int node, double *volume )
{
  AdjIterator it;
  int cell, nodes[4];
  double local_volume;

  *volume = DBL_MAX;

  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    cell = adjItem(it);
    gridCell( grid, cell, nodes);
    local_volume = gridVolume(grid, nodes);
    if ( local_volume < *volume ) *volume = local_volume;
  }

  return grid;
}

Grid *gridGemAR( Grid *grid, double *ar ){
  int i, nodes[4];

  *ar = 2.0;

  for ( i = 0 ; i < gridNGem(grid) ; i++ ){
    gridCell(grid, gridGem(grid,i), nodes);
    *ar = MIN(*ar,gridAR( grid, nodes ));
  }

  return grid;
}

double gridCostValid(Grid *grid, int *nodes )
{
  int nodes_on_surface;
  if ( !gridValidNode(grid, nodes[0]) || 
       !gridValidNode(grid, nodes[1]) ||
       !gridValidNode(grid, nodes[2]) ||
       !gridValidNode(grid, nodes[3]) ) return -4.0;
  
  if ( (gridCostConstraint(grid)&gridCOST_CNST_AREAUV) ||
       (gridCostConstraint(grid)&gridCOST_CNST_VALID)  ) {
    nodes_on_surface = 0;
    if ( gridGeometryFace(grid, nodes[0]) ) nodes_on_surface++;
    if ( gridGeometryFace(grid, nodes[1]) ) nodes_on_surface++;
    if ( gridGeometryFace(grid, nodes[2]) ) nodes_on_surface++;
    if ( gridGeometryFace(grid, nodes[3]) ) nodes_on_surface++;
    if ( nodes_on_surface > 1 ) {
      if ( (gridCostConstraint(grid)&gridCOST_CNST_VALID) &&
	   ( gridMinCellJacDet2(grid,nodes) <= GRID_VOLUME_TOL ) ) return -3.0;
      if ( ( nodes_on_surface > 2 ) &&
	   (gridCostConstraint(grid)&gridCOST_CNST_AREAUV) &&
	   ( gridMinCellFaceAreaUV(grid,nodes) <= GRID_AREA_TOL ) ) return -2.0;
    }
  }

  if (gridCostConstraint(grid)&gridCOST_CNST_VOLUME) {
    if ( gridVolume(grid, nodes ) <= GRID_VOLUME_TOL ) return -1.0;
  }

  return 0.0;
}

double gridAR(Grid *grid, int *nodes )
{
  double xyz1[3], xyz2[3], xyz3[3], xyz4[3]; 
  double *p1, *p2, *p3, *p4; 
  int i;
  double *m0, *m1, *m2, *m3; 
  double map0[6], map1[6], map2[6], map3[6];
  double m[6], j[9];
  double aspect;

  double valid;

  valid = gridCostValid(grid, nodes );
  if ( -0.5 > valid ) return valid;
  
  if ( gridCOST_FCN_EDGE_LENGTH == gridCostFunction(grid) )
    return gridEdgeRatioCost(grid, nodes);

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid,nodes[0]);
    m1 = gridMapPointer(grid,nodes[1]);
    m2 = gridMapPointer(grid,nodes[2]);
    m3 = gridMapPointer(grid,nodes[3]);

    for (i=0;i<6;i++) m[i]=0.25*(m0[i]+m1[i]+m2[i]+m3[i]);
  } else {
    if ( (grid!=gridMap(grid, nodes[0], map0) ) ||
	 (grid!=gridMap(grid, nodes[1], map1) ) ||
	 (grid!=gridMap(grid, nodes[2], map2) ) ||
	 (grid!=gridMap(grid, nodes[3], map3) ) ) {
      return -999.0;
    }
    for (i=0;i<6;i++) m[i]=0.25*(map0[i]+map1[i]+map2[i]+map3[i]);
  }

  if (grid != gridConvertMetricToJacobian(grid, m, j) ) {
    printf("%s: %d: gridAR: gridConvertMetricToJacobian NULL\n",
	   __FILE__,__LINE__);
    printf("nodes %d %d %d %d\n",nodes[0],nodes[1],nodes[2],nodes[3]);
    printf("map %e %e %e %e %e %e\n",m[0],m[1],m[2],m[3],m[4],m[5]);
    return -999.0;
  }
  
  p1 = gridNodeXYZPointer(grid,nodes[0]);
  p2 = gridNodeXYZPointer(grid,nodes[1]);
  p3 = gridNodeXYZPointer(grid,nodes[2]);
  p4 = gridNodeXYZPointer(grid,nodes[3]);

  if ( gridCOST_FCN_CONFORMITY == gridCostFunction(grid) )
    return gridCellMetricConformity( p1, p2, p3, p4, m );

  xyz1[0] = j[0] * p1[0] + j[1] * p1[1] + j[2] * p1[2]; 
  xyz1[1] = j[3] * p1[0] + j[4] * p1[1] + j[5] * p1[2]; 
  xyz1[2] = j[6] * p1[0] + j[7] * p1[1] + j[8] * p1[2]; 

  xyz2[0] = j[0] * p2[0] + j[1] * p2[1] + j[2] * p2[2]; 
  xyz2[1] = j[3] * p2[0] + j[4] * p2[1] + j[5] * p2[2]; 
  xyz2[2] = j[6] * p2[0] + j[7] * p2[1] + j[8] * p2[2]; 

  xyz3[0] = j[0] * p3[0] + j[1] * p3[1] + j[2] * p3[2]; 
  xyz3[1] = j[3] * p3[0] + j[4] * p3[1] + j[5] * p3[2]; 
  xyz3[2] = j[6] * p3[0] + j[7] * p3[1] + j[8] * p3[2]; 

  xyz4[0] = j[0] * p4[0] + j[1] * p4[1] + j[2] * p4[2]; 
  xyz4[1] = j[3] * p4[0] + j[4] * p4[1] + j[5] * p4[2]; 
  xyz4[2] = j[6] * p4[0] + j[7] * p4[1] + j[8] * p4[2]; 

  switch ( gridCostFunction(grid) ) {
  case gridCOST_FCN_MEAN_RATIO:
    aspect = gridCellMeanRatio( xyz1, xyz2, xyz3, xyz4 ); break;
  case gridCOST_FCN_ASPECT_RATIO:
    aspect = gridCellAspectRatio( xyz1, xyz2, xyz3, xyz4 ); break;
  default:
    printf("%s: %d: error Cost Function %d not supported.\n",__FILE__,__LINE__,
	   gridCostFunction(grid));
    return -1.0;
  }

  return aspect;

}

double gridEdgeRatioCost(Grid *grid, int *nodes )
{
  double err[6], worstErr;
  int edge;

  err[0] = gridEdgeRatioError(grid, nodes[0], nodes[1] );
  err[1] = gridEdgeRatioError(grid, nodes[0], nodes[2] );
  err[2] = gridEdgeRatioError(grid, nodes[0], nodes[3] );
  err[3] = gridEdgeRatioError(grid, nodes[1], nodes[2] );
  err[4] = gridEdgeRatioError(grid, nodes[1], nodes[3] );
  err[5] = gridEdgeRatioError(grid, nodes[2], nodes[3] );
  
  worstErr = 0.0;
  for (edge=0;edge<6;edge++) {
    if (err[edge] < -0.5) return err[edge];
    worstErr=MAX(worstErr,err[edge]);
  }
  return 1.0/(1.0+worstErr);
}

Grid *gridCellMetricConformityFD( Grid *grid,
				  double *xyz0, double *xyz1, 
				  double *xyz2, double *xyz3,
				  double *requested_metric,
				  double *cost, double *dCostdx )
{
  double xyz[3];
  double delta;
  
  delta = 1.0e-7;

  xyz[0] = xyz0[0];
  xyz[1] = xyz0[1];
  xyz[2] = xyz0[2];

  (*cost) = gridCellMetricConformity( xyz, xyz1, xyz2, xyz3,
				      requested_metric );
  if ( *cost < -0.5 ) return NULL;
  
  xyz[0] = xyz0[0] + delta;
  xyz[1] = xyz0[1];
  xyz[2] = xyz0[2];

  dCostdx[0] = gridCellMetricConformity( xyz, xyz1, xyz2, xyz3,
					 requested_metric );
  if ( dCostdx[0] < -0.5 ) return NULL;
  dCostdx[0] = (dCostdx[0] - (*cost)) / delta;
  
  xyz[0] = xyz0[0];
  xyz[1] = xyz0[1] + delta;
  xyz[2] = xyz0[2];

  dCostdx[1] = gridCellMetricConformity( xyz, xyz1, xyz2, xyz3,
					 requested_metric );
  if ( dCostdx[1] < -0.5 ) return NULL;
  dCostdx[1] = (dCostdx[1] - (*cost)) / delta;
  
  xyz[0] = xyz0[0];
  xyz[1] = xyz0[1];
  xyz[2] = xyz0[2] + delta;

  dCostdx[2] = gridCellMetricConformity( xyz, xyz1, xyz2, xyz3,
					 requested_metric );
  if ( dCostdx[2] < -0.5 ) return NULL;
  dCostdx[2] = (dCostdx[2] - (*cost)) / delta;
  
  return grid;
}

double gridCellMetricConformity( double *xyz0, double *xyz1, 
				 double *xyz2, double *xyz3,
				 double *requested_metric )
{
  double implied_metric[6];
  double minv[6];
  double mm1[9], mm2[9];
  double rt[9];
  double norm;
  int i;

  if ( !gridImpliedMetric( xyz0, xyz1, xyz2, xyz3, implied_metric ) )
    return -1.0;

  if ( !gridInverseM( requested_metric, minv ) )
    return -1.0;
  
  gridMatrixMultiplyM( minv, implied_metric, mm1 );

  if ( !gridInverseM( implied_metric, minv ) )
    return -1.0;

  gridMatrixMultiplyM( minv, requested_metric, mm2 );

  for ( i = 0 ; i < 9 ; i++ ) rt[i] = mm1[i] + mm2[i];
  rt[0] -= 2.0;
  rt[4] -= 2.0;
  rt[8] -= 2.0;

  norm = 0.0;
  for ( i = 0 ; i < 9 ; i++ ) norm += rt[i]*rt[i];
  norm = sqrt(norm);

  norm = 1.0/(1.0+norm);

  return norm;
}

double gridCellAspectRatio( double *xyz1, double *xyz2, 
			    double *xyz3, double *xyz4 )
{
  double x1, x2, x3, x4; 
  double y1, y2, y3, y4; 
  double z1, z2, z3, z4; 
  double s1, s2, s3, s4, det;
  double xr, yr, zr;
  double circ;
  double nx1, ny1, nz1, rmag1;
  double nx2, ny2, nz2, rmag2;
  double nx3, ny3, nz3, rmag3;
  double nx4, ny4, nz4, rmag4;
  double xins;

  x1 = xyz1[0];
  y1 = xyz1[1];
  z1 = xyz1[2];

  x2 = xyz2[0];
  y2 = xyz2[1];
  z2 = xyz2[2];

  x3 = xyz3[0];
  y3 = xyz3[1];
  z3 = xyz3[2];

  x4 = xyz4[0];
  y4 = xyz4[1];
  z4 = xyz4[2];

  det = (x4-x1)*((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
      + (y4-y1)*((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
      + (z4-z1)*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
  s1 = ((x2*x2 + y2*y2 + z2*z2) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  s2 = ((x3*x3 + y3*y3 + z3*z3) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  s3 = ((x4*x4 + y4*y4 + z4*z4) - (x1*x1 + y1*y1 + z1*z1))/2.0;
  xr  =(s3     *((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
       	+ (y4-y1)*(s2     *(z2-z1) - s1     *(z3-z1))
	+ (z4-z1)*(s1     *(y3-y1) - s2     *(y2-y1)))/det;
  yr  =((x4-x1)*(s1     *(z3-z1) - s2     *(z2-z1))
	+ s3     *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
	+ (z4-z1)*((x2-x1)*s2      - (x3-x1)*s1     ))/det;
  zr  =((x4-x1)*((y2-y1)*s2      - (y3-y1)*s1     )
	+ (y4-y1)*((x3-x1)*s1      - (x2-x1)*s2     )
	+ s3     *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)))/det;
  circ = sqrt((x1-xr)*(x1-xr) + (y1-yr)*(y1-yr) + (z1-zr)*(z1-zr));

  /* Get the in-circle */

  nx1 = (y3 - y1)*(z2 - z1) - (y2 - y1)*(z3 - z1);
  ny1 =-(x3 - x1)*(z2 - z1) + (x2 - x1)*(z3 - z1);
  nz1 = (x3 - x1)*(y2 - y1) - (x2 - x1)*(y3 - y1);
  rmag1 = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);

  nx2 = (y2 - y1)*(z4 - z1) - (y4 - y1)*(z2 - z1);
  ny2 =-(x2 - x1)*(z4 - z1) + (x4 - x1)*(z2 - z1);
  nz2 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);
  rmag2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);

  nx3 = (y4 - y1)*(z3 - z1) - (y3 - y1)*(z4 - z1);
  ny3 =-(x4 - x1)*(z3 - z1) + (x3 - x1)*(z4 - z1);
  nz3 = (x4 - x1)*(y3 - y1) - (x3 - x1)*(y4 - y1);
  rmag3 = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3);

  nx4 = (y3 - y2)*(z4 - z2) - (y4 - y2)*(z3 - z2);
  ny4 =-(x3 - x2)*(z4 - z2) + (x4 - x2)*(z3 - z2);
  nz4 = (x3 - x2)*(y4 - y2) - (x4 - x2)*(y3 - y2);
  rmag4 = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4);
  nx1 = nx1/rmag1;
  ny1 = ny1/rmag1;
  nz1 = nz1/rmag1;
  nx2 = nx2/rmag2;
  ny2 = ny2/rmag2;
  nz2 = nz2/rmag2;
  nx3 = nx3/rmag3;
  ny3 = ny3/rmag3;
  nz3 = nz3/rmag3;
  nx4 = nx4/rmag4;
  ny4 = ny4/rmag4;
  nz4 = nz4/rmag4;
  det= -(nx3*ny2*nz1) + nx4*ny2*nz1 + nx2*ny3*nz1 - nx4*ny3*nz1
    -nx2*ny4*nz1 + nx3*ny4*nz1 + nx3*ny1*nz2 - nx4*ny1*nz2
    -nx1*ny3*nz2 + nx4*ny3*nz2 + nx1*ny4*nz2 - nx3*ny4*nz2
    -nx2*ny1*nz3 + nx4*ny1*nz3 + nx1*ny2*nz3 - nx4*ny2*nz3
    -nx1*ny4*nz3 + nx2*ny4*nz3 + nx2*ny1*nz4 - nx3*ny1*nz4
    -nx1*ny2*nz4 + nx3*ny2*nz4 + nx1*ny3*nz4 - nx2*ny3*nz4;
  s1 = nx1*x1 + ny1*y1 + nz1*z1;
  s2 = nx2*x1 + ny2*y1 + nz2*z1;
  s3 = nx3*x1 + ny3*y1 + nz3*z1;
  s4 = nx4*x4 + ny4*y4 + nz4*z4;
  xins = (nx4*ny3*nz2*s1 - nx3*ny4*nz2*s1 - nx4*ny2*nz3*s1 +
	  nx2*ny4*nz3*s1 +
	  nx3*ny2*nz4*s1 - nx2*ny3*nz4*s1 - nx4*ny3*nz1*s2 +
	  nx3*ny4*nz1*s2 +
	  nx4*ny1*nz3*s2 - nx1*ny4*nz3*s2 - nx3*ny1*nz4*s2 +
	  nx1*ny3*nz4*s2 +
	  nx4*ny2*nz1*s3 - nx2*ny4*nz1*s3 - nx4*ny1*nz2*s3 +
	  nx1*ny4*nz2*s3 +
	  nx2*ny1*nz4*s3 - nx1*ny2*nz4*s3 - nx3*ny2*nz1*s4 +
	  nx2*ny3*nz1*s4 +
	  nx3*ny1*nz2*s4 - nx1*ny3*nz2*s4 - nx2*ny1*nz3*s4 +
	  nx1*ny2*nz3*s4)/det;

  return xins/circ*3.0;
}

Grid *gridNodeARDerivative (Grid *grid, int node, double *ar, double *dARdx )
{
  AdjIterator it;
  int nodes[4], orientedNodes[4];
  double local_ar, local_dARdx[3];

  *ar = 1.1;
  dARdx[0] = DBL_MAX;
  dARdx[1] = DBL_MAX;
  dARdx[2] = DBL_MAX;

  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    gridCell(grid,adjItem(it),nodes);
    orientedNodes[0] = node;
    if (node == nodes[0]){
      orientedNodes[1] = nodes[1];
    }else{
      orientedNodes[1] = nodes[0];
    }
    gridOrient( grid, nodes, orientedNodes);
    if ( grid != gridCellARDerivative( grid, orientedNodes, 
				      &local_ar, local_dARdx ) ) {
      *ar = 0.0;
      dARdx[0] = DBL_MAX;
      dARdx[1] = DBL_MAX;
      dARdx[2] = DBL_MAX;
      return NULL;
    }
    if ( local_ar < *ar ) {
      *ar = local_ar;
      dARdx[0] = local_dARdx[0];
      dARdx[1] = local_dARdx[1];
      dARdx[2] = local_dARdx[2];
    }
  }

  return grid;
}
Grid *gridCellARDerivative(Grid *grid, int *nodes, double *ar, double *dARdx )
{
  double xyz1[3], xyz2[3], xyz3[3], xyz4[3]; 
  int i;
  double *m0, *m1, *m2, *m3; 
  double map0[6], map1[6], map2[6], map3[6];
  double m[6], j[9];

  if ( gridCOST_FCN_EDGE_LENGTH == gridCostFunction(grid) )
    return gridCellRatioErrorDerivative(grid, nodes, ar, dARdx );

  if (grid != gridNodeXYZ(grid,nodes[0],xyz1) ) return NULL;
  if (grid != gridNodeXYZ(grid,nodes[1],xyz2) ) return NULL;
  if (grid != gridNodeXYZ(grid,nodes[2],xyz3) ) return NULL;
  if (grid != gridNodeXYZ(grid,nodes[3],xyz4) ) return NULL;

  if ( gridMapPointerAllocated(grid) ) {
    m0 = gridMapPointer(grid,nodes[0]);
    m1 = gridMapPointer(grid,nodes[1]);
    m2 = gridMapPointer(grid,nodes[2]);
    m3 = gridMapPointer(grid,nodes[3]);

    for (i=0;i<6;i++) m[i]=0.25*(m0[i]+m1[i]+m2[i]+m3[i]);
  } else {
    gridMap(grid, nodes[0], map0);
    gridMap(grid, nodes[1], map1);
    gridMap(grid, nodes[2], map2);
    gridMap(grid, nodes[3], map3);
    for (i=0;i<6;i++) m[i]=0.25*(map0[i]+map1[i]+map2[i]+map3[i]);
  }

  if ( gridCOST_FCN_CONFORMITY == gridCostFunction(grid) ) {
    if ( grid != gridCellMetricConformityFD( grid, xyz1, xyz2, xyz3, xyz4, 
					     m, ar, dARdx ) ) return NULL;
    return grid;
  }

  if (grid != gridConvertMetricToJacobian(grid, m, j) ) {
    printf("%s: %d: gridCellARDerivative: gridConvertMetricToJacobian NULL\n",
	   __FILE__,__LINE__);
  }

  gridMapXYZWithJ(j, &xyz1[0], &xyz1[1], &xyz1[2]);
  gridMapXYZWithJ(j, &xyz2[0], &xyz2[1], &xyz2[2]);
  gridMapXYZWithJ(j, &xyz3[0], &xyz3[1], &xyz3[2]);
  gridMapXYZWithJ(j, &xyz4[0], &xyz4[1], &xyz4[2]);

  switch ( gridCostFunction(grid) ) {
  case gridCOST_FCN_MEAN_RATIO:
    gridCellMeanRatioDerivative( xyz1, xyz2, xyz3, xyz4, ar, dARdx); break;
  case gridCOST_FCN_ASPECT_RATIO:
    gridCellAspectRatioDerivative( xyz1, xyz2, xyz3, xyz4, ar, dARdx); break;
  default:
    printf("%s: %d: error Cost Function %d not supported.\n",__FILE__,__LINE__,
	   gridCostFunction(grid));
    return NULL;
  }

  return grid;
}

Grid *gridCellRatioErrorDerivative(Grid *grid, int *nodes, 
				   double *cost, double *dCostdx )
{
  int edge, worstEdge;
  double err, worstErr;
  double xyz0[3], xyz1[3], direction[3];
  double ratio;

  *cost = gridEdgeRatioCost(grid,nodes);
  dCostdx[0] = 0.0;
  dCostdx[1] = 0.0;
  dCostdx[2] = 0.0;

  worstErr = -1.0;
  worstEdge = EMPTY;
  for (edge=0;edge<3;edge++){
    err = gridEdgeRatioError(grid,nodes[0], nodes[edge+1]);
    if (err >= worstErr) {
      worstEdge=edge;
      worstErr =err;
    }
  }

  if (worstEdge==EMPTY) return NULL;
  
  if (grid!=gridNodeXYZ(grid,nodes[0],xyz0))return NULL;
  if (grid!=gridNodeXYZ(grid,nodes[worstEdge+1],xyz1))return NULL;
  gridSubtractVector(xyz1,xyz0,direction);
  ratio = gridEdgeRatio(grid,nodes[0], nodes[worstEdge+1]);
  if (ratio>=1.0) {
    gridVectorScale(direction,ratio-1.0);
  }else{
    gridVectorScale(direction,-(1.0-ratio));
  }    
  
  gridVectorCopy(dCostdx,direction);

  return grid;
}

void gridCellAspectRatioDerivative( double *xyz1, double *xyz2, 
				    double *xyz3, double *xyz4,
				    double *ar, double *dARdx)
{
  double x1, x2, x3, x4; 
  double y1, y2, y3, y4; 
  double z1, z2, z3, z4; 
  double s1, s2, s3, s4, det;
  double xr, yr, zr;
  double circ;
  double nx1, ny1, nz1, rmag1;
  double nx2, ny2, nz2, rmag2;
  double nx3, ny3, nz3, rmag3;
  double nx4, ny4, nz4, rmag4;
  double xins;

  double s1_dx, s2_dx, s3_dx, s4_dx, det_dx;
  double s1_dy, s2_dy, s3_dy, s4_dy, det_dy;
  double s1_dz, s2_dz, s3_dz, s4_dz, det_dz;
  double xr_dx, yr_dx, zr_dx;
  double xr_dy, yr_dy, zr_dy;
  double xr_dz, yr_dz, zr_dz;
  double circ_dx, circ_dy, circ_dz; 

  double nx1_dx, ny1_dx, nz1_dx, rmag1_dx;
  double nx1_dy, ny1_dy, nz1_dy, rmag1_dy;
  double nx1_dz, ny1_dz, nz1_dz, rmag1_dz;

  double nx2_dx, ny2_dx, nz2_dx, rmag2_dx;
  double nx2_dy, ny2_dy, nz2_dy, rmag2_dy;
  double nx2_dz, ny2_dz, nz2_dz, rmag2_dz;

  double nx3_dx, ny3_dx, nz3_dx, rmag3_dx;
  double nx3_dy, ny3_dy, nz3_dy, rmag3_dy;
  double nx3_dz, ny3_dz, nz3_dz, rmag3_dz;

  double nx4_dx, ny4_dx, nz4_dx;
  double nx4_dy, ny4_dy, nz4_dy;
  double nx4_dz, ny4_dz, nz4_dz;

  double xins_dx, xins_dy, xins_dz;

  x1 = xyz1[0];
  y1 = xyz1[1];
  z1 = xyz1[2];

  x2 = xyz2[0];
  y2 = xyz2[1];
  z2 = xyz2[2];

  x3 = xyz3[0];
  y3 = xyz3[1];
  z3 = xyz3[2];

  x4 = xyz4[0];
  y4 = xyz4[1];
  z4 = xyz4[2];

        det = (x4-x1)*((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
            + (y4-y1)*((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
            + (z4-z1)*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

	det_dx = -((y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
   	       + (y4-y1)*(-(z2-z1) + (z3-z1))
	       + (z4-z1)*(-(y3-y1) + (y2-y1));

        det_dy = (x4-x1)*(-(z3-z1) +(z2-z1))
            - ((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
            + (z4-z1)*(-(x2-x1) + (x3-x1));

        det_dz = (x4-x1)*(-(y2-y1) + (y3-y1))
            + (y4-y1)*(-(x3-x1) + (x2-x1))
            -((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

         s1 = ((x2*x2 + y2*y2 + z2*z2) - (x1*x1 + y1*y1 + z1*z1))/2.0;
         s2 = ((x3*x3 + y3*y3 + z3*z3) - (x1*x1 + y1*y1 + z1*z1))/2.0;
         s3 = ((x4*x4 + y4*y4 + z4*z4) - (x1*x1 + y1*y1 + z1*z1))/2.0;

	 s1_dx = -x1;
	 s1_dy = -y1;
	 s1_dz = -z1;
	 s2_dx = -x1;
	 s2_dy = -y1;
	 s2_dz = -z1;
	 s3_dx = -x1;
	 s3_dy = -y1;
	 s3_dz = -z1;


         xr  = (   s3    * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
               + (y4-y1) * (    s2  *(z2-z1) -   s1   *(z3-z1))
               + (z4-z1) * (    s1  *(y3-y1) -   s2   *(y2-y1)) ) / det;


	 xr_dx = (  s3_dx  * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
                 + (y4-y1) * (  s2_dx *(z2-z1) -  s1_dx *(z3-z1))
                 + (z4-z1) * (  s1_dx *(y3-y1) -  s2_dx *(y2-y1)) ) ;

	 xr_dx = det * xr_dx - (xr*det) * det_dx ;
	 xr_dx = xr_dx / det / det;

	 xr_dy = (   s3    * (        -(z3-z1) +         (z2-z1))
                 +  s3_dy  * ( (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1))
                 + (y4-y1) * (  s2_dy *(z2-z1) -  s1_dy *(z3-z1))
                 -           (    s2  *(z2-z1) -   s1   *(z3-z1))
                 + (z4-z1) * ( s1_dy  *(y3-y1) -  s2_dy *(y2-y1)) 
                 + (z4-z1) * (        -s1      +        s2     ) ) ;

	 xr_dy = det * xr_dy - (xr*det) * det_dy ;
	 xr_dy = xr_dy / det / det;


         xr_dz = (   s3    * (           -(y2-y1) +            (y3-y1))
	         +  s3_dz  * (    (y2-y1)*(z3-z1) -    (y3-y1)*(z2-z1))
                 + (y4-y1) * (    -s2             +    s1             )
                 + (y4-y1) * (    s2_dz  *(z2-z1) -   s1_dz   *(z3-z1))
                 + (z4-z1) * (    s1_dz  *(y3-y1) -   s2_dz   *(y2-y1)) 
                 -           (    s1     *(y3-y1) -   s2      *(y2-y1)) );
	 xr_dz = det * xr_dz - (xr*det) * det_dz ;
	 xr_dz = xr_dz / det / det;

         yr  =((x4-x1)*(s1     *(z3-z1) - s2     *(z2-z1))
             + s3     *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
             + (z4-z1)*((x2-x1)*s2      - (x3-x1)*s1     ))/det;

         yr_dx = ( (x4-x1)*(s1_dx  *(z3-z1) - s2_dx  *(z2-z1))
		 -         (s1     *(z3-z1) - s2     *(z2-z1))
                 + s3     *(       -(z2-z1)          +(z3-z1))
                 + s3_dx  *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                 + (z4-z1)*(       -s2               +s1     )
                 + (z4-z1)*((x2-x1)*s2_dx   - (x3-x1)*s1_dx  ));

	 yr_dx = det * yr_dx - (yr*det) * det_dx ;
	 yr_dx = yr_dx / det / det;

         yr_dy  = ( (x4-x1)*(s1_dy  *(z3-z1) - s2_dy  *(z2-z1))
                  +  s3_dy *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                  + (z4-z1)*((x2-x1)*s2_dy   - (x3-x1)*s1_dy  ) );

	 yr_dy = det * yr_dy - (yr*det) * det_dy ;
	 yr_dy = yr_dy / det / det;


         yr_dz = ( (x4-x1)*(  -s1           + s2             )
		 + (x4-x1)*( s1_dz *(z3-z1) - s2_dz  *(z2-z1))
                 + s3     *(-(x3-x1)        + (x2-x1)        )
                 + s3_dz  *((x3-x1)*(z2-z1) - (x2-x1)*(z3-z1))
                 + (z4-z1)*((x2-x1)*s2_dz   - (x3-x1)*s1_dz  )
                 -         ((x2-x1)*s2      - (x3-x1)*s1     ) );

	 yr_dz = det * yr_dz - (yr*det) * det_dz ;
	 yr_dz = yr_dz / det / det;

         zr  =((x4-x1)*((y2-y1)*s2      - (y3-y1)*s1     )
             + (y4-y1)*((x3-x1)*s1      - (x2-x1)*s2     )
             + s3     *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)))/det;

         zr_dx = ( (x4-x1)*((y2-y1)*s2_dx   - (y3-y1)*s1_dx  )
		 -         ((y2-y1)*s2      - (y3-y1)*s1     )
                 + (y4-y1)*(       -s1      +         s2     )
                 + (y4-y1)*((x3-x1)*s1_dx   - (x2-x1)*s2_dx  )
                 + s3     *(       -(y3-y1) +         (y2-y1))
                 + s3_dx  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dx = det * zr_dx - (zr*det) * det_dx ;
	 zr_dx = zr_dx / det / det;

         zr_dy = ( (x4-x1)*((y2-y1)*s2_dy   - (y3-y1)*s1_dy  )
		 + (x4-x1)*(       -s2      +         s1     )
                 + (y4-y1)*((x3-x1)*s1_dy   - (x2-x1)*s2_dy  )
                 -         ((x3-x1)*s1      - (x2-x1)*s2     )
                 + s3     *(-(x2-x1)        + (x3-x1)        )
                 + s3_dy  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dy = det * zr_dy - (zr*det) * det_dy ;
	 zr_dy = zr_dy / det / det;

         zr_dz = ( (x4-x1)*((y2-y1)*s2_dz   - (y3-y1)*s1_dz  )
                 + (y4-y1)*((x3-x1)*s1_dz   - (x2-x1)*s2_dz  )
                 + s3_dz  *((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );

	 zr_dz = det * zr_dz - (zr*det) * det_dz ;
	 zr_dz = zr_dz / det / det;

         circ = sqrt((x1-xr)*(x1-xr) + (y1-yr)*(y1-yr) + (z1-zr)*(z1-zr));

	 circ_dx = 2.0*(x1-xr)*(1.0-xr_dx) 
	         + 2.0*(y1-yr)*(   -yr_dx) 
	         + 2.0*(z1-zr)*(   -zr_dx);

	 circ_dx = 0.5 / circ * circ_dx;

	 circ_dy = 2.0*(x1-xr)*(   -xr_dy) 
	         + 2.0*(y1-yr)*(1.0-yr_dy) 
	         + 2.0*(z1-zr)*(   -zr_dy);

	 circ_dy = 0.5 / circ * circ_dy;

	 circ_dz = 2.0*(x1-xr)*(   -xr_dz) 
	         + 2.0*(y1-yr)*(   -yr_dz) 
	         + 2.0*(z1-zr)*(1.0-zr_dz);

	 circ_dz = 0.5 / circ * circ_dz;
 
  /* Get the in-circle */

         nx1 = (y3 - y1)*(z2 - z1) - (y2 - y1)*(z3 - z1);
         ny1 =-(x3 - x1)*(z2 - z1) + (x2 - x1)*(z3 - z1);
         nz1 = (x3 - x1)*(y2 - y1) - (x2 - x1)*(y3 - y1);

	 nx1_dx = 0.0;
	 ny1_dx = (z2 - z1) - (z3 - z1);
         nz1_dx =-(y2 - y1) + (y3 - y1);

	 nx1_dy = -(z2 - z1) + (z3 - z1);
	 ny1_dy = 0.0;
         nz1_dy = -(x3 - x1) + (x2 - x1);
 
	 nx1_dz = -(y3 - y1) + (y2 - y1);
	 ny1_dz =  (x3 - x1) - (x2 - x1);
         nz1_dz = 0.0;

         nx2 = (y2 - y1)*(z4 - z1) - (y4 - y1)*(z2 - z1);
         ny2 =-(x2 - x1)*(z4 - z1) + (x4 - x1)*(z2 - z1);
         nz2 = (x2 - x1)*(y4 - y1) - (x4 - x1)*(y2 - y1);

         nx2_dx =  0.0;
         ny2_dx =  (z4 - z1) - (z2 - z1);
         nz2_dx = -(y4 - y1) + (y2 - y1);

         nx2_dy = -(z4 - z1) + (z2 - z1);
         ny2_dy =  0.0;
         nz2_dy = -(x2 - x1) + (x4 - x1);

         nx2_dz = -(y2 - y1) + (y4 - y1);
         ny2_dz =  (x2 - x1) - (x4 - x1);
         nz2_dz =  0.0;

         nx3 = (y4 - y1)*(z3 - z1) - (y3 - y1)*(z4 - z1);
         ny3 =-(x4 - x1)*(z3 - z1) + (x3 - x1)*(z4 - z1);
         nz3 = (x4 - x1)*(y3 - y1) - (x3 - x1)*(y4 - y1);

         nx3_dx =  0.0;
         ny3_dx =  (z3 - z1) - (z4 - z1);
         nz3_dx = -(y3 - y1) + (y4 - y1);

         nx3_dy = -(z3 - z1) + (z4 - z1);
         ny3_dy = 0.0;
         nz3_dy = -(x4 - x1) + (x3 - x1);

         nx3_dz =-(y4 - y1) + (y3 - y1);
         ny3_dz = (x4 - x1) - (x3 - x1);
         nz3_dz = 0.0;

         nx4 = (y3 - y2)*(z4 - z2) - (y4 - y2)*(z3 - z2);
         ny4 =-(x3 - x2)*(z4 - z2) + (x4 - x2)*(z3 - z2);
         nz4 = (x3 - x2)*(y4 - y2) - (x4 - x2)*(y3 - y2);

         rmag1 = sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);

	 rmag1_dx = (nx1*nx1_dx + ny1*ny1_dx + nz1*nz1_dx) / rmag1;
	 rmag1_dy = (nx1*nx1_dy + ny1*ny1_dy + nz1*nz1_dy) / rmag1;
	 rmag1_dz = (nx1*nx1_dz + ny1*ny1_dz + nz1*nz1_dz) / rmag1;

         rmag2 = sqrt(nx2*nx2 + ny2*ny2 + nz2*nz2);

	 rmag2_dx = (nx2*nx2_dx + ny2*ny2_dx + nz2*nz2_dx) / rmag2;
	 rmag2_dy = (nx2*nx2_dy + ny2*ny2_dy + nz2*nz2_dy) / rmag2;
	 rmag2_dz = (nx2*nx2_dz + ny2*ny2_dz + nz2*nz2_dz) / rmag2;

         rmag3 = sqrt(nx3*nx3 + ny3*ny3 + nz3*nz3);

	 rmag3_dx = (nx3*nx3_dx + ny3*ny3_dx + nz3*nz3_dx) / rmag3;
	 rmag3_dy = (nx3*nx3_dy + ny3*ny3_dy + nz3*nz3_dy) / rmag3;
	 rmag3_dz = (nx3*nx3_dz + ny3*ny3_dz + nz3*nz3_dz) / rmag3;

         rmag4 = sqrt(nx4*nx4 + ny4*ny4 + nz4*nz4);
	
	 /* note that the derivatives procede the normalization */
	 /* because I need the un-nornmalized value for deriavtives */
 
         nx1_dx = ( rmag1 * nx1_dx - nx1 * rmag1_dx ) / rmag1 / rmag1;
         ny1_dx = ( rmag1 * ny1_dx - ny1 * rmag1_dx ) / rmag1 / rmag1;
         nz1_dx = ( rmag1 * nz1_dx - nz1 * rmag1_dx ) / rmag1 / rmag1;

         nx1_dy = ( rmag1 * nx1_dy - nx1 * rmag1_dy ) / rmag1 / rmag1;
         ny1_dy = ( rmag1 * ny1_dy - ny1 * rmag1_dy ) / rmag1 / rmag1;
         nz1_dy = ( rmag1 * nz1_dy - nz1 * rmag1_dy ) / rmag1 / rmag1;

         nx1_dz = ( rmag1 * nx1_dz - nx1 * rmag1_dz ) / rmag1 / rmag1;
         ny1_dz = ( rmag1 * ny1_dz - ny1 * rmag1_dz ) / rmag1 / rmag1;
         nz1_dz = ( rmag1 * nz1_dz - nz1 * rmag1_dz ) / rmag1 / rmag1;

         nx1 = nx1/rmag1;
         ny1 = ny1/rmag1;
         nz1 = nz1/rmag1;

         nx2_dx = ( rmag2 * nx2_dx - nx2 * rmag2_dx ) / rmag2 / rmag2;
         ny2_dx = ( rmag2 * ny2_dx - ny2 * rmag2_dx ) / rmag2 / rmag2;
         nz2_dx = ( rmag2 * nz2_dx - nz2 * rmag2_dx ) / rmag2 / rmag2;

         nx2_dy = ( rmag2 * nx2_dy - nx2 * rmag2_dy ) / rmag2 / rmag2;
         ny2_dy = ( rmag2 * ny2_dy - ny2 * rmag2_dy ) / rmag2 / rmag2;
         nz2_dy = ( rmag2 * nz2_dy - nz2 * rmag2_dy ) / rmag2 / rmag2;

         nx2_dz = ( rmag2 * nx2_dz - nx2 * rmag2_dz ) / rmag2 / rmag2;
         ny2_dz = ( rmag2 * ny2_dz - ny2 * rmag2_dz ) / rmag2 / rmag2;
         nz2_dz = ( rmag2 * nz2_dz - nz2 * rmag2_dz ) / rmag2 / rmag2;

         nx2 = nx2/rmag2;
         ny2 = ny2/rmag2;
         nz2 = nz2/rmag2;

         nx3_dx = ( rmag3 * nx3_dx - nx3 * rmag3_dx ) / rmag3 / rmag3;
         ny3_dx = ( rmag3 * ny3_dx - ny3 * rmag3_dx ) / rmag3 / rmag3;
         nz3_dx = ( rmag3 * nz3_dx - nz3 * rmag3_dx ) / rmag3 / rmag3;

         nx3_dy = ( rmag3 * nx3_dy - nx3 * rmag3_dy ) / rmag3 / rmag3;
         ny3_dy = ( rmag3 * ny3_dy - ny3 * rmag3_dy ) / rmag3 / rmag3;
         nz3_dy = ( rmag3 * nz3_dy - nz3 * rmag3_dy ) / rmag3 / rmag3;

         nx3_dz = ( rmag3 * nx3_dz - nx3 * rmag3_dz ) / rmag3 / rmag3;
         ny3_dz = ( rmag3 * ny3_dz - ny3 * rmag3_dz ) / rmag3 / rmag3;
         nz3_dz = ( rmag3 * nz3_dz - nz3 * rmag3_dz ) / rmag3 / rmag3;

         nx3 = nx3/rmag3;
         ny3 = ny3/rmag3;
         nz3 = nz3/rmag3;

         nx4_dx = 0.0;
         ny4_dx = 0.0;
         nz4_dx = 0.0;

         nx4_dy = 0.0;
         ny4_dy = 0.0;
         nz4_dy = 0.0;

         nx4_dz = 0.0;
         ny4_dz = 0.0;
         nz4_dz = 0.0;

         nx4 = nx4/rmag4;
         ny4 = ny4/rmag4;
         nz4 = nz4/rmag4;

         det= -nx3*ny2*nz1 + nx4*ny2*nz1 + nx2*ny3*nz1 - nx4*ny3*nz1
              -nx2*ny4*nz1 + nx3*ny4*nz1 + nx3*ny1*nz2 - nx4*ny1*nz2
              -nx1*ny3*nz2 + nx4*ny3*nz2 + nx1*ny4*nz2 - nx3*ny4*nz2
              -nx2*ny1*nz3 + nx4*ny1*nz3 + nx1*ny2*nz3 - nx4*ny2*nz3
              -nx1*ny4*nz3 + nx2*ny4*nz3 + nx2*ny1*nz4 - nx3*ny1*nz4
              -nx1*ny2*nz4 + nx3*ny2*nz4 + nx1*ny3*nz4 - nx2*ny3*nz4;

         det_dx = 
	   - nx3_dx*ny2*nz1 - nx3*ny2_dx*nz1 - nx3*ny2*nz1_dx 
	   + nx4_dx*ny2*nz1 + nx4*ny2_dx*nz1 + nx4*ny2*nz1_dx 
	   + nx2_dx*ny3*nz1 + nx2*ny3_dx*nz1 + nx2*ny3*nz1_dx 
	   - nx4_dx*ny3*nz1 - nx4*ny3_dx*nz1 - nx4*ny3*nz1_dx
	   - nx2_dx*ny4*nz1 - nx2*ny4_dx*nz1 - nx2*ny4*nz1_dx 
	   + nx3_dx*ny4*nz1 + nx3*ny4_dx*nz1 + nx3*ny4*nz1_dx
	   + nx3_dx*ny1*nz2 + nx3*ny1_dx*nz2 + nx3*ny1*nz2_dx 
	   - nx4_dx*ny1*nz2 - nx4*ny1_dx*nz2 - nx4*ny1*nz2_dx
	   - nx1_dx*ny3*nz2 - nx1*ny3_dx*nz2 - nx1*ny3*nz2_dx 
	   + nx4_dx*ny3*nz2 + nx4*ny3_dx*nz2 + nx4*ny3*nz2_dx 
	   + nx1_dx*ny4*nz2 + nx1*ny4_dx*nz2 + nx1*ny4*nz2_dx 
	   - nx3_dx*ny4*nz2 - nx3*ny4_dx*nz2 - nx3*ny4*nz2_dx
	   - nx2_dx*ny1*nz3 - nx2*ny1_dx*nz3 - nx2*ny1*nz3_dx 
	   + nx4_dx*ny1*nz3 + nx4*ny1_dx*nz3 + nx4*ny1*nz3_dx 
	   + nx1_dx*ny2*nz3 + nx1*ny2_dx*nz3 + nx1*ny2*nz3_dx 
	   - nx4_dx*ny2*nz3 - nx4*ny2_dx*nz3 - nx4*ny2*nz3_dx
	   - nx1_dx*ny4*nz3 - nx1*ny4_dx*nz3 - nx1*ny4*nz3_dx 
	   + nx2_dx*ny4*nz3 + nx2*ny4_dx*nz3 + nx2*ny4*nz3_dx 
	   + nx2_dx*ny1*nz4 + nx2*ny1_dx*nz4 + nx2*ny1*nz4_dx 
	   - nx3_dx*ny1*nz4 - nx3*ny1_dx*nz4 - nx3*ny1*nz4_dx
	   - nx1_dx*ny2*nz4 - nx1*ny2_dx*nz4 - nx1*ny2*nz4_dx 
	   + nx3_dx*ny2*nz4 + nx3*ny2_dx*nz4 + nx3*ny2*nz4_dx 
	   + nx1_dx*ny3*nz4 + nx1*ny3_dx*nz4 + nx1*ny3*nz4_dx 
	   - nx2_dx*ny3*nz4 - nx2*ny3_dx*nz4 - nx2*ny3*nz4_dx;

         det_dy = 
	   - nx3_dy*ny2*nz1 - nx3*ny2_dy*nz1 - nx3*ny2*nz1_dy 
	   + nx4_dy*ny2*nz1 + nx4*ny2_dy*nz1 + nx4*ny2*nz1_dy 
	   + nx2_dy*ny3*nz1 + nx2*ny3_dy*nz1 + nx2*ny3*nz1_dy 
	   - nx4_dy*ny3*nz1 - nx4*ny3_dy*nz1 - nx4*ny3*nz1_dy
	   - nx2_dy*ny4*nz1 - nx2*ny4_dy*nz1 - nx2*ny4*nz1_dy 
	   + nx3_dy*ny4*nz1 + nx3*ny4_dy*nz1 + nx3*ny4*nz1_dy
	   + nx3_dy*ny1*nz2 + nx3*ny1_dy*nz2 + nx3*ny1*nz2_dy 
	   - nx4_dy*ny1*nz2 - nx4*ny1_dy*nz2 - nx4*ny1*nz2_dy
	   - nx1_dy*ny3*nz2 - nx1*ny3_dy*nz2 - nx1*ny3*nz2_dy 
	   + nx4_dy*ny3*nz2 + nx4*ny3_dy*nz2 + nx4*ny3*nz2_dy 
	   + nx1_dy*ny4*nz2 + nx1*ny4_dy*nz2 + nx1*ny4*nz2_dy 
	   - nx3_dy*ny4*nz2 - nx3*ny4_dy*nz2 - nx3*ny4*nz2_dy
	   - nx2_dy*ny1*nz3 - nx2*ny1_dy*nz3 - nx2*ny1*nz3_dy 
	   + nx4_dy*ny1*nz3 + nx4*ny1_dy*nz3 + nx4*ny1*nz3_dy 
	   + nx1_dy*ny2*nz3 + nx1*ny2_dy*nz3 + nx1*ny2*nz3_dy 
	   - nx4_dy*ny2*nz3 - nx4*ny2_dy*nz3 - nx4*ny2*nz3_dy
	   - nx1_dy*ny4*nz3 - nx1*ny4_dy*nz3 - nx1*ny4*nz3_dy 
	   + nx2_dy*ny4*nz3 + nx2*ny4_dy*nz3 + nx2*ny4*nz3_dy 
	   + nx2_dy*ny1*nz4 + nx2*ny1_dy*nz4 + nx2*ny1*nz4_dy 
	   - nx3_dy*ny1*nz4 - nx3*ny1_dy*nz4 - nx3*ny1*nz4_dy
	   - nx1_dy*ny2*nz4 - nx1*ny2_dy*nz4 - nx1*ny2*nz4_dy 
	   + nx3_dy*ny2*nz4 + nx3*ny2_dy*nz4 + nx3*ny2*nz4_dy 
	   + nx1_dy*ny3*nz4 + nx1*ny3_dy*nz4 + nx1*ny3*nz4_dy 
	   - nx2_dy*ny3*nz4 - nx2*ny3_dy*nz4 - nx2*ny3*nz4_dy;

         det_dz = 
	   - nx3_dz*ny2*nz1 - nx3*ny2_dz*nz1 - nx3*ny2*nz1_dz 
	   + nx4_dz*ny2*nz1 + nx4*ny2_dz*nz1 + nx4*ny2*nz1_dz 
	   + nx2_dz*ny3*nz1 + nx2*ny3_dz*nz1 + nx2*ny3*nz1_dz 
	   - nx4_dz*ny3*nz1 - nx4*ny3_dz*nz1 - nx4*ny3*nz1_dz
	   - nx2_dz*ny4*nz1 - nx2*ny4_dz*nz1 - nx2*ny4*nz1_dz 
	   + nx3_dz*ny4*nz1 + nx3*ny4_dz*nz1 + nx3*ny4*nz1_dz
	   + nx3_dz*ny1*nz2 + nx3*ny1_dz*nz2 + nx3*ny1*nz2_dz 
	   - nx4_dz*ny1*nz2 - nx4*ny1_dz*nz2 - nx4*ny1*nz2_dz
	   - nx1_dz*ny3*nz2 - nx1*ny3_dz*nz2 - nx1*ny3*nz2_dz 
	   + nx4_dz*ny3*nz2 + nx4*ny3_dz*nz2 + nx4*ny3*nz2_dz 
	   + nx1_dz*ny4*nz2 + nx1*ny4_dz*nz2 + nx1*ny4*nz2_dz 
	   - nx3_dz*ny4*nz2 - nx3*ny4_dz*nz2 - nx3*ny4*nz2_dz
	   - nx2_dz*ny1*nz3 - nx2*ny1_dz*nz3 - nx2*ny1*nz3_dz 
	   + nx4_dz*ny1*nz3 + nx4*ny1_dz*nz3 + nx4*ny1*nz3_dz 
	   + nx1_dz*ny2*nz3 + nx1*ny2_dz*nz3 + nx1*ny2*nz3_dz 
	   - nx4_dz*ny2*nz3 - nx4*ny2_dz*nz3 - nx4*ny2*nz3_dz
	   - nx1_dz*ny4*nz3 - nx1*ny4_dz*nz3 - nx1*ny4*nz3_dz 
	   + nx2_dz*ny4*nz3 + nx2*ny4_dz*nz3 + nx2*ny4*nz3_dz 
	   + nx2_dz*ny1*nz4 + nx2*ny1_dz*nz4 + nx2*ny1*nz4_dz 
	   - nx3_dz*ny1*nz4 - nx3*ny1_dz*nz4 - nx3*ny1*nz4_dz
	   - nx1_dz*ny2*nz4 - nx1*ny2_dz*nz4 - nx1*ny2*nz4_dz 
	   + nx3_dz*ny2*nz4 + nx3*ny2_dz*nz4 + nx3*ny2*nz4_dz 
	   + nx1_dz*ny3*nz4 + nx1*ny3_dz*nz4 + nx1*ny3*nz4_dz 
	   - nx2_dz*ny3*nz4 - nx2*ny3_dz*nz4 - nx2*ny3*nz4_dz;

         s1 = nx1*x1 + ny1*y1 + nz1*z1;

         s1_dx = nx1_dx*x1 + nx1 + ny1_dx*y1       + nz1_dx*z1;
         s1_dy = nx1_dy*x1       + ny1_dy*y1 + ny1 + nz1_dy*z1;
         s1_dz = nx1_dz*x1       + ny1_dz*y1       + nz1_dz*z1 + nz1;

         s2 = nx2*x1 + ny2*y1 + nz2*z1;

         s2_dx = nx2_dx*x1 + nx2 + ny2_dx*y1       + nz2_dx*z1;
         s2_dy = nx2_dy*x1       + ny2_dy*y1 + ny2 + nz2_dy*z1;
         s2_dz = nx2_dz*x1       + ny2_dz*y1       + nz2_dz*z1 + nz2;

         s3 = nx3*x1 + ny3*y1 + nz3*z1;

         s3_dx = nx3_dx*x1 + nx3 + ny3_dx*y1       + nz3_dx*z1;
         s3_dy = nx3_dy*x1       + ny3_dy*y1 + ny3 + nz3_dy*z1;
         s3_dz = nx3_dz*x1       + ny3_dz*y1       + nz3_dz*z1 + nz3;

         s4 = nx4*x4 + ny4*y4 + nz4*z4;

         s4_dx = 0.0;
         s4_dy = 0.0;
         s4_dz = 0.0;

         xins = nx4*ny3*nz2*s1 - nx3*ny4*nz2*s1 - nx4*ny2*nz3*s1 +
                nx2*ny4*nz3*s1 +
                nx3*ny2*nz4*s1 - nx2*ny3*nz4*s1 - nx4*ny3*nz1*s2 +
                nx3*ny4*nz1*s2 +
                nx4*ny1*nz3*s2 - nx1*ny4*nz3*s2 - nx3*ny1*nz4*s2 +
                nx1*ny3*nz4*s2 +
                nx4*ny2*nz1*s3 - nx2*ny4*nz1*s3 - nx4*ny1*nz2*s3 +
                nx1*ny4*nz2*s3 +
                nx2*ny1*nz4*s3 - nx1*ny2*nz4*s3 - nx3*ny2*nz1*s4 +
                nx2*ny3*nz1*s4 +
                nx3*ny1*nz2*s4 - nx1*ny3*nz2*s4 - nx2*ny1*nz3*s4 +
                nx1*ny2*nz3*s4;

         xins_dx = 
  nx4_dx*ny3*nz2*s1 + nx4*ny3_dx*nz2*s1 + nx4*ny3*nz2_dx*s1 + nx4*ny3*nz2*s1_dx 
- nx3_dx*ny4*nz2*s1 - nx3*ny4_dx*nz2*s1 - nx3*ny4*nz2_dx*s1 - nx3*ny4*nz2*s1_dx 
- nx4_dx*ny2*nz3*s1 - nx4*ny2_dx*nz3*s1 - nx4*ny2*nz3_dx*s1 - nx4*ny2*nz3*s1_dx 
+ nx2_dx*ny4*nz3*s1 + nx2*ny4_dx*nz3*s1 + nx2*ny4*nz3_dx*s1 + nx2*ny4*nz3*s1_dx 
+ nx3_dx*ny2*nz4*s1 + nx3*ny2_dx*nz4*s1 + nx3*ny2*nz4_dx*s1 + nx3*ny2*nz4*s1_dx 
- nx2_dx*ny3*nz4*s1 - nx2*ny3_dx*nz4*s1 - nx2*ny3*nz4_dx*s1 - nx2*ny3*nz4*s1_dx 
- nx4_dx*ny3*nz1*s2 - nx4*ny3_dx*nz1*s2 - nx4*ny3*nz1_dx*s2 - nx4*ny3*nz1*s2_dx 
+ nx3_dx*ny4*nz1*s2 + nx3*ny4_dx*nz1*s2 + nx3*ny4*nz1_dx*s2 + nx3*ny4*nz1*s2_dx 
+ nx4_dx*ny1*nz3*s2 + nx4*ny1_dx*nz3*s2 + nx4*ny1*nz3_dx*s2 + nx4*ny1*nz3*s2_dx 
- nx1_dx*ny4*nz3*s2 - nx1*ny4_dx*nz3*s2 - nx1*ny4*nz3_dx*s2 - nx1*ny4*nz3*s2_dx 
- nx3_dx*ny1*nz4*s2 - nx3*ny1_dx*nz4*s2 - nx3*ny1*nz4_dx*s2 - nx3*ny1*nz4*s2_dx 
+ nx1_dx*ny3*nz4*s2 + nx1*ny3_dx*nz4*s2 + nx1*ny3*nz4_dx*s2 + nx1*ny3*nz4*s2_dx 
+ nx4_dx*ny2*nz1*s3 + nx4*ny2_dx*nz1*s3 + nx4*ny2*nz1_dx*s3 + nx4*ny2*nz1*s3_dx 
- nx2_dx*ny4*nz1*s3 - nx2*ny4_dx*nz1*s3 - nx2*ny4*nz1_dx*s3 - nx2*ny4*nz1*s3_dx 
- nx4_dx*ny1*nz2*s3 - nx4*ny1_dx*nz2*s3 - nx4*ny1*nz2_dx*s3 - nx4*ny1*nz2*s3_dx 
+ nx1_dx*ny4*nz2*s3 + nx1*ny4_dx*nz2*s3 + nx1*ny4*nz2_dx*s3 + nx1*ny4*nz2*s3_dx 
+ nx2_dx*ny1*nz4*s3 + nx2*ny1_dx*nz4*s3 + nx2*ny1*nz4_dx*s3 + nx2*ny1*nz4*s3_dx 
- nx1_dx*ny2*nz4*s3 - nx1*ny2_dx*nz4*s3 - nx1*ny2*nz4_dx*s3 - nx1*ny2*nz4*s3_dx 
- nx3_dx*ny2*nz1*s4 - nx3*ny2_dx*nz1*s4 - nx3*ny2*nz1_dx*s4 - nx3*ny2*nz1*s4_dx 
+ nx2_dx*ny3*nz1*s4 + nx2*ny3_dx*nz1*s4 + nx2*ny3*nz1_dx*s4 + nx2*ny3*nz1*s4_dx 
+ nx3_dx*ny1*nz2*s4 + nx3*ny1_dx*nz2*s4 + nx3*ny1*nz2_dx*s4 + nx3*ny1*nz2*s4_dx 
- nx1_dx*ny3*nz2*s4 - nx1*ny3_dx*nz2*s4 - nx1*ny3*nz2_dx*s4 - nx1*ny3*nz2*s4_dx 
- nx2_dx*ny1*nz3*s4 - nx2*ny1_dx*nz3*s4 - nx2*ny1*nz3_dx*s4 - nx2*ny1*nz3*s4_dx 
+ nx1_dx*ny2*nz3*s4 + nx1*ny2_dx*nz3*s4 + nx1*ny2*nz3_dx*s4 + nx1*ny2*nz3*s4_dx;

         xins_dy = 
  nx4_dy*ny3*nz2*s1 + nx4*ny3_dy*nz2*s1 + nx4*ny3*nz2_dy*s1 + nx4*ny3*nz2*s1_dy 
- nx3_dy*ny4*nz2*s1 - nx3*ny4_dy*nz2*s1 - nx3*ny4*nz2_dy*s1 - nx3*ny4*nz2*s1_dy 
- nx4_dy*ny2*nz3*s1 - nx4*ny2_dy*nz3*s1 - nx4*ny2*nz3_dy*s1 - nx4*ny2*nz3*s1_dy 
+ nx2_dy*ny4*nz3*s1 + nx2*ny4_dy*nz3*s1 + nx2*ny4*nz3_dy*s1 + nx2*ny4*nz3*s1_dy 
+ nx3_dy*ny2*nz4*s1 + nx3*ny2_dy*nz4*s1 + nx3*ny2*nz4_dy*s1 + nx3*ny2*nz4*s1_dy 
- nx2_dy*ny3*nz4*s1 - nx2*ny3_dy*nz4*s1 - nx2*ny3*nz4_dy*s1 - nx2*ny3*nz4*s1_dy 
- nx4_dy*ny3*nz1*s2 - nx4*ny3_dy*nz1*s2 - nx4*ny3*nz1_dy*s2 - nx4*ny3*nz1*s2_dy 
+ nx3_dy*ny4*nz1*s2 + nx3*ny4_dy*nz1*s2 + nx3*ny4*nz1_dy*s2 + nx3*ny4*nz1*s2_dy 
+ nx4_dy*ny1*nz3*s2 + nx4*ny1_dy*nz3*s2 + nx4*ny1*nz3_dy*s2 + nx4*ny1*nz3*s2_dy 
- nx1_dy*ny4*nz3*s2 - nx1*ny4_dy*nz3*s2 - nx1*ny4*nz3_dy*s2 - nx1*ny4*nz3*s2_dy 
- nx3_dy*ny1*nz4*s2 - nx3*ny1_dy*nz4*s2 - nx3*ny1*nz4_dy*s2 - nx3*ny1*nz4*s2_dy 
+ nx1_dy*ny3*nz4*s2 + nx1*ny3_dy*nz4*s2 + nx1*ny3*nz4_dy*s2 + nx1*ny3*nz4*s2_dy 
+ nx4_dy*ny2*nz1*s3 + nx4*ny2_dy*nz1*s3 + nx4*ny2*nz1_dy*s3 + nx4*ny2*nz1*s3_dy 
- nx2_dy*ny4*nz1*s3 - nx2*ny4_dy*nz1*s3 - nx2*ny4*nz1_dy*s3 - nx2*ny4*nz1*s3_dy 
- nx4_dy*ny1*nz2*s3 - nx4*ny1_dy*nz2*s3 - nx4*ny1*nz2_dy*s3 - nx4*ny1*nz2*s3_dy 
+ nx1_dy*ny4*nz2*s3 + nx1*ny4_dy*nz2*s3 + nx1*ny4*nz2_dy*s3 + nx1*ny4*nz2*s3_dy 
+ nx2_dy*ny1*nz4*s3 + nx2*ny1_dy*nz4*s3 + nx2*ny1*nz4_dy*s3 + nx2*ny1*nz4*s3_dy 
- nx1_dy*ny2*nz4*s3 - nx1*ny2_dy*nz4*s3 - nx1*ny2*nz4_dy*s3 - nx1*ny2*nz4*s3_dy 
- nx3_dy*ny2*nz1*s4 - nx3*ny2_dy*nz1*s4 - nx3*ny2*nz1_dy*s4 - nx3*ny2*nz1*s4_dy 
+ nx2_dy*ny3*nz1*s4 + nx2*ny3_dy*nz1*s4 + nx2*ny3*nz1_dy*s4 + nx2*ny3*nz1*s4_dy 
+ nx3_dy*ny1*nz2*s4 + nx3*ny1_dy*nz2*s4 + nx3*ny1*nz2_dy*s4 + nx3*ny1*nz2*s4_dy 
- nx1_dy*ny3*nz2*s4 - nx1*ny3_dy*nz2*s4 - nx1*ny3*nz2_dy*s4 - nx1*ny3*nz2*s4_dy 
- nx2_dy*ny1*nz3*s4 - nx2*ny1_dy*nz3*s4 - nx2*ny1*nz3_dy*s4 - nx2*ny1*nz3*s4_dy 
+ nx1_dy*ny2*nz3*s4 + nx1*ny2_dy*nz3*s4 + nx1*ny2*nz3_dy*s4 + nx1*ny2*nz3*s4_dy;

         xins_dz = 
  nx4_dz*ny3*nz2*s1 + nx4*ny3_dz*nz2*s1 + nx4*ny3*nz2_dz*s1 + nx4*ny3*nz2*s1_dz 
- nx3_dz*ny4*nz2*s1 - nx3*ny4_dz*nz2*s1 - nx3*ny4*nz2_dz*s1 - nx3*ny4*nz2*s1_dz 
- nx4_dz*ny2*nz3*s1 - nx4*ny2_dz*nz3*s1 - nx4*ny2*nz3_dz*s1 - nx4*ny2*nz3*s1_dz 
+ nx2_dz*ny4*nz3*s1 + nx2*ny4_dz*nz3*s1 + nx2*ny4*nz3_dz*s1 + nx2*ny4*nz3*s1_dz 
+ nx3_dz*ny2*nz4*s1 + nx3*ny2_dz*nz4*s1 + nx3*ny2*nz4_dz*s1 + nx3*ny2*nz4*s1_dz 
- nx2_dz*ny3*nz4*s1 - nx2*ny3_dz*nz4*s1 - nx2*ny3*nz4_dz*s1 - nx2*ny3*nz4*s1_dz 
- nx4_dz*ny3*nz1*s2 - nx4*ny3_dz*nz1*s2 - nx4*ny3*nz1_dz*s2 - nx4*ny3*nz1*s2_dz 
+ nx3_dz*ny4*nz1*s2 + nx3*ny4_dz*nz1*s2 + nx3*ny4*nz1_dz*s2 + nx3*ny4*nz1*s2_dz 
+ nx4_dz*ny1*nz3*s2 + nx4*ny1_dz*nz3*s2 + nx4*ny1*nz3_dz*s2 + nx4*ny1*nz3*s2_dz 
- nx1_dz*ny4*nz3*s2 - nx1*ny4_dz*nz3*s2 - nx1*ny4*nz3_dz*s2 - nx1*ny4*nz3*s2_dz 
- nx3_dz*ny1*nz4*s2 - nx3*ny1_dz*nz4*s2 - nx3*ny1*nz4_dz*s2 - nx3*ny1*nz4*s2_dz 
+ nx1_dz*ny3*nz4*s2 + nx1*ny3_dz*nz4*s2 + nx1*ny3*nz4_dz*s2 + nx1*ny3*nz4*s2_dz 
+ nx4_dz*ny2*nz1*s3 + nx4*ny2_dz*nz1*s3 + nx4*ny2*nz1_dz*s3 + nx4*ny2*nz1*s3_dz 
- nx2_dz*ny4*nz1*s3 - nx2*ny4_dz*nz1*s3 - nx2*ny4*nz1_dz*s3 - nx2*ny4*nz1*s3_dz 
- nx4_dz*ny1*nz2*s3 - nx4*ny1_dz*nz2*s3 - nx4*ny1*nz2_dz*s3 - nx4*ny1*nz2*s3_dz 
+ nx1_dz*ny4*nz2*s3 + nx1*ny4_dz*nz2*s3 + nx1*ny4*nz2_dz*s3 + nx1*ny4*nz2*s3_dz 
+ nx2_dz*ny1*nz4*s3 + nx2*ny1_dz*nz4*s3 + nx2*ny1*nz4_dz*s3 + nx2*ny1*nz4*s3_dz 
- nx1_dz*ny2*nz4*s3 - nx1*ny2_dz*nz4*s3 - nx1*ny2*nz4_dz*s3 - nx1*ny2*nz4*s3_dz 
- nx3_dz*ny2*nz1*s4 - nx3*ny2_dz*nz1*s4 - nx3*ny2*nz1_dz*s4 - nx3*ny2*nz1*s4_dz 
+ nx2_dz*ny3*nz1*s4 + nx2*ny3_dz*nz1*s4 + nx2*ny3*nz1_dz*s4 + nx2*ny3*nz1*s4_dz 
+ nx3_dz*ny1*nz2*s4 + nx3*ny1_dz*nz2*s4 + nx3*ny1*nz2_dz*s4 + nx3*ny1*nz2*s4_dz 
- nx1_dz*ny3*nz2*s4 - nx1*ny3_dz*nz2*s4 - nx1*ny3*nz2_dz*s4 - nx1*ny3*nz2*s4_dz 
- nx2_dz*ny1*nz3*s4 - nx2*ny1_dz*nz3*s4 - nx2*ny1*nz3_dz*s4 - nx2*ny1*nz3*s4_dz 
+ nx1_dz*ny2*nz3*s4 + nx1*ny2_dz*nz3*s4 + nx1*ny2*nz3_dz*s4 + nx1*ny2*nz3*s4_dz;

	 xins_dx = ( det * xins_dx - xins * det_dx ) / det / det ;

	 xins_dy = ( det * xins_dy - xins * det_dy ) / det / det ;

	 xins_dz = ( det * xins_dz - xins * det_dz ) / det / det ;

	 xins = xins / det;

	 dARdx[0] = ( circ * xins_dx - xins * circ_dx ) / circ / circ * 3.0;
	 dARdx[1] = ( circ * xins_dy - xins * circ_dy ) / circ / circ * 3.0;
	 dARdx[2] = ( circ * xins_dz - xins * circ_dz ) / circ / circ * 3.0;

	 *ar = xins/circ*3.0;

}

double gridMinVolume( Grid *grid )
{
  double min_volume;
  int total_count;
  gridMinVolumeAndCount( grid, &min_volume, &total_count );
  return min_volume;
}

Grid *gridMinVolumeAndCount( Grid *grid, double *min_volume, int *total_count )
{
  int cellId, nodes[4];
  double volume, minVol;
  int count;

  minVol = 999.0;
  count = 0;
  for (cellId=0;cellId<gridMaxCell(grid);cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) ){
      volume = gridVolume(grid, nodes);
      if (volume < GRID_VOLUME_TOL ) {
	count++;
      }
      minVol = MIN(minVol, volume);
    }
  *min_volume = minVol;
  *total_count = count;
  return grid;
}

GridBool gridNegCellAroundNode( Grid *grid, int node )
{
  int nodes[4];
  AdjIterator it;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ) {
    gridCell( grid, adjItem(it), nodes );
    if (gridVolume(grid, nodes) <= 0.0) return TRUE;
  }

  return FALSE;
}

GridBool gridNegCellAroundNodeExceptGem( Grid *grid, int node )
{
  int igem, cellId, nodes[4];
  GridBool inGem;
  AdjIterator it;

  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ) {
    cellId = adjItem(it);
    inGem = FALSE;
    for ( igem =0; !inGem && igem < gridNGem(grid) ; igem++)
      inGem = inGem || (cellId == gridGem(grid,igem));
    if ( !inGem ) {
      gridCell( grid, cellId, nodes );
      if (gridVolume(grid, nodes) <= 0.0) return TRUE;
    }
  }

  return FALSE;
}

double gridMinARAroundNodeExceptGem( Grid *grid, int node )
{
  int igem, cellId, nodes[4];
  double minAR;
  GridBool inGem;
  AdjIterator it;

  minAR = 2.0;

  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ) {
    cellId = adjItem(it);
    inGem = FALSE;
    for ( igem =0; !inGem && igem < gridNGem(grid) ; igem++)
      inGem = inGem || (cellId == gridGem(grid,igem));
    if ( !inGem ) {
      gridCell( grid, cellId, nodes );
      minAR = MIN(minAR,gridAR(grid, nodes));
    }
  }

  return minAR;
}

double gridMinARAroundNodeExceptGemRecon( Grid *grid, int node, int becomes )
{
  int igem, i, cellId, nodes[4];
  double minAR;
  GridBool inGem;
  AdjIterator it;

  minAR = 2.0;

  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ) {
    cellId = adjItem(it);
    inGem = FALSE;
    for ( igem =0; !inGem && igem < gridNGem(grid) ; igem++)
      inGem = inGem || (cellId == gridGem(grid,igem));
    if ( !inGem ) {
      gridCell( grid, cellId, nodes );
      for (i=0;i<4;i++) if (nodes[i]==node) nodes[i]=becomes;
      minAR = MIN(minAR,gridAR(grid, nodes));
    }
  }

  return minAR;
}

double gridMinAR( Grid *grid )
{
  int cellId, nodes[4];
  double minAR;
  minAR = 999.0;
  for (cellId=0;cellId<gridMaxCell(grid);cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      minAR = MIN(minAR, gridAR(grid, nodes) );
  return minAR;
}

double gridMinThawedAR( Grid *grid )
{
  int cellId, nodes[4];
  double minAR;
  minAR = 999.0;
  for (cellId=0;cellId<gridMaxCell(grid);cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) &&
	 ( !gridNodeFrozen(grid,nodes[0]) ||
	   !gridNodeFrozen(grid,nodes[1]) ||
	   !gridNodeFrozen(grid,nodes[2]) ||
	   !gridNodeFrozen(grid,nodes[3]) )  )
      minAR = MIN(minAR, gridAR(grid, nodes) );
  return minAR;
}

GridBool gridRightHandedFace(Grid *grid, int face ){
  int cell;
  int faceNodes[4], faceId;
  int cellNodes[4];

  if ( grid != gridFace(grid, face, faceNodes, &faceId ) ) return FALSE;

  cell = gridFindCellWithFace(grid, face );
  if (grid != gridCell(grid, cell, cellNodes) ) {
    printf("gridRightHandedFace: %s: %d: no cell - face %d id %d n %d %d %d\n",
	   __FILE__, __LINE__,
	   face, faceId, faceNodes[0], faceNodes[1], faceNodes[2]);
    return FALSE;
  }

  faceNodes[3] = cellNodes[0] + cellNodes[1] + cellNodes[2] + cellNodes[3]
      - faceNodes[0] - faceNodes[1] - faceNodes[2];

  return (gridVolume(grid, faceNodes) > 0.0);
}

GridBool gridRightHandedBoundary( Grid *grid )
{
  int face, nodes[3], faceId;
  GridBool rightHanded;
  rightHanded = TRUE;

  for (face=0;face<gridMaxFace(grid);face++)
    if ( grid == gridFace(grid, face, nodes, &faceId ) )
      rightHanded = gridRightHandedFace(grid, face) && rightHanded;

  return rightHanded;
}

Grid *gridFlipLeftHandedFaces( Grid *grid )
{
  int face;

  for (face=0;face<gridMaxFace(grid);face++)
    if ( !gridRightHandedFace(grid, face) )
      gridFlipFace(grid, face );

  return grid;
}

double gridFaceArea(Grid *grid, int n0, int n1, int n2 )
{
  double xyz0[3], xyz1[3], xyz2[3];
  double edge1[3], edge2[3], norm[3], length; 
  
  if (grid != gridNodeXYZ(grid,n0,xyz0) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n1,xyz1) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n2,xyz2) ) return -1.0;

  gridSubtractVector(xyz1,xyz0,edge1);
  gridSubtractVector(xyz2,xyz0,edge2);

  gridCrossProduct(edge1, edge2, norm);

  length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  return  0.5*length;
}

GridBool gridRightHandedBoundaryUV( Grid *grid )
{
  int face, nodes[3], faceId;
  GridBool rightHanded;
  rightHanded = TRUE;
  
  for (face=0;face<gridMaxFace(grid);face++)
    if ( grid == gridFace(grid, face, nodes, &faceId ) )
      rightHanded = (gridFaceAreaUV(grid, face)>GRID_AREA_TOL) && rightHanded;
  
  return rightHanded;
}

double gridFaceAreaUV(Grid *grid, int face)
{
  int faceId, nodes[3];
  double uv0[2], uv1[2], uv2[2];

  if (grid != gridFace(grid,face,nodes,&faceId)) {
    printf("%s: %d: Invalid Face number %d\n",__FILE__,__LINE__,face);
    return DBL_MAX;  
  }
  gridNodeUV( grid, nodes[0], faceId, uv0);
  gridNodeUV( grid, nodes[1], faceId, uv1);
  gridNodeUV( grid, nodes[2], faceId, uv2);

  if ( DBL_MAX == uv0[0] ) return DBL_MAX;
  if ( DBL_MAX == uv0[1] ) return DBL_MAX;
  if ( DBL_MAX == uv1[0] ) return DBL_MAX;
  if ( DBL_MAX == uv1[1] ) return DBL_MAX;
  if ( DBL_MAX == uv2[0] ) return DBL_MAX;
  if ( DBL_MAX == uv2[1] ) return DBL_MAX;

  return gridFaceAreaUVDirect(uv0, uv1, uv2, faceId);
}

double gridFaceAreaUVDirect(double *uv0,  double *uv1,  double *uv2,
			    int faceId)
{
  double edge0[3], edge1[3], norm[3];
  int vol = 1;
  double area;

  /* Using 3D Vector math to save typing */
  edge0[2] = edge1[2] = 0.0;

  edge0[0] = uv1[0]-uv0[0];
  edge0[1] = uv1[1]-uv0[1];

  edge1[0] = uv2[0]-uv0[0];
  edge1[1] = uv2[1]-uv0[1];

  gridCrossProduct(edge0,edge1,norm);

  /* CAPrI Normal points outside, triangles oriented inside wrt Face */
  if( !CADGeom_ReversedSurfaceNormal(vol,faceId) )
    area = (-0.5*norm[2]);
  else
    area = (0.5*norm[2]);

  return area;
}

Grid *gridMinFaceAreaUV(Grid *grid, int node, double *min_area)
{
  AdjIterator it;
  int face;
  double local_area;

  *min_area = DBL_MAX;

  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it); 
	it = adjNext(it) ){
    face = adjItem(it);
    local_area = gridFaceAreaUV(grid, face);
    if ( local_area < *min_area ) *min_area = local_area;
  }

  return grid;
}

double gridMinCellFaceAreaUV(Grid *grid, int *nodes)
{
  int face;
  double local_area, min_area;

  min_area = DBL_MAX;

  face = gridFindFace( grid, nodes[1], nodes[3], nodes[2] );
  if (EMPTY != face){
    local_area = gridFaceAreaUV(grid, face);
    if ( local_area < min_area ) min_area = local_area;
  }
  face = gridFindFace( grid, nodes[0], nodes[2], nodes[3] );
  if (EMPTY != face){
    local_area = gridFaceAreaUV(grid, face);
    if ( local_area < min_area ) min_area = local_area;
  }
  face = gridFindFace( grid, nodes[0], nodes[3], nodes[1] );
  if (EMPTY != face){
    local_area = gridFaceAreaUV(grid, face);
    if ( local_area < min_area ) min_area = local_area;
  }
  face = gridFindFace( grid, nodes[0], nodes[1], nodes[2] );
  if (EMPTY != face){
    local_area = gridFaceAreaUV(grid, face);
    if ( local_area < min_area ) min_area = local_area;
  }

  return min_area;
}

double gridMinGridFaceAreaUV(Grid *grid)
{
  int face, nodes[3], faceId;
  double min_area;

  min_area = DBL_MAX;

  for ( face = 0 ; face < gridMaxFace(grid) ; face++ ) {
    if (grid == gridFace(grid, face, nodes, &faceId ) ) {
      min_area = MIN(min_area, gridFaceAreaUV(grid, face) );
    }
  }

  return min_area;
}

double gridFaceAR(Grid *grid, int n0, int n1, int n2 )
{
  double xyz0[3], xyz1[3], xyz2[3];
  double e0[3], e1[3], e2[3];
  double l0, l1, l2;
  double n[3], a28, p, ar;

  if (grid != gridNodeXYZ(grid,n0,xyz0) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n1,xyz1) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n2,xyz2) ) return -1.0;

  gridSubtractVector(xyz1,xyz0,e0);
  gridSubtractVector(xyz2,xyz1,e1);
  gridSubtractVector(xyz0,xyz2,e2);

  gridCrossProduct(e1, e2, n);

  l0 = sqrt(gridDotProduct(e0,e0));
  l1 = sqrt(gridDotProduct(e1,e1));
  l2 = sqrt(gridDotProduct(e2,e2));

  a28 = 4.0*gridDotProduct(n,n);

  p = l0+l1+l2;

  ar = a28/(p*l0*l1*l2);

  return ar;

}

double gridMinFaceMR( Grid *grid )
{
  int face, nodes[3], faceId;
  double minMR;
  minMR = 999.0;
  for (face=0;face<gridMaxFace(grid);face++) {
    if ( grid == gridFace(grid,face,nodes,&faceId) ) {
      minMR = MIN(minMR, gridFaceMR(grid, nodes[0], nodes[1], nodes[2]) );
    }
  }
  return minMR;
}

double gridMinThawedFaceMR( Grid *grid )
{
  int face, nodes[3], faceId;
  double minMR;
  minMR = 999.0;
  for (face=0;face<gridMaxFace(grid);face++) {
    if ( grid == gridFace(grid,face,nodes,&faceId) &&
	 ( !gridNodeFrozen(grid,nodes[0]) ||
	   !gridNodeFrozen(grid,nodes[1]) ||
	   !gridNodeFrozen(grid,nodes[2]) ) ) {
      minMR = MIN(minMR, gridFaceMR(grid, nodes[0], nodes[1], nodes[2]) );
    }
  }
  return minMR;
}

#define SQRT3 (1.73205080756888)

double gridFaceMR(Grid *grid, int n0, int n1, int n2 )
{
  double xyz1[3], xyz2[3], xyz3[3]; 
  double x1, x2, x3; 
  double y1, y2, y3; 
  double z1, z2, z3;
  int i;
  double m[6], m0[6], m1[6], m2[6], j[9];
  
  double ex1, ey1, ez1;
  double ex2, ey2, ez2;
  double ex3, ey3, ez3;
  double nx, ny, nz;
  double l1, l2, l3;
  double a,d;
  double r;

  if (grid != gridNodeXYZ(grid,n0,xyz1) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n1,xyz2) ) return -1.0;
  if (grid != gridNodeXYZ(grid,n2,xyz3) ) return -1.0;

  x1 = xyz1[0];
  y1 = xyz1[1];
  z1 = xyz1[2];

  x2 = xyz2[0];
  y2 = xyz2[1];
  z2 = xyz2[2];

  x3 = xyz3[0];
  y3 = xyz3[1];
  z3 = xyz3[2];
 
  gridMap(grid,n0,m0);
  gridMap(grid,n1,m1);
  gridMap(grid,n2,m2);

  for (i=0;i<6;i++) m[i]=0.333333333333333*(m0[i]+m1[i]+m2[i]);
  if (grid != gridConvertMetricToJacobian(grid, m, j) ) {
    printf("%s: %d: gridFaceMR: gridConvertMetricToJacobian NULL\n",
	   __FILE__,__LINE__);
    printf("m \n %25.18f %25.18f %25.18f \n %25.18f %25.18f %25.18f \n %25.18f %25.18f %25.18f\n",m[0],m[1],m[2],m[1],m[3],m[4],m[2],m[4],m[5]);
  }

  gridMapXYZWithJ(j, &x1, &y1, &z1);
  gridMapXYZWithJ(j, &x2, &y2, &z2);
  gridMapXYZWithJ(j, &x3, &y3, &z3);

  ex1 = x2-x1;
  ey1 = y2-y1;
  ez1 = z2-z1;

  ex2 = x3-x2;
  ey2 = y3-y2;
  ez2 = z3-z2;

  ex3 = x1-x3;
  ey3 = y1-y3;
  ez3 = z1-z3;

  nx = ey1*ez2 - ez1*ey2; 

  ny = ez1*ex2 - ex1*ez2; 

  nz = ex1*ey2 - ey1*ex2; 

  l1 = ex1*ex1 + ey1*ey1 + ez1*ez1;

  l2 = ex2*ex2 + ey2*ey2 + ez2*ez2;

  l3 = ex3*ex3 + ey3*ey3 + ez3*ez3;

  a = sqrt(nx*nx + ny*ny + nz*nz);

  a = 0.5*a;

  d = l1+l2+l3;

  r = a/d;

  return 4*SQRT3*r;

}

Grid *gridFaceMRDerivative(Grid *grid, int *nodes, double *mr, double *dMRdx )
{
  double xyz1[3], xyz2[3], xyz3[3]; 
  double x1, x2, x3; 
  double y1, y2, y3; 
  double z1, z2, z3;
  int i;
  double m[6], m0[6], m1[6], m2[6], j[9];
  
  if (grid != gridNodeXYZ(grid,nodes[0],xyz1) ) return NULL;
  if (grid != gridNodeXYZ(grid,nodes[1],xyz2) ) return NULL;
  if (grid != gridNodeXYZ(grid,nodes[2],xyz3) ) return NULL;

  x1 = xyz1[0];
  y1 = xyz1[1];
  z1 = xyz1[2];

  x2 = xyz2[0];
  y2 = xyz2[1];
  z2 = xyz2[2];

  x3 = xyz3[0];
  y3 = xyz3[1];
  z3 = xyz3[2];
 
  gridMap(grid,nodes[0],m0);
  gridMap(grid,nodes[1],m1);
  gridMap(grid,nodes[2],m2);

  for (i=0;i<6;i++) m[i]=0.333333333333333*(m0[i]+m1[i]+m2[i]);
  if (grid != gridConvertMetricToJacobian(grid, m, j) ) {
    printf("%s: %d: griFaceMRDerivative: gridConvertMetricToJacobian NULL\n",
	   __FILE__,__LINE__);
  }
  gridConvertMetricToJacobian( grid, m, j );

  gridMapXYZWithJ(j, &x1, &y1, &z1);
  gridMapXYZWithJ(j, &x2, &y2, &z2);
  gridMapXYZWithJ(j, &x3, &y3, &z3);

  FaceMRDerivative(x1,y1,z1,x2,y2,z2,x3,y3,z3,mr,dMRdx );
  return grid;
}

void FaceMRDerivative(double x1, double y1, double z1,
		      double x2, double y2, double z2,
		      double x3, double y3, double z3,
		      double *mr, double *dMRdx  )
{

  double ex1, ey1, ez1;
  double ex2, ey2, ez2;
  double ex3, ey3, ez3;

  double ex1_dx, ey1_dy, ez1_dz;
  double ex3_dx, ey3_dy, ez3_dz;

  double nx, ny, nz;
  double nx_dy,nx_dz;
  double ny_dx,ny_dz;
  double nz_dx,nz_dy;
  double l1, l2, l3;
  double l1_dx, l1_dy, l1_dz;
  double l3_dx, l3_dy, l3_dz;
  double a,d;
  double a_dx, a_dy, a_dz;
  double d_dx, d_dy, d_dz;
  double r;
  double r_dx, r_dy, r_dz;

  ex1 = x2-x1;
  ey1 = y2-y1;
  ez1 = z2-z1;

  ex1_dx = -1.0;
  ey1_dy = -1.0;
  ez1_dz = -1.0;

  ex2 = x3-x2;
  ey2 = y3-y2;
  ez2 = z3-z2;

  ex3 = x1-x3;
  ey3 = y1-y3;
  ez3 = z1-z3;

  ex3_dx = 1.0;
  ey3_dy = 1.0;
  ez3_dz = 1.0;

  nx = ey1*ez2 - ez1*ey2; 
  nx_dy =  ey1_dy * ez2;
  nx_dz = -ez1_dz*ey2; 

  ny = ez1*ex2 - ex1*ez2; 
  ny_dx = -ex1_dx*ez2; 
  ny_dz =  ez1_dz*ex2; 

  nz = ex1*ey2 - ey1*ex2; 
  nz_dx =  ex1_dx*ey2; 
  nz_dy = -ey1_dy*ex2; 

  l1 = ex1*ex1 + ey1*ey1 + ez1*ez1;
  l1_dx = 2*ex1_dx*ex1;
  l1_dy = 2*ey1_dy*ey1;
  l1_dz = 2*ez1_dz*ez1;

  l2 = ex2*ex2 + ey2*ey2 + ez2*ez2;

  l3 = ex3*ex3 + ey3*ey3 + ez3*ez3;
  l3_dx = 2*ex3_dx*ex3;
  l3_dy = 2*ey3_dy*ey3;
  l3_dz = 2*ez3_dz*ez3;

  a = sqrt(nx*nx + ny*ny + nz*nz);
  a_dx = (ny*ny_dx + nz*nz_dx)/a;
  a_dy = (nx*nx_dy + nz*nz_dy)/a;
  a_dz = (nx*nx_dz + ny*ny_dz)/a; 

  a = 0.5*a;
  a_dx = 0.5*a_dx;
  a_dy = 0.5*a_dy;
  a_dz = 0.5*a_dz;

  d = l1+l2+l3;
  d_dx = l1_dx + l3_dx;
  d_dy = l1_dy + l3_dy;
  d_dz = l1_dz + l3_dz;

  r = a/d;
  r_dx = ( d * a_dx - a * d_dx ) / d / d;
  r_dy = ( d * a_dy - a * d_dy ) / d / d;
  r_dz = ( d * a_dz - a * d_dz ) / d / d;

  *mr = 4*SQRT3*r;

  dMRdx[0] = 4*SQRT3*r_dx;
  dMRdx[1] = 4*SQRT3*r_dy;
  dMRdx[2] = 4*SQRT3*r_dz;

}

Grid *gridNodeFaceMR(Grid *grid, int node, double *mr )
{
  AdjIterator it;
  int nodes[3], faceId;
  double local_mr;

  *mr = 1.0;

  if ( gridNegCellAroundNode( grid, node ) ) {
    *mr = -1.0;
    return grid;
  }

  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it); 
	it = adjNext(it) ){
    gridFace(grid,adjItem(it),nodes,&faceId);
    local_mr = gridFaceMR(grid, nodes[0], nodes[1], nodes[2]);
    if ( local_mr < *mr ) *mr = local_mr;
  }

  return grid;
}

Grid *gridNodeFaceMRDerivative (Grid *grid, int node, double *mr, double *dMRdx )
{
  AdjIterator it;
  int originalNodes[3], orientedNodes[3], faceId, index, inode;
  double local_mr, local_dMRdx[3];

  *mr = 1.1;
  dMRdx[0] = DBL_MAX;
  dMRdx[1] = DBL_MAX;
  dMRdx[2] = DBL_MAX;

  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    gridFace(grid,adjItem(it),originalNodes,&faceId);
    orientedNodes[0] = node;
    index = 1;
    for (inode=0; inode<3;inode++) {
      if ( node != originalNodes[inode]){
	orientedNodes[index] = originalNodes[inode];
	index++;
      }
    }

    if ( grid != gridFaceMRDerivative(grid, orientedNodes, 
				      &local_mr, local_dMRdx ) ) {
      *mr = 0.0;
      dMRdx[0] = DBL_MAX;
      dMRdx[1] = DBL_MAX;
      dMRdx[2] = DBL_MAX;
      return NULL;
    }
    if ( local_mr < *mr ) {
      *mr = local_mr;
      dMRdx[0] = local_dMRdx[0];
      dMRdx[1] = local_dMRdx[1];
      dMRdx[2] = local_dMRdx[2];
    }
  }

  return grid;
}

double gridCellMeanRatio( double *xyz0, double *xyz1, double *xyz2, double *xyz3 )
{
  double edge1[3], edge2[3], edge3[3];
  double edge4[3], edge5[3], edge6[3];
  double norm[3], volume, mr; 

  gridSubtractVector( xyz1, xyz0, edge1);
  gridSubtractVector( xyz2, xyz0, edge2);
  gridSubtractVector( xyz3, xyz0, edge3);
  gridSubtractVector( xyz2, xyz1, edge4);
  gridSubtractVector( xyz3, xyz1, edge5);
  gridSubtractVector( xyz3, xyz2, edge6);

  gridCrossProduct( edge5, edge4, norm );

  volume = -gridDotProduct(norm,edge1)/6.0;

  if ( volume < GRID_VOLUME_TOL ) return -1.0;

  mr = 12.0 * pow(9.0*volume*volume,1.0/3.0) /
    ( gridDotProduct(edge1,edge1) +
      gridDotProduct(edge2,edge2) + 
      gridDotProduct(edge3,edge3) +
      gridDotProduct(edge4,edge4) +
      gridDotProduct(edge5,edge5) +
      gridDotProduct(edge6,edge6) );

  return mr;
}

void gridCellMeanRatioDerivative( double *xyz0, double *xyz1, 
				  double *xyz2, double *xyz3,
				  double *mr, double *dMRdx)
{
  double edge1[3], edge2[3], edge3[3];
  double edge4[3], edge5[3], edge6[3];
  double norm[3];
  double volume, dVdx, dVdy, dVdz;
  double cuberoot;
  double i, didx, didy, didz;
  double n, dndx, dndy, dndz;
  double d, dddx, dddy, dddz;
  double d2;

  gridSubtractVector( xyz1, xyz0, edge1);
  gridSubtractVector( xyz2, xyz0, edge2);
  gridSubtractVector( xyz3, xyz0, edge3);
  gridSubtractVector( xyz2, xyz1, edge4);
  gridSubtractVector( xyz3, xyz1, edge5);
  gridSubtractVector( xyz3, xyz2, edge6);

  gridCrossProduct( edge5, edge4, norm );

  volume = -gridDotProduct(norm,edge1)/6.0;

  dVdx = norm[0]/6.0;
  dVdy = norm[1]/6.0;
  dVdz = norm[2]/6.0;

  i = 9.0*volume*volume;

  didx = 9.0 * 2.0 * volume * dVdx;
  didy = 9.0 * 2.0 * volume * dVdy;
  didz = 9.0 * 2.0 * volume * dVdz;

  n = 12.0 * pow(i,1.0/3.0);

  cuberoot = 4.0 * pow(i,-2.0/3.0);

  dndx = cuberoot * didx;
  dndy = cuberoot * didy;
  dndz = cuberoot * didz;

  d = gridDotProduct(edge1,edge1)
    + gridDotProduct(edge2,edge2)
    + gridDotProduct(edge3,edge3)  
    + gridDotProduct(edge4,edge4)
    + gridDotProduct(edge5,edge5)
    + gridDotProduct(edge6,edge6) ;

  dddx = -2.0 * ( edge1[0] + edge2[0] + edge3[0]);
  dddy = -2.0 * ( edge1[1] + edge2[1] + edge3[1]);
  dddz = -2.0 * ( edge1[2] + edge2[2] + edge3[2]);

  *mr = n/d;

  d2 = 1.0/(d*d);

  dMRdx[0] = (d*dndx-n*dddx)*d2;
  dMRdx[1] = (d*dndy-n*dddy)*d2;
  dMRdx[2] = (d*dndz-n*dddz)*d2;
}

Grid *gridCollapseCost(Grid *grid, int node0, int node1, double *currentCost, 
		       double *node0Cost, double *node1Cost)
{
  double xyz0[3], xyz1[3];
  double map0[6], map1[6];
  int faceId0, faceId1;
  double node0Id0uv[2], node1Id0uv[2];
  double node0Id1uv[2], node1Id1uv[2];
  int edge, edgeId;
  double t0, t1;
  int gap0, gap1;
  GridBool onBoundary;

  *currentCost = 0.0;
  *node0Cost = 0.0;
  *node1Cost = 0.0;

  if ( NULL == gridNodeXYZ( grid, node0, xyz0) ) return NULL;
  if ( NULL == gridNodeXYZ( grid, node1, xyz1) ) return NULL;
  if ( NULL == gridMap( grid, node0, map0) ) return NULL;
  if ( NULL == gridMap( grid, node1, map1) ) return NULL;

  if (grid!=gridEquator(grid,node0,node1)) return NULL;
  onBoundary = !gridContinuousEquator(grid);

  *currentCost = MIN( gridMinARAroundNodeExceptGem(grid,node0),
		      gridMinARAroundNodeExceptGem(grid,node1) );

  faceId0 = faceId1 = EMPTY;
  node0Id0uv[0] = node0Id0uv[1] = DBL_MAX;
  node1Id0uv[0] = node1Id0uv[1] = DBL_MAX;
  node0Id1uv[0] = node0Id1uv[1] = DBL_MAX;
  node1Id1uv[0] = node1Id1uv[1] = DBL_MAX;
  edge = edgeId = EMPTY;
  t0 = t1 = DBL_MAX;
  if ( onBoundary ) {
    gap0 = gridEqu(grid,0);
    gap1 = gridEqu(grid,gridNGem(grid));
    faceId0 = gridFaceId(grid, node0, node1, gap0 );
    faceId1 = gridFaceId(grid, node0, node1, gap1 );
    if ( faceId0 == EMPTY || faceId1 == EMPTY ) {
      printf("%s: %d:empty gap faces\n",__FILE__,__LINE__);
      return NULL;
    }
    gridNodeUV(grid,node0,faceId0,node0Id0uv);
    gridNodeUV(grid,node1,faceId0,node1Id0uv);
    gridNodeUV(grid,node0,faceId1,node0Id1uv);
    gridNodeUV(grid,node1,faceId1,node1Id1uv);
    edge = gridFindEdge(grid,node0,node1);
    if ( edge != EMPTY ) {
      edgeId = gridEdgeId(grid,node0,node1);
      gridNodeT(grid, node0, edgeId, &t0 );
      gridNodeT(grid, node1, edgeId, &t1 );
    }
  }

  gridSetNodeXYZ( grid, node1, xyz0);
  gridSetMap( grid, node1, map0[0],map0[1],map0[2],map0[3],map0[4],map0[5]);
  if ( EMPTY != faceId0)
    gridSetNodeUV( grid, node1, faceId0, node0Id0uv[0],  node0Id0uv[1] );
  if ( EMPTY != faceId1)
    gridSetNodeUV( grid, node1, faceId1, node0Id1uv[0],  node0Id1uv[1] );
  if ( edgeId != EMPTY ) gridSetNodeT(grid, node1, edgeId, t0 );

  *node0Cost = MIN( gridMinARAroundNodeExceptGem(grid,node0),
		    gridMinARAroundNodeExceptGem(grid,node1) );

  gridSetNodeXYZ( grid, node1, xyz1);
  gridSetMap( grid, node1, map1[0],map1[1],map1[2],map1[3],map1[4],map1[5]);
  if ( EMPTY != faceId0)
    gridSetNodeUV( grid, node1, faceId0, node1Id0uv[0],  node1Id0uv[1] );
  if ( EMPTY != faceId1)
    gridSetNodeUV( grid, node1, faceId1, node1Id1uv[0],  node1Id1uv[1] );
  if ( edgeId != EMPTY ) gridSetNodeT(grid, node1, edgeId, t1 );


  gridSetNodeXYZ( grid, node0, xyz1);
  gridSetMap( grid, node0, map1[0],map1[1],map1[2],map1[3],map1[4],map1[5]);
  if ( EMPTY != faceId0)
    gridSetNodeUV( grid, node0, faceId0, node1Id0uv[0],  node1Id0uv[1] );
  if ( EMPTY != faceId1)
    gridSetNodeUV( grid, node0, faceId1, node1Id1uv[0],  node1Id1uv[1] );
  if ( edgeId != EMPTY ) gridSetNodeT(grid, node0, edgeId, t1 );

  *node1Cost = MIN( gridMinARAroundNodeExceptGem(grid,node0),
		    gridMinARAroundNodeExceptGem(grid,node1) );

  gridSetNodeXYZ( grid, node0, xyz0);
  gridSetMap( grid, node0, map0[0],map0[1],map0[2],map0[3],map0[4],map0[5]);
  if ( EMPTY != faceId0)
    gridSetNodeUV( grid, node0, faceId0, node0Id0uv[0],  node0Id0uv[1] );
  if ( EMPTY != faceId1)
    gridSetNodeUV( grid, node0, faceId1, node0Id1uv[0],  node0Id1uv[1] );
  if ( edgeId != EMPTY ) gridSetNodeT(grid, node0, edgeId, t0 );

  return grid;
}

Grid *gridSplitCost(Grid *grid, int node0, int node1, 
		    double *currentCost, double *splitCost )
{
  double xyz0[3], xyz1[3], newxyz[3];
  int newnode, nodes0[4], nodes1[4];

  double ratio = 0.5;

  int inode, i;

  if ( NULL == gridNodeXYZ( grid, node0, xyz0) ) return NULL;
  if ( NULL == gridNodeXYZ( grid, node1, xyz1) ) return NULL;

  if (grid!=gridMakeGem(grid,node0,node1)) return NULL;
  if (grid!=gridGemAR(grid,currentCost)) return NULL;

  for (inode = 0 ; inode < 3 ; inode++) 
    newxyz[inode] = (1-ratio)*xyz0[inode] + ratio*xyz1[inode]; 

  newnode = gridAddNode(grid, newxyz[0], newxyz[1], newxyz[2] );
  if ( newnode == EMPTY ) return NULL;

  if ( grid != gridInterpolateMap2(grid,node0,node1,ratio,newnode ) ) {
    gridRemoveNode(grid,newnode);
    return NULL;
  }

  (*splitCost) = 10.0;

  for ( i = 0 ; i < gridNGem(grid) ; i++ ){
    gridCell(grid, gridGem(grid,i), nodes0);
    gridCell(grid, gridGem(grid,i), nodes1);
    for (inode = 0 ; inode < 4 ; inode++) {
      if (nodes0[inode] == node0) nodes0[inode]=newnode;
      if (nodes1[inode] == node1) nodes1[inode]=newnode;
    }
    *splitCost = MIN(*splitCost,gridAR( grid, nodes0 ));
    *splitCost = MIN(*splitCost,gridAR( grid, nodes1 ));
  }

  gridRemoveNode(grid,newnode);

  return grid;
}

int parse_key_value_pair( char *line, char *key, int *val );
int parse_key_value_pair( char *line, char *key, int *val )
{
  char *token;
  char *token_start;
  int length;

  token_start = strstr( line, key );
  if ( NULL == token_start )
    {
      printf("no key found\n");
      return 1;
    }
  token_start += strlen( key );
  length = strspn( token_start, "1234567890" );
  if ( 0 == length )
    {
      printf("no int value found\n");
      return 1;
    }
  token = (char *)malloc( (length+1)*sizeof(char) );
  strncpy( token, token_start, length );
  token[length] = '\0';
  *val = atoi(token);
  free( token );
  return 0;
}

Grid *gridSpacingFromTecplot(Grid *grid, char *filename )
{
  FILE *file;
#define LINE_SIZE (1025)
  char line[LINE_SIZE];
  int nnode, ntet;
  int node;
  double x,y,z,spacing;
  file = fopen(filename,"r");
  if (NULL == file) return NULL;

  while ( !feof(file) )
    {
      if (NULL == fgets( line, LINE_SIZE, file ))  return NULL;
      if ( !strncmp(line,"vari",4) || !strncmp(line,"VARI",4) )
	{
	  printf("%s", line);
	}
      if ( !strncmp(line,"zone",4) || !strncmp(line,"ZONE",4) )
	{
	  printf("%s", line);
 	  if( parse_key_value_pair( line, "I=", &nnode ) ) return NULL;
	  if( parse_key_value_pair( line, "J=", &ntet ) ) return NULL;
	  printf("%d nodes %d tets\n", nnode, ntet);
	  if (nnode != gridNNode(grid) || ntet != gridNCell(grid)) 
	    {
	      printf("%s: %d: gridSpacingFromTecplot: sizes nodes %d != %d tets %d != %d \n",
		     __FILE__,__LINE__,
		     nnode,gridNNode(grid),ntet, gridNCell(grid));
	      return NULL;
	    }
	  for (node =0 ; node < nnode ; node++)
	    {
	      if ( 4 != fscanf(file,"%lf %lf %lf %lf",&x,&y,&z,&spacing) )
		{
		  printf("%s: %d: gridSpacingFromTecplot: node %d read failed \n",
			 __FILE__,__LINE__,node);
		  return NULL;
		}
	      if (grid != gridSetSpacing(grid,node,spacing) )
		{
		  printf("%s: %d: gridSpacingFromTecplot: node %d set spc failed \n",
			 __FILE__,__LINE__,node);
		  return NULL;
		}
	    }
	  return grid;
	}
    }

  return grid;
}

