
/* Compute the intersection of geometry elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include "intersect.h"
#include "gridmath.h"

int intersectSide( double *vertex0, double *vertex1, double *vertex2,
                   double *node )
{
  double edge1[3], edge2[3], normal[3], direction[3];
  double height;
  gridSubtractVector(vertex1, vertex0, edge1);
  gridSubtractVector(vertex2, vertex0, edge2);
  gridCrossProduct(edge1, edge2, normal);
  gridSubtractVector(node, vertex0, direction);
  height = gridDotProduct(normal,direction);
  if ( height == 0 ) return 0;
  if ( height > 0 ) return 1;
  return -1;
}

bool intersectTriangleNode( double *vertex0, double *vertex1, double *vertex2,
			    double *node)
{
  int i;
  double side[3], otherSide[3], target[3];
  double triangleNormal[3], targetNormal[3];

  gridSubtractVector(vertex1, vertex0, side);

  gridSubtractVector(vertex2, vertex0, otherSide);
  gridCrossProduct(side, otherSide, triangleNormal);

  gridSubtractVector(node, vertex0, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  gridSubtractVector(vertex2, vertex1, side);
  gridSubtractVector(node, vertex1, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  gridSubtractVector(vertex0, vertex2, side);
  gridSubtractVector(node, vertex2, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  return TRUE;
}

bool intersectTriangleSegment(double *vertex0, double *vertex1, double *vertex2,
			      double *node0, double *node1)
{
  int i;
  double edge1[3], edge2[3], normal[3];
  double dir0[3], dir1[3];
  double h0, h1;
  bool coplanar;
  double denom, intersection[3];

  gridSubtractVector(vertex1, vertex0, edge1);
  gridSubtractVector(vertex2, vertex0, edge2);
  gridCrossProduct(edge1, edge2, normal);
  gridSubtractVector(node0, vertex0, dir0);  
  gridSubtractVector(node1, vertex0, dir1);  

  h0 = gridDotProduct(normal,dir0);
  h1 = gridDotProduct(normal,dir1);

  if (h0>0 && h1>0) return FALSE;
  if (h0<0 && h1<0) return FALSE;

  coplanar = (h0==0 && h1==0);

  if (!coplanar) {
    denom = 1/(h0-h1);
    for(i=0;i<3;i++) intersection[i] = ( h0*dir1[i] - h1*dir0[i] )*denom;
    return intersectTriangleNode(vertex0, vertex1, vertex2, intersection);
  }else{
    if (intersectTriangleNode(vertex0, vertex1, vertex2, node0) ) return TRUE;
    if (intersectTriangleNode(vertex0, vertex1, vertex2, node1) ) return TRUE;
    return FALSE;
  }
}

bool intersectTetSegment(double *vertex0, double *vertex1, 
			 double *vertex2, double *vertex3,
			 double *node0, double *node1)
{

  return TRUE;
}
