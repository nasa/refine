
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
