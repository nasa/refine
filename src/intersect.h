
/* Compute the intersection of geometry elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef INTERSECT_H
#define INTERSECT_H

#include "master_header.h"

BEGIN_C_DECLORATION

int intersectSide( double *vertex0, double *vertex1, double *vertex2,
                   double *node );

bool intersectTriangleNode(double *vertex0, double *vertex1, double *vertex2,
			   double *node );

bool intersectTriangleSegment(double *vertex0, double *vertex1, double *vertex2,
			      double *node0, double *node1 );

bool intersectTetSegment(double *vertex0, double *vertex1, 
			 double *vertex2, double *vertex3,
			 double *node0, double *node1 );

END_C_DECLORATION

#endif /* INTERSECT_H */
