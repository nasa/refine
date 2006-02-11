
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

#include "refine_defs.h"

BEGIN_C_DECLORATION

GridBool intersectAbove( double *vertex0, double *vertex1, double *vertex2,
		     double *node );

GridBool intersectTriangleNode(double *vertex0, double *vertex1, double *vertex2,
			   double *node );

GridBool intersectTriangleSegment(double *vertex0, double *vertex1, double *vertex2,
			      double *node0, double *node1, double *ratio );

GridBool intersectInsideTet(double *vertex0, double *vertex1, 
			double *vertex2, double *vertex3,
			double *node );

GridBool intersectTetSegment(double *vertex0, double *vertex1, 
			 double *vertex2, double *vertex3,
			 double *node0, double *node1 );

GridBool intersectTetTet(double *vertex0, double *vertex1, 
		     double *vertex2, double *vertex3,
		     double *node0, double *node1,
		     double *node2, double *node3 );

END_C_DECLORATION

#endif /* INTERSECT_H */
