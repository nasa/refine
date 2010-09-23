
/* Compute the intersection of geometry elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#ifndef INTERSECT_H
#define INTERSECT_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

/* TRUE if 0.0 < segment0,segment1 < 1.0 unless colinear */
GridBool intersectSegmentSegment(double *segment0_node0,double *segment0_node1,
				 double *segment1_node0,double *segment1_node1,
				 double *segment0, double *segment1 );

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
