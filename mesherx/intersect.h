
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
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
