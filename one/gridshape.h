
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

#ifndef GRIDSHAPE_H
#define GRIDSHAPE_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridCurvedEdgeMidpoint(Grid *,int node0, int node1, double *midpoint);

double gridMinCellJacDet2(Grid *, int *nodes);
Grid *gridMinCellJacDetDeriv2(Grid *, int *nodes,
			      double *determinate, double *dDetdx);
Grid *gridNodeMinCellJacDet2(Grid *, int node, double *determinate );

Grid *gridPlotMinDeterminateAtSurface(Grid *);

Grid *gridWriteTecplotCurvedGeom(Grid *, char *filename );
Grid *gridWriteTecplotCellJacDet(Grid *, int cell, char *filename );

Grid *gridShapeJacobian1(Grid *,
			 double *n0, double *n1, double *n2, double *n3,
			 double *jacobian );
Grid *gridShapeJacobian2(Grid *,
			 double *n0, double *n1, double *n2, double *n3,
			 double *e01, double *e02, double *e03,
			 double *e12, double *e13, double *e23,
			 double *where, double *jacobian );
double gridShapeJacobianDet2(Grid *,
			     double *n0, double *n1, double *n2, double *n3,
			     double *e01, double *e02, double *e03,
			     double *e12, double *e13, double *e23,
			     double *where);

Grid *gridShapeJacobianDetDeriv2(Grid *,
				 double *n0, double *n1, double *n2, double *n3,
				 double *e01, double *e02, double *e03,
				 double *e12, double *e13, double *e23,
				 double *where,
				 double *determinate, double *dDetdx);

END_C_DECLORATION

#endif /* GRIDSHAPE_H */
