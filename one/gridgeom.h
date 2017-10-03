
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

#ifndef GRIDGEOM_H
#define GRIDGEOM_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridGeomStartOnly( Grid *, char *url, char *modler, char *project );

Grid *gridParallelGeomLoad( Grid *, char *url, char *modler, char *project );
Grid *gridParallelGeomSave( Grid *, char *project );

Grid *gridUpdateEdgeGrid( Grid *, int edgeId, int nCurveNode,
			  double *xyz, double *t );

int gridFaceEdgeCount( int faceId );
Grid *gridFaceEdgeLocal2Global( Grid *, int faceId, 
				int faceEdgeCount, int *local2global );

Grid *gridUpdateGeometryFace( Grid *, int faceId, 
			      int nnode, double *xyz, double *uv, 
			      int nface, int *f2n );

Grid *gridCreateShellFromFaces( Grid * );

END_C_DECLORATION

#endif /* GRIDGEOM_H */
