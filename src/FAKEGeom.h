
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef CADGEOM_H
#define CADGEOM_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

bool CADGeom_NearestOnEdge(int vol, int edgeId, 
			   double *xyz, double *t, double *xyznew);
bool CADGeom_NearestOnFace(int vol, int faceId, 
			   double *xyz, double *uv, double *xyznew);

bool CADGeom_PointOnEdge(int vol, int edgeId, 
			 double t, double *xyz, 
			 int derivativeFlag, double *dt, double *dtdt );
bool CADGeom_PointOnFace(int vol, int faceId, 
			 double *uv, double *xyz, 
			 int derivativeFlag, double *du, double *dv,
			 double *dudu, double *dudv, double *dvdv );

END_C_DECLORATION

#endif /* CADGEOM_H */
