
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include "FAKEGeom.h"

GridBool CADGeom_NearestOnEdge(int vol, int edgeId, 
			   double *xyz, double *t, double *xyznew)
{
  *t = xyz[0];
  xyznew[0] = xyz[0];
  xyznew[1] = 0.0;
  xyznew[2] = 0.0;

  return TRUE;
}

GridBool CADGeom_NearestOnFace(int vol, int faceId, 
			   double *xyz, double *uv, double *xyznew)
{
  uv[0] = xyz[0]+10.0;
  uv[1] = xyz[1]+20.0;
  xyznew[0] = xyz[0];
  xyznew[1] = xyz[1];
  xyznew[2] = 0.0;

  return TRUE;
}

GridBool CADGeom_PointOnEdge(int vol, int edgeId, 
			 double t, double *xyz, 
			 int derivativeFlag, double *dt, double *dtdt )
{
  xyz[0] = t;
  xyz[1] = 0.0;
  xyz[2] = 0.0;

  if (derivativeFlag > 0){
    dt[0] = 1.0;
    dt[1] = 0.0;
    dt[2] = 0.0;
  }

  if (derivativeFlag > 1){
    dtdt[0] = 0.0;
    dtdt[1] = 0.0;
    dtdt[2] = 0.0;
  }

  return TRUE;
}

GridBool CADGeom_PointOnFace(int vol, int faceId, 
			 double *uv, double *xyz, 
			 int derivativeFlag, double *du, double *dv,
			 double *dudu, double *dudv, double *dvdv )
{
  xyz[0] = uv[0]-10.0;
  xyz[1] = uv[1]-20.0;
  xyz[2] = 0.0;

  if (derivativeFlag > 0){
    du[0] = 1.0;
    du[1] = 0.0;
    du[2] = 0.0;
    dv[0] = 0.0;
    dv[1] = 1.0;
    dv[2] = 0.0;
  }

  if (derivativeFlag > 1){
    dudu[0] = 0.0;
    dudu[1] = 0.0;
    dudu[2] = 0.0;
    dudv[0] = 0.0;
    dudv[1] = 0.0;
    dudv[2] = 0.0;
    dvdv[0] = 0.0;
    dvdv[1] = 0.0;
    dvdv[2] = 0.0;
  }

  return TRUE;
}

GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal)
{
  xyz[0] = uv[0]-10.0;
  xyz[1] = uv[1]-20.0;
  xyz[2] = 0.0;

  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 1.0;

  return TRUE;
}
