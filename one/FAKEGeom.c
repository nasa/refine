
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

  normal[0] =  0.0;
  normal[1] =  0.0;
  normal[2] = -1.0;

  return TRUE;
}
