
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include "CADGeom/CADGeom.h"

bool CADGeom_NearestOnEdge(int vol, int edge, 
			   double *xyz, double *t, double *xyznew)
{
  *t = xyz[0];
  xyznew[0] = xyz[0];
  xyznew[1] = 0.0;
  xyznew[2] = 0.0;

  return TRUE;
}
