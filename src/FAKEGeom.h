
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef CADGEOM_H
#define CADGEOM_H

#include "master_header.h"

BEGIN_C_DECLORATION

bool CADGeom_NearestOnEdge(int vol, int edge, 
			   double *xyz, double *t, double *xyznew);

END_C_DECLORATION

#endif /* CADGEOM_H */
