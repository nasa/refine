
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

GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal);

END_C_DECLORATION

#endif /* CADGEOM_H */
