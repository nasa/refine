
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDGEOM_H
#define GRIDGEOM_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridParallelGeomLoad( Grid *, char *project );
Grid *gridParallelGeomSave( Grid *, char *project );
Grid *gridGeomSize( Grid *, int *nGeomNode, int *nGeomEdge, int *nGeomFace );

END_C_DECLORATION

#endif /* GRIDGEOM_H */
