
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDCAD_H
#define GRIDCAD_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridProjectNodeToEdge(Grid *g, int node, int edgeId );
Grid *gridSafeProject(Grid *g, int node );

END_C_DECLORATION

#endif /* GRIDCAD_H */
