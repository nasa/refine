
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDINSERT_H
#define GRIDINSERT_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION
Grid *gridThrash(Grid *g );
Grid *gridSplitEdge(Grid *g, int n0, int n1 );

END_C_DECLORATION

#endif /* GRIDINSERT_H */
