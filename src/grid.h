
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRID_H
#define GRID_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Grid Grid;


Grid *gridCreate(long nnodes);
long gridNNodes(Grid *n);
void gridFree(Grid *n);

END_C_DECLORATION

#endif /* GRID_H */
