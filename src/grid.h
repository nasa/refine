
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
long gridNNodes(Grid *g);
long gridNodeDeg(Grid *g, long nodeIndex);
Grid *gridRegisterCell(Grid *g, long n0, long n1, long n2, long n3);
void gridFree(Grid *g);

END_C_DECLORATION

#endif /* GRID_H */
