
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


Grid *gridCreate(long nnode, long ncell);
long gridNNode(Grid *g);
long gridNCell(Grid *g);
long gridNodeDeg(Grid *g, long nodeIndex);
Grid *gridRegisterNodeCell(Grid *g, long nodeIndex, long cellIndex);
void gridFirstNodeCell(Grid *g, long nodeIndex);
void gridNextNodeCell(Grid *g);
long gridCurrectNodeCell(Grid *g);
int gridMoreNodeCell(Grid *g);

void gridFree(Grid *g);

END_C_DECLORATION

#endif /* GRID_H */
