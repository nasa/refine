
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


Grid *gridCreate(int nnode, int ncell, int nlist);
int gridNNode(Grid *g);
int gridNCell(Grid *g);
int gridNodeDeg(Grid *g, int nodeIndex);
int gridCellExists(Grid *g, int nodeIndex, int cellIndex);
Grid *gridRegisterNodeCell(Grid *g, int nodeIndex, int cellIndex);
Grid *gridRemoveNodeCell(Grid *g, int nodeIndex, int cellIndex);
void gridFirstNodeCell(Grid *g, int nodeIndex);
void gridNextNodeCell(Grid *g);
int gridCurrentNodeCell(Grid *g);
int gridValidNodeCell(Grid *g);
int gridMoreNodeCell(Grid *g);

Grid *gridAddCell(Grid *g, int n0, int n1, int n2, int n3 );
Grid *gridGetGem(Grid *g, int n0, int n1, int maxgem,int *ngem, int *gem );

Grid *gridDump(Grid *g);

void gridFree(Grid *g);

END_C_DECLORATION

#endif /* GRID_H */
