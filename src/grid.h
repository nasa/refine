
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

Grid *gridCreate(int nnode, int ncell);
void gridFree(Grid *g);

int gridMaxNode(Grid *g);
int gridNNode(Grid *g);
int gridMaxCell(Grid *g);
int gridNCell(Grid *g);
int gridNodeDeg(Grid *g, int nodeIndex);
bool gridCellExists(Grid *g, int nodeIndex, int cellIndex);
Grid *gridRegisterNodeCell(Grid *g, int nodeIndex, int cellIndex);
Grid *gridRemoveNodeCell(Grid *g, int nodeIndex, int cellIndex);
void gridFirstNodeCell(Grid *g, int nodeIndex);
void gridNextNodeCell(Grid *g);
int gridCurrentNodeCell(Grid *g);
bool gridValidNodeCell(Grid *g);
bool gridMoreNodeCell(Grid *g);

Grid *gridAddCell(Grid *g, int n0, int n1, int n2, int n3 );
Grid *gridRemoveCell(Grid *g, int cellId );
Grid *gridCell(Grid *g, int cellId, int *nodes );

Grid *gridMakeGem(Grid *g, int n0, int n1 );
int gridNGem(Grid *g );
int gridGem(Grid *g, int index );

Grid *gridOrient(Grid *g, int *cell, int *nodes );
Grid *gridEquator(Grid *g, int n0, int n1 );
int gridNEqu(Grid *g );
int gridEqu(Grid *g, int index );
Grid *gridSwap(Grid *g, int n0, int n1 );

int gridAddNode(Grid *g, double x, double y, double z );

double gridVolume(Grid *g, int *nodes );
double gridAR(Grid *g, int *nodes );
double gridMinVolume(Grid *g);

END_C_DECLORATION

#endif /* GRID_H */
