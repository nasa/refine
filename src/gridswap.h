
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDSWAP_H
#define GRIDSWAP_H

#include "master_header.h"
#include "grid.h"

Grid *gridSwapFace(Grid *g, int n0, int n1, int n2 );
Grid *gridSwapEdge(Grid *g, int n0, int n1 );
Grid *gridSwapNearNode(Grid *g, int node );
Grid *gridSwapNearNodeExceptBoundary(Grid *g, int node );
Grid *gridSwap(Grid *g );
Grid *gridSwapCellFaceArea(Grid *g, int cell );
Grid *gridSwapFaceArea(Grid *g );

#endif /* GRIDSWAP_H */
