
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

Grid *gridSwapEdge(Grid *g, int n0, int n1 );
Grid *gridSwap(Grid *g );

#endif /* GRIDSWAP_H */
