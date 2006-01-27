
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDSWAP_H
#define GRIDSWAP_H

#include "refine_defs.h"
#include "grid.h"
#include "queue.h"

Grid *gridSwapFace(Grid *g, Queue *q, int n0, int n1, int n2 );
Grid *gridSwapEdge(Grid *g, Queue *q, int n0, int n1 );
Grid *gridSwapNearNode(Grid *g, int node, double limit );
Grid *gridSwapNearNodeExceptBoundary(Grid *g, int node );
/* improvementLimit is set to default if less than 0.0 */
Grid *gridSwap(Grid *g, double improvementLimit );
Grid *gridRemoveTwoFaceCell(Grid *g, Queue *q, int cell );
Grid *gridRemoveThreeFaceCell(Grid *g, Queue *q, int cell );

#endif /* GRIDSWAP_H */
