
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMOVE_H
#define GRIDMOVE_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

typedef struct GridMove GridMove;
struct GridMove {
  Grid *grid;
  void *gridRubyVALUEusedForGC;
  double *displacement;
};

GridMove *gridmoveCreate(Grid *);
Grid *gridmoveGrid(GridMove *);
void gridmoveFree(GridMove *);
void gridmoveNodeSorter(void *gridmove, int maxnode, int *o2n);
void gridmoveReallocator(void *gridmove, int reallocType, 
			 int lastSize, int newSize);

GridMove *gridmoveDisplace(GridMove *, int node, double *displace);
GridMove *gridmoveDisplacement(GridMove *, int node, double *displacement);

END_C_DECLORATION

#endif /* GRIDMOVE_H */
