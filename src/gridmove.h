
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
  GridBool *specified;
};

GridMove *gridmoveCreate(Grid *);
Grid *gridmoveGrid(GridMove *);
void gridmoveFree(GridMove *);
void gridmovePack(void *voidGridMove, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n);
void gridmoveSortNode(void *voidGridMove, int maxnode, int *o2n);
void gridmoveReallocator(void *voidGridMove, int reallocType, 
			 int lastSize, int newSize);
void gridmoveGridHasBeenFreed(void *voidGridMove );

GridMove *gridmoveDisplace(GridMove *, int node, double *displace);
GridMove *gridmoveDisplacement(GridMove *, int node, double *displacement);
GridBool gridmoveSpecified(GridMove *, int node);
GridMove *gridmoveMove(GridMove *);
GridMove *gridmoveSprings(GridMove *, int *nsprings, int **springs);

END_C_DECLORATION

#endif /* GRIDMOVE_H */
