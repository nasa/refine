
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "grid.h"

struct Grid {
  long nnode;
  long ncell;
  long *firstcell;
  long *celllist;
};

Grid* gridCreate(long nnode, long ncell)
{
  long i, nlist;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->nnode = nnode;
  grid->firstcell = malloc(grid->nnode * sizeof(long));
  for (i=0;i < grid->nnode; i++ ) grid->firstcell[i] = 0;

  grid->ncell = ncell;
  nlist = (grid->ncell*4+grid->nnode)+1;
  grid->celllist = malloc( nlist * sizeof(long));
  for (i=0;i < nlist; i++ ) grid->celllist[i] = -(i+1);
  grid->celllist[nlist-1] =0;
  
 
  return  grid;
}

long gridNNode(Grid *grid)
{
  return grid->nnode;
}

long gridNCell(Grid *grid)
{
  return grid->ncell;
}

long* gridDEBUGcelllist(Grid *grid)
{
  return grid->celllist;
}

long gridNodeDeg(Grid *grid, long id)
{
  return grid->firstcell[id];
}

Grid* gridRegisterNodeCell(Grid *grid, long id)
{
  grid->firstcell[id]++;
  return grid;
}

void gridFirstNodeCell(Grid *grid, long id)
{
}
void gridNextNodeCell(Grid *grid)
{
}
int gridLastNodeCell(Grid *grid)
{
  return 0;
}

void gridFree(Grid *grid)
{
  free(grid->firstcell);
  free(grid);
}
