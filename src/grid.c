
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
  long nnodes;
  long *firstcell
};

Grid* grid_create(long n)
{
  long i;
  Grid *grid;
  grid = malloc(sizeof(Grid));
  grid->nnodes =n;
  grid->firstcell = malloc(grid->nnodes * sizeof(long));
  for (i=0;i < grid->nnodes; i++ ) grid->firstcell[i] = 0;
  return  grid;
}

long grid_nnodes(Grid *grid)
{
  return grid->nnodes;
}

long *grid_firstcell(Grid *grid, long id)
{
  return grid->firstcell[id];
}

void grid_free(Grid *grid)
{
  free(grid->firstcell);
  free(grid);
}
