
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
  int count;
  int *firstcell
};

Grid* grid_create(int n)
{
  int i;
  Grid *grid;
  grid = malloc(sizeof(Grid));
  grid->count =n;
  grid->firstcell = malloc(grid->count * sizeof(int));
  for (i=1;i < grid->count; i++ ) grid->firstcell[i] = NULL;
  return  grid;
}

int grid_count(Grid *grid)
{
  return grid->count;
}

int *grid_firstcell(Grid *grid, int id)
{
  return grid->firstcell[id];
}

void grid_free(Grid *grid)
{
  free(grid->firstcell);
  free(grid);
}
