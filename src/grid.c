
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "grid.h"

struct Grid {
  long nnode;
  long ncell;
  long *firstcell;
  long currentcell;
  long *celllist;
};

//#define EBUG

Grid* gridCreate(long nnode, long ncell, long nlist)
{
  long i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->nnode = nnode;
  grid->ncell = ncell;
  /* pad one for [0] and one to terminate */
  if (nlist < 1) nlist = (grid->ncell*4+grid->nnode)+2;

  grid->firstcell = malloc(grid->nnode * sizeof(long));
  for (i=0;i < grid->nnode; i++ ) grid->firstcell[i] = nlist-1;

  grid->celllist = malloc( nlist * sizeof(long));
  for (i=0;i < nlist; i++ ) grid->celllist[i] = -(i+1);
  grid->celllist[nlist-1] =0;
  grid->currentcell=nlist-1;
 
#ifdef EBUG
  printf("\n Cre n %d c %d l %d\n", 
  nnode,ncell,nlist);
  fflush(stdout);
#endif

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

long gridNodeDeg(Grid *grid, long id)
{
  int n;
  n =0;
  for ( gridFirstNodeCell(grid,id); 
	gridMoreNodeCell(grid); 
	gridNextNodeCell(grid)) n++;
  return n;
}

Grid* gridRegisterNodeCell(Grid *grid, long nodeId, long cellId)
{
  long entry, terminator, nextOpen;

  entry = -grid->celllist[0];
  if (entry == 0 ) return NULL;
  if ( grid->celllist[entry] == 0 ) return NULL;
  terminator = entry+1;
  //  terminator = -grid->celllist[entry]; make sure this is same as ^
  if ( grid->celllist[terminator] == 0 ) return NULL;
  nextOpen = -grid->celllist[terminator];
 
#ifdef EBUG
  printf("\n Reg n %d c %3d e %d t %d n %d", 
  nodeId, cellId, entry, terminator, nextOpen);
  fflush(stdout);
#endif

  grid->firstcell[nodeId]=entry;
  grid->celllist[entry]=cellId+1;
  grid->celllist[terminator]=0;
  grid->celllist[0]=-nextOpen;

  return grid;
}

Grid* gridRemoveNodeCell(Grid *grid, long nodeId, long cellId)
{
  long cellIndex;
  cellIndex = EMPTY;

  for ( gridFirstNodeCell(grid,nodeId); 
	gridMoreNodeCell(grid); 
	gridNextNodeCell(grid)) {
    if (gridCurrentNodeCell(grid)==cellId) 
      cellIndex = grid->currentcell;
  }

  if (cellIndex == EMPTY) return NULL;

  grid->celllist[cellIndex] = -(cellIndex+1);

  return grid;
}

void gridFirstNodeCell(Grid *grid, long nodeId)
{
  grid->currentcell = grid->firstcell[nodeId];
  while (grid->celllist[grid->currentcell] < 0){
    grid->currentcell = -grid->celllist[grid->currentcell];
  }
}
void gridNextNodeCell(Grid *grid)
{
  grid->currentcell++;
  while (grid->celllist[grid->currentcell] < 0){
    grid->currentcell = -grid->celllist[grid->currentcell];
  }
}
long gridCurrentNodeCell(Grid *grid)
{
  return grid->celllist[grid->currentcell]-1;
}
int gridMoreNodeCell(Grid *grid)
{
  return (gridCurrentNodeCell(grid) != EMPTY);
}

void gridFree(Grid *grid)
{
  free(grid->firstcell);
  free(grid);
}
