
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
  long maxcell;
  long ncell;
  long nlist;
  long *firstcell;
  long currentcell;
  long *celllist;
  long *c2n;
};

//#define EBUG

Grid* gridCreate(long nnode, long maxcell, long nlist)
{
  long i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->nnode = nnode;
  grid->maxcell = maxcell;
  grid->ncell = 0;
  grid->c2n = malloc(4 * grid->maxcell * sizeof(long));

  grid->nlist = nlist;

  /* pad one for [0] and one to terminate */
  if (grid->nlist < 1) grid->nlist = (grid->maxcell*4+grid->nnode)+2;

  grid->firstcell = malloc(grid->nnode * sizeof(long));
  /* I could set first cell to zero to terminate when realloc used */
  for (i=0;i < grid->nnode; i++ ) grid->firstcell[i] = grid->nlist-1; 

  grid->celllist = malloc( grid->nlist * sizeof(long));
  for (i=0;i < grid->nlist; i++ ) grid->celllist[i] = -(i+1);
  grid->celllist[grid->nlist-1] =0;
  grid->currentcell=grid->nlist-1;
 
#ifdef EBUG
  printf("\n Cre n %d c %d l %d\n", 
  grid->nnode,grid->maxcell,grid->nlist);
  fflush(stdout);
#endif

  return  grid;
}

void gridFree(Grid *grid)
{
  free(grid->c2n);
  free(grid->firstcell);
  free(grid);
}

Grid *gridDump(Grid *grid)
{
  int i;
  printf("\n");
  for ( i=0; i<grid->nnode; i++) printf("fc%5d -> %5d\n",i,grid->firstcell[i]);
  for ( i=0; i<grid->nlist; i++) printf("cl%5d -> %5d\n",i,grid->celllist[i]);
  fflush(stdout);
  return grid;
}

long gridNNode(Grid *grid)
{
  return grid->nnode;
}

long gridMaxCell(Grid *grid)
{
  return grid->maxcell;
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
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) n++;
  return n;
}

int gridCellExists(Grid *grid, long nodeId, long cellId)
{
  int exist;
  exist = (0==1);
  for ( gridFirstNodeCell(grid,nodeId); 
	!exist && gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) 
    exist = (cellId == gridCurrentNodeCell(grid));
  return exist;
}

Grid* gridRegisterNodeCell(Grid *grid, long nodeId, long cellId)
{
  long firstAvailable, oldTerminator, entry, terminator, nextOpen;

#ifdef EBUG
  gridDump(grid);
#endif

  firstAvailable = -grid->celllist[0];
  if ( firstAvailable == 0 ) return NULL; /* list full */

  for ( gridFirstNodeCell(grid,nodeId); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid));

  oldTerminator = grid->currentcell;

  if (oldTerminator == (firstAvailable-1) ) {
    entry = oldTerminator;
    grid->celllist[entry] = -firstAvailable;
  }else{
    entry = firstAvailable;
  }

  if ( grid->celllist[entry] == 0 ) return NULL;
  terminator = entry+1;
  nextOpen = -grid->celllist[terminator];
  if ( grid->celllist[terminator] == 0 ) nextOpen = terminator;
 
#ifdef EBUG
  printf("\n FReg n %d c %3d e %d t %d n %d", 
  nodeId, cellId, entry, terminator, nextOpen);
  fflush(stdout);
#endif

  grid->celllist[terminator]=0;
  if (grid->celllist[grid->firstcell[nodeId]]==0)
    grid->firstcell[nodeId]=entry;
  grid->celllist[entry]=cellId+1;
  grid->celllist[0]=-nextOpen;

  return grid;
}

Grid* gridRemoveNodeCell(Grid *grid, long nodeId, long cellId)
{
  long cellIndex;
  cellIndex = EMPTY;

  for ( gridFirstNodeCell(grid,nodeId); 
	gridValidNodeCell(grid); 
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
  if ( !gridValidNodeCell(grid) ) return;
  grid->currentcell++;
  while (grid->celllist[grid->currentcell] < 0){
    grid->currentcell = -grid->celllist[grid->currentcell];
  }
}

long gridCurrentNodeCell(Grid *grid)
{
  return grid->celllist[grid->currentcell]-1;
}

int gridValidNodeCell(Grid *grid)
{
  return (gridCurrentNodeCell(grid) != EMPTY);
}

int gridMoreNodeCell(Grid *grid)
{
  long next;

  if ( !gridValidNodeCell(grid) ) return 0;
  next = grid->currentcell + 1;
  while (grid->celllist[next] < 0){
    next = -grid->celllist[next];
  }
  return (grid->celllist[next] != 0);
}

Grid *gridAddCell(Grid *grid, long n0, long n1, long n2, long n3)
{
  long cellId,icell;
  cellId = 0;
  icell = grid->ncell;
  grid->ncell++;
  
  grid->c2n[0+4*icell] = n0;
  grid->c2n[1+4*icell] = n1;
  grid->c2n[2+4*icell] = n2;
  grid->c2n[3+4*icell] = n3;
  
  gridRegisterNodeCell( grid, n0, cellId );
  gridRegisterNodeCell( grid, n1, cellId );
  gridRegisterNodeCell( grid, n2, cellId );
  gridRegisterNodeCell( grid, n3, cellId );
  
  return grid;
}

Grid *gridGetGem(Grid *grid, long n0, long n1, int maxgem, int *ngem, int *gem )
{
  int cellId, inode;
  *ngem = 0;
  
  for ( gridFirstNodeCell(grid,n0); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) {
    cellId = gridCurrentNodeCell(grid);
    for ( inode = 0 ; inode < 4 ; inode++) {
      if ( grid->c2n[inode+4*cellId] == n1 ) {
	(*ngem)++;
	if ( *ngem > (maxgem-1) ) return NULL;
	gem[*ngem-1] = cellId;
      }
    }
  }

  return grid;
}

