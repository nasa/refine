
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
  int nnode;
  int maxcell;
  int ncell;
  int nlist;
  int *firstcell;
  int currentcell;
  int firstblankcell;
  int *celllist;
  int *c2n;
};

//#define EBUG

Grid* gridCreate(int nnode, int maxcell, int nlist)
{
  int i;
  Grid *grid;

  grid = malloc(sizeof(Grid));

  grid->nnode = nnode;
  grid->maxcell = maxcell;
  grid->ncell = 0;
  grid->c2n = malloc(4 * grid->maxcell * sizeof(int));

  grid->nlist = nlist;

  /* pad one for [0] and one to terminate */
  if (grid->nlist < 1) grid->nlist = (grid->maxcell*4+grid->nnode)+2;

  grid->firstcell = malloc(grid->nnode * sizeof(int));
  /* I could set first cell to zero to terminate when realloc used */
  for (i=0;i < grid->nnode; i++ ) grid->firstcell[i] = 0; 

  grid->celllist = malloc( grid->nlist * sizeof(int));
  for (i=0;i < grid->nlist; i++ ) grid->celllist[i] = -(i+1);
  grid->celllist[grid->nlist-1] =0;
  grid->celllist[0] =0;
  grid->firstblankcell = 1;
  grid->currentcell=grid->nlist-1;
 
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
  printf("\nfirst blank %d\n",grid->firstblankcell);
  for ( i=0; i<grid->nnode; i++) printf("fc%5d -> %5d\n",i,grid->firstcell[i]);
  for ( i=0; i<grid->nlist; i++) printf("cl%5d -> %5d\n",i,grid->celllist[i]);
  fflush(stdout);
  return grid;
}

int gridNNode(Grid *grid)
{
  return grid->nnode;
}

int gridMaxCell(Grid *grid)
{
  return grid->maxcell;
}

int gridNCell(Grid *grid)
{
  return grid->ncell;
}

int gridNodeDeg(Grid *grid, int id)
{
  int n;
  n =0;
  for ( gridFirstNodeCell(grid,id); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) n++;
  return n;
}

int gridCellExists(Grid *grid, int nodeId, int cellId)
{
  int exist;
  exist = (0==1);
  for ( gridFirstNodeCell(grid,nodeId); 
	!exist && gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) 
    exist = (cellId == gridCurrentNodeCell(grid));
  return exist;
}

Grid* gridRegisterNodeCell(Grid *grid, int nodeId, int cellId)
{
  int firstAvailable, oldTerminator, entry, terminator, nextOpen;

  if (nodeId >= grid->nnode) return NULL;

  if ( grid->firstcell[nodeId] == 0 ) {
    entry = grid->firstblankcell;
    if ( entry == 0 ) return NULL;
    if ( grid->celllist[entry] == 0 ) return NULL;
    terminator = -grid->celllist[entry];
    if ( terminator == 0 ) return NULL;
    nextOpen = -grid->celllist[terminator];
    grid->firstcell[nodeId] = entry;
    grid->celllist[entry] = cellId+1;
    grid->celllist[terminator] = 0;
    grid->firstblankcell = nextOpen;
  }else{

    for ( gridFirstNodeCell(grid,nodeId);
	  gridValidNodeCell(grid);
	  gridNextNodeCell(grid));  // replace a removed here

    oldTerminator = grid->currentcell;

    entry = grid->firstblankcell;
    if ( entry == 0 ) return NULL;
    
    if (entry == (oldTerminator+1) ) { // append to current list
      grid->celllist[oldTerminator] = cellId+1;
      grid->firstblankcell = -grid->celllist[entry];
      grid->celllist[entry] = 0;      
    }else{
      if ( grid->celllist[entry] == 0 ) return NULL;
      terminator = entry+1;
      grid->celllist[oldTerminator] = -entry;      
      grid->celllist[entry] = cellId+1;
      grid->firstblankcell = -grid->celllist[terminator];      
      grid->celllist[terminator] = 0;
    }

  }

  return grid;
}

Grid* gridRemoveNodeCell(Grid *grid, int nodeId, int cellId)
{
  int cellIndex;
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

void gridFirstNodeCell(Grid *grid, int nodeId)
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

int gridCurrentNodeCell(Grid *grid)
{
  return grid->celllist[grid->currentcell]-1;
}

int gridValidNodeCell(Grid *grid)
{
  return (gridCurrentNodeCell(grid) != EMPTY);
}

int gridMoreNodeCell(Grid *grid)
{
  int next;

  if ( !gridValidNodeCell(grid) ) return 0;
  next = grid->currentcell + 1;
  while (grid->celllist[next] < 0){
    next = -grid->celllist[next];
  }
  return (grid->celllist[next] != 0);
}

Grid *gridAddCell(Grid *grid, int n0, int n1, int n2, int n3)
{
  int cellId,icell;
  cellId = grid->ncell;
  grid->ncell++;
  
  grid->c2n[0+4*cellId] = n0;
  grid->c2n[1+4*cellId] = n1;
  grid->c2n[2+4*cellId] = n2;
  grid->c2n[3+4*cellId] = n3;
  
  gridRegisterNodeCell( grid, n0, cellId );
  gridRegisterNodeCell( grid, n1, cellId );
  gridRegisterNodeCell( grid, n2, cellId );
  gridRegisterNodeCell( grid, n3, cellId );
  
  return grid;
}

Grid *gridPack(Grid *grid)
{
  int nodeId, current, deg, i;
  current = 1;
  for ( nodeId=0 ; nodeId<grid->nnode ; nodeId++) {
    deg = gridNodeDeg( grid, nodeId );
    if ( deg == 0 ) {
      grid->firstcell[nodeId] = 0;
    }else{
      grid->firstcell[nodeId] = current;
      for ( gridFirstNodeCell(grid,nodeId); 
	    gridValidNodeCell(grid); 
	    gridNextNodeCell(grid)){
	grid->celllist[current] = gridCurrentNodeCell(grid)+1;
	current++;
      }
      grid->celllist[current] = 0;
      current++;
    }
  }
  grid->firstblankcell = current;
  for (i=current;i < grid->nlist; i++ ) grid->celllist[i] = -(i+1);
  grid->celllist[grid->nlist-1] =0;
  return grid;
}

Grid *gridGetGem(Grid *grid, int n0, int n1, int maxgem, int *ngem, int *gem )
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

