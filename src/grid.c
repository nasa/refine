
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

typedef struct N2C N2C;

struct N2C {
  int id;
  N2C *next;
};

struct Grid {
  int nnode;
  int maxcell;
  int ncell;
  int nlist;
  N2C **first;
  N2C *current;
  N2C *blank;
  N2C *n2c;
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

  if (grid->nlist < 1) grid->nlist = grid->maxcell*4;

  grid->first = (N2C **)malloc(grid->nnode * sizeof(N2C *));

  for (i=0;i < grid->nnode; i++ ) grid->first[i] = NULL; 
  grid->n2c = (N2C *)malloc( grid->nlist * sizeof(N2C));
  for (i=0;i < grid->nlist-1; i++ ) { // pointer majic?
    grid->n2c[i].id   = EMPTY;
    grid->n2c[i].next = &grid->n2c[i+1];
  }
  grid->n2c[grid->nlist-1].id   = EMPTY;
  grid->n2c[grid->nlist-1].next = NULL;
  grid->blank   = grid->n2c;
  grid->current = NULL;
  return  grid;
}

void gridFree(Grid *grid)
{
  free(grid->c2n);
  free(grid->first);
  free(grid->n2c);
  free(grid);
}

Grid *gridDump(Grid *grid)
{
  printf("\n Dump Not impl.\n");
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

bool gridCellExists(Grid *grid, int nodeId, int cellId)
{
  bool exist;
  exist = FALSE;
  for ( gridFirstNodeCell(grid,nodeId); 
	!exist && gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) 
    exist = (cellId == gridCurrentNodeCell(grid));
  return exist;
}

Grid* gridRegisterNodeCell(Grid *grid, int nodeId, int cellId)
{
  N2C *new;
  if (nodeId > grid->nnode-1) return NULL;
  if (grid->blank == NULL) return NULL;
  new = grid->blank;
  grid->blank = grid->blank->next;
  new->next = grid->first[nodeId];
  grid->first[nodeId] = new;
  new->id = cellId;

  return grid;
}

Grid* gridRemoveNodeCell(Grid *grid, int nodeId, int cellId)
{
  N2C *remove, *previous;
  remove = NULL;

  for ( gridFirstNodeCell(grid,nodeId); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) {
    if (gridCurrentNodeCell(grid)==cellId) 
      remove = grid->current;
  }

  if (remove == NULL) return NULL;
 
  previous = NULL;

  for ( gridFirstNodeCell(grid,nodeId); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) {
    if (grid->current != NULL && grid->current->next == remove) 
      previous = grid->current;
  }

  if ( previous == NULL ) {
    grid->first[nodeId] = remove->next;
  }else{
    previous->next = remove->next;
  }

  remove->id = EMPTY;
  remove->next = grid->blank;
  grid->blank = remove;

  return grid;
}

void gridFirstNodeCell(Grid *grid, int nodeId)
{
  grid->current = grid->first[nodeId];
}

void gridNextNodeCell(Grid *grid)
{
  if ( grid->current != NULL) grid->current = grid->current->next;
}

int gridCurrentNodeCell(Grid *grid)
{
  if (grid->current == NULL ) return EMPTY;
  return grid->current->id;
}

bool gridValidNodeCell(Grid *grid)
{
  return (grid->current != NULL);
}

bool gridMoreNodeCell(Grid *grid)
{
  return ( (grid->current != NULL) && (grid->current->next != NULL) );
}

Grid *gridAddCell(Grid *grid, int n0, int n1, int n2, int n3)
{
  int cellId,icell;
  if (grid->ncell >= grid->maxcell) return NULL;
  cellId = grid->ncell;
  grid->ncell++;
  
  grid->c2n[0+4*cellId] = n0;
  grid->c2n[1+4*cellId] = n1;
  grid->c2n[2+4*cellId] = n2;
  grid->c2n[3+4*cellId] = n3;
  
  if ( NULL == gridRegisterNodeCell( grid, n0, cellId ) ) return NULL;
  if ( NULL == gridRegisterNodeCell( grid, n1, cellId ) ) return NULL;
  if ( NULL == gridRegisterNodeCell( grid, n2, cellId ) ) return NULL;
  if ( NULL == gridRegisterNodeCell( grid, n3, cellId ) ) return NULL;
  
  return grid;
}

Grid *gridPack(Grid *grid)
{
  printf("\n Pack Not impl.\n");
  fflush(stdout);
  return grid;
}

Grid *gridGem(Grid *grid, int n0, int n1, int maxgem, int *ngem, int *gem )
{
  int cellId, inode, i;
  int c0,c1,c2,c3;
  *ngem = 0;
  
#define SAFE_NGEM_INC (*ngem)++;if(*ngem>(maxgem-1)){*ngem=0;return(NULL);}

  for ( gridFirstNodeCell(grid,n0); 
	gridValidNodeCell(grid); 
	gridNextNodeCell(grid)) {
    cellId = gridCurrentNodeCell(grid);
    c0 = grid->c2n[0+4*cellId];
    c1 = grid->c2n[1+4*cellId];
    c2 = grid->c2n[2+4*cellId];
    c3 = grid->c2n[3+4*cellId];

    /* 0 leads */
    if ( n0 == c0 && n1 == c1 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c0;
      gem[1+4*(*ngem-1)] = c1;
      gem[2+4*(*ngem-1)] = c2;
      gem[3+4*(*ngem-1)] = c3;
    }

    if ( n0 == c0 && n1 == c2 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c0;
      gem[1+4*(*ngem-1)] = c2;
      gem[2+4*(*ngem-1)] = c3;
      gem[3+4*(*ngem-1)] = c1;
    }

    if ( n0 == c0 && n1 == c3 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c0;
      gem[1+4*(*ngem-1)] = c3;
      gem[2+4*(*ngem-1)] = c1;
      gem[3+4*(*ngem-1)] = c2;
    }

    /* 1 leads */
    if ( n0 == c1 && n1 == c0 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c1;
      gem[1+4*(*ngem-1)] = c0;
      gem[2+4*(*ngem-1)] = c3;
      gem[3+4*(*ngem-1)] = c2;
    }

    if ( n0 == c2 && n1 == c0 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c2;
      gem[1+4*(*ngem-1)] = c0;
      gem[2+4*(*ngem-1)] = c1;
      gem[3+4*(*ngem-1)] = c3;
    }

    if ( n0 == c3 && n1 == c0 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c3;
      gem[1+4*(*ngem-1)] = c0;
      gem[2+4*(*ngem-1)] = c2;
      gem[3+4*(*ngem-1)] = c1;
    }

    /* 2 leads */
    if ( n0 == c2 && n1 == c3 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c2;
      gem[1+4*(*ngem-1)] = c3;
      gem[2+4*(*ngem-1)] = c0;
      gem[3+4*(*ngem-1)] = c1;
    }

    if ( n0 == c3 && n1 == c1 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c3;
      gem[1+4*(*ngem-1)] = c1;
      gem[2+4*(*ngem-1)] = c0;
      gem[3+4*(*ngem-1)] = c2;
    }

    if ( n0 == c1 && n1 == c2 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c1;
      gem[1+4*(*ngem-1)] = c2;
      gem[2+4*(*ngem-1)] = c0;
      gem[3+4*(*ngem-1)] = c3;
    }

    /* 3 leads */
    if ( n0 == c3 && n1 == c2 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c3;
      gem[1+4*(*ngem-1)] = c2;
      gem[2+4*(*ngem-1)] = c1;
      gem[3+4*(*ngem-1)] = c0;
    }

    if ( n0 == c1 && n1 == c3 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c1;
      gem[1+4*(*ngem-1)] = c3;
      gem[2+4*(*ngem-1)] = c2;
      gem[3+4*(*ngem-1)] = c0;
    }

    if ( n0 == c2 && n1 == c1 ) {
      SAFE_NGEM_INC;
      gem[0+4*(*ngem-1)] = c2;
      gem[1+4*(*ngem-1)] = c1;
      gem[2+4*(*ngem-1)] = c3;
      gem[3+4*(*ngem-1)] = c0;
    }


  }

  return grid;
}

Grid *gridEquator(Grid *grid, int n0, int n1, int maxequ, int *nequ, int *equ )
{
#define MAXGEM 200
  int ngem;
  int gem[MAXGEM];

  *nequ = 0;

  if ( NULL == gridGetGem( grid, n0, n1, MAXGEM, &ngem, gem ) ) return NULL;

  return grid;
}

