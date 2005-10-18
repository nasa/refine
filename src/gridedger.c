/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridedger.h"

GridEdger *gridedgerCreate( Grid *grid )
{
  int i;
  GridEdger *gm;
  gm = malloc(sizeof(GridEdger));
  gm->grid = grid;

  gridAttachPacker( grid, gridedgerPack, (void *)gm );
  gridAttachNodeSorter( grid, gridedgerSortNode, (void *)gm );
  gridAttachReallocator( grid, gridedgerReallocator, (void *)gm );
  gridAttachFreeNotifier( grid, gridedgerGridHasBeenFreed, (void *)gm );

  return gm;
}

Grid *gridedgerGrid(GridEdger *gm)
{
  return gm->grid;
}

void gridedgerFree(GridEdger *gm)
{
  if (NULL != gm->grid) { 
    gridDetachPacker( gm->grid );
    gridDetachNodeSorter( gm->grid );
    gridDetachReallocator( gm->grid );
    gridDetachFreeNotifier( gm->grid );
  }
  free(gm);
}

void gridedgerPack(void *voidGridEdger, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerSortNode(void *voidGridEdger, int maxnode, int *o2n)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerReallocator(void *voidGridEdger, int reallocType, 
			 int lastSize, int newSize)
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
}

void gridedgerGridHasBeenFreed(void *voidGridEdger )
{
  GridEdger *gm = (GridEdger *)voidGridEdger;
  gm->grid = NULL;
}

