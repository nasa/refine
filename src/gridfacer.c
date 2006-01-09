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
#include <limits.h>
#include <values.h>

#include "queue.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "gridfacer.h"

GridFacer *gridfacerCreate( Grid *grid, int faceId )
{
  GridFacer *ge;
  ge = malloc(sizeof(GridFacer));
  ge->grid = grid;

  ge->faceId = faceId;

  gridAttachPacker( grid, gridfacerPack, (void *)ge );
  gridAttachNodeSorter( grid, gridfacerSortNode, (void *)ge );
  gridAttachReallocator( grid, gridfacerReallocator, (void *)ge );
  gridAttachFreeNotifier( grid, gridfacerGridHasBeenFreed, (void *)ge );

  return ge;
}

Grid *gridfacerGrid(GridFacer *ge)
{
  return ge->grid;
}

void gridfacerFree(GridFacer *ge)
{
  if (NULL != ge->grid) { 
    gridDetachPacker( ge->grid );
    gridDetachNodeSorter( ge->grid );
    gridDetachReallocator( ge->grid );
    gridDetachFreeNotifier( ge->grid );
  }
  free(ge);
}

void gridfacerPack(void *voidGridFacer, 
		   int nnode, int maxnode, int *nodeo2n,
		   int ncell, int maxcell, int *cello2n,
		   int nface, int maxface, int *faceo2n,
		   int nedge, int maxedge, int *edgeo2n)
{
  int i;
  GridFacer *ge = (GridFacer *)voidGridFacer;
}

void gridfacerSortNode(void *voidGridFacer, int maxnode, int *o2n)
{
  int i;
  GridFacer *ge = (GridFacer *)voidGridFacer;
}

void gridfacerReallocator(void *voidGridFacer, int reallocType, 
			  int lastSize, int newSize)
{

}

void gridfacerGridHasBeenFreed(void *voidGridFacer )
{
  GridFacer *ge = (GridFacer *)voidGridFacer;
  ge->grid = NULL;
}

int gridfacerFaceId(GridFacer *ge)
{
  return ge->faceId;
}

