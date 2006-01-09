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
  GridFacer *gf;
  gf = malloc(sizeof(GridFacer));
  gf->grid = grid;

  gf->faceId = faceId;

  gridAttachPacker( grid, gridfacerPack, (void *)gf );
  gridAttachNodeSorter( grid, gridfacerSortNode, (void *)gf );
  gridAttachReallocator( grid, gridfacerReallocator, (void *)gf );
  gridAttachFreeNotifier( grid, gridfacerGridHasBeenFreed, (void *)gf );

  return gf;
}

Grid *gridfacerGrid(GridFacer *gf)
{
  return gf->grid;
}

void gridfacerFree(GridFacer *gf)
{
  if (NULL != gf->grid) { 
    gridDetachPacker( gf->grid );
    gridDetachNodeSorter( gf->grid );
    gridDetachReallocator( gf->grid );
    gridDetachFreeNotifier( gf->grid );
  }
  free(gf);
}

void gridfacerPack(void *voidGridFacer, 
		   int nnode, int maxnode, int *nodeo2n,
		   int ncell, int maxcell, int *cello2n,
		   int nface, int maxface, int *faceo2n,
		   int nedge, int maxedge, int *edgeo2n)
{
  int i;
  GridFacer *gf = (GridFacer *)voidGridFacer;
}

void gridfacerSortNode(void *voidGridFacer, int maxnode, int *o2n)
{
  int i;
  GridFacer *gf = (GridFacer *)voidGridFacer;
}

void gridfacerReallocator(void *voidGridFacer, int reallocType, 
			  int lastSize, int newSize)
{

}

void gridfacerGridHasBeenFreed(void *voidGridFacer )
{
  GridFacer *gf = (GridFacer *)voidGridFacer;
  gf->grid = NULL;
}

int gridfacerFaceId(GridFacer *gf)
{
  return gf->faceId;
}

GridFacer *gridfacerExamine(GridFacer *gf)
{
  int face;
  int nodes[3];
  int faceId;

  double map[6];
  double d[3], e[3], v0[3], v1[3], v2[3];
  double h0, h1, h2;

  Grid *grid = gridfacerGrid( gf );

  for ( face = 0 ; face < gridMaxFace(grid) ; face++ ) {
    if ( grid == gridFace( grid, face, nodes, &faceId ) ) {
      if ( gridfacerFaceId(gf) == faceId ) {
	
	gridMap(grid, nodes[0], map);
	gridTriDiag3x3(map, d, e, v0, v1, v2);
	if ( !gridEigTriDiag3x3(d, e, v0, v1, v2 )) {
	  printf("%s: %d: gridfacerExamine: gridEigTriDiag3x3 FAILED.\n",
		 __FILE__,__LINE__);
	  return NULL;
	}
	/* the new EigTriDiag should be ortho-normal */
	gridEigOrtho3x3( v0, v1, v2 );

	h0 = 1.0/sqrt(d[0]);
	h1 = 1.0/sqrt(d[1]);
	h2 = 1.0/sqrt(d[2]);

	printf("face %d %f %f %f\n",face,h0,h1,h2);
      }
    }
  }
  return gf; 
}
