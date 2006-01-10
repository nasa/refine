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
  int face;
  int nodes[3];
  int id;

  GridFacer *gf;
  gf = malloc(sizeof(GridFacer));
  gf->grid = grid;

  gf->faceId = faceId;

  gridAttachPacker( grid, gridfacerPack, (void *)gf );
  gridAttachNodeSorter( grid, gridfacerSortNode, (void *)gf );
  gridAttachReallocator( grid, gridfacerReallocator, (void *)gf );
  gridAttachFreeNotifier( grid, gridfacerGridHasBeenFreed, (void *)gf );

  gf->nedge = 0;
  gf->maxedge = 0;
  gf->e2n = NULL;

  for ( face = 0 ; face < gridMaxFace(grid) ; face++ ) {
    if ( grid == gridFace( grid, face, nodes, &faceId ) ) {
      if ( gridfacerFaceId(gf) == faceId ) {
	gridfacerAddUniqueEdge( gf, nodes[0], nodes[1] );
	gridfacerAddUniqueEdge( gf, nodes[1], nodes[2] );
	gridfacerAddUniqueEdge( gf, nodes[2], nodes[0] );
      }
    }
  }

  return gf;
}

Grid *gridfacerGrid(GridFacer *gf)
{
  return gf->grid;
}

void gridfacerFree(GridFacer *gf)
{
  if (NULL != gf->e2n) free(gf->e2n);
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
  for ( i = 0 ; i < 2*gridfacerEdges(gf) ; i++ ) {
    gf->e2n[i] = nodeo2n[gf->e2n[i]];
  }
}

void gridfacerSortNode(void *voidGridFacer, int maxnode, int *o2n)
{
  int i;
  GridFacer *gf = (GridFacer *)voidGridFacer;
  for ( i = 0 ; i < 2*gridfacerEdges(gf) ; i++ ) {
    gf->e2n[i] = o2n[gf->e2n[i]];
  }
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

int gridfacerEdges(GridFacer *gf)
{
  return gf->nedge;
}

GridFacer *gridfacerAddUniqueEdge(GridFacer *gf, int nodeA, int nodeB)
{
  int node0, node1;
  int edge;
  int chunk;

  node0 = MIN(nodeA,nodeB);
  node1 = MAX(nodeA,nodeB);

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    if ( (node0 == gf->e2n[0+2*edge]) && 
	 (node1 == gf->e2n[1+2*edge]) ) return gf; 
  }

  gf->nedge++;
  if ( gf->nedge >= gf->maxedge ) {
    chunk = 2000;
    gf->maxedge += chunk;
    if ( NULL == gf->e2n ) {
      gf->e2n = (int *)malloc(2 * gf->maxedge * sizeof(int));
    } else {
      gf->e2n = (int *)realloc( gf->e2n, 2 * gf->maxedge * sizeof(int));
    }
  }
  edge = gf->nedge - 1;

  gf->e2n[0+2*edge]=node0; 
  gf->e2n[1+2*edge]=node1;

  return gf;
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
