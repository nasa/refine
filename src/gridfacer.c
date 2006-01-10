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
  int edge;
  int nodes[2];
  double ratio;
  Grid *grid = gridfacerGrid( gf );

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    nodes[0] = gf->e2n[0+2*edge];
    nodes[1] = gf->e2n[1+2*edge];
    if ( 0 < gridParentGeometry(grid,nodes[0],nodes[1] ) ) {
      ratio = gridEdgeRatio(grid,nodes[0],nodes[1]);
      printf("edge%8d ratio%8.4f\n",edge,ratio);
    }
  }
  return gf; 
}
GridFacer *gridfacerSwap(GridFacer *gf)
{
  int edge;
  int node0, node1;
  int face0, face1;
  int nodes[3], faceId;
  int node2, node3;
  double current_ratio;
  double swapped_ratio;

  Grid *grid = gridfacerGrid( gf );

  edge = 0;
  while ( edge < gridfacerEdges(gf) ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 > gridParentGeometry(grid,node0,node1 ) ) {
      edge++;
      continue;
    }
    face0 = gridFindFaceWithNodesUnless(grid, node0, node1, EMPTY);
    if (EMPTY == face0) {
      printf("%s: %d: face0 EMPTY.\n",__FILE__,__LINE__);
      return NULL;
    }
    face1 = gridFindFaceWithNodesUnless(grid, node0, node1, face0);
    if (EMPTY == face1) {
      printf("%s: %d: face1 EMPTY for interior face.\n",__FILE__,__LINE__);
      return NULL;
    }
    gridFace(grid,face0,nodes,&faceId);
    node2 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;
    gridFace(grid,face1,nodes,&faceId);
    node3 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;

    current_ratio = gridEdgeRatio(grid,node0,node1);
    swapped_ratio = gridEdgeRatio(grid,node2,node3);
    printf("edge%8d current%8.4f swapped%8.4f\n",
	   edge,current_ratio,swapped_ratio);
    edge++;
  }
  return gf; 
}
