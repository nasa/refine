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
#include "gridmove.h"

GridMove *gridmoveCreate( Grid *grid )
{
  int i;
  GridMove *gm;
  gm = malloc(sizeof(GridMove));
  gm->grid = grid;

  gridAttachNodeSorter( grid, gridmoveNodeSorter, gm );
  gridAttachReallocator( grid, gridmoveReallocator, gm );

  gm->displacement = malloc(3*gridMaxNode(grid)*sizeof(double));
  for (i=0;i<3*gridMaxNode(grid);i++) gm->displacement[i] = 0.0;

  return gm;
}

Grid *gridmoveGrid(GridMove *gm)
{
  return gm->grid;
}

void gridmoveFree(GridMove *gm)
{
  free(gm->displacement);
  gridDetachNodeSorter( gm->grid );
  gridDetachReallocator( gm->grid );
  free(gm);
}

void gridmoveNodeSorter(void *voidGridMove, int maxnode, int *o2n)
{
  GridMove *gm = (GridMove *)voidGridMove;
  int node, ixyz;
  double *temp_xyz;
  temp_xyz = malloc( maxnode * sizeof(double) );
  for ( ixyz = 0; ixyz < 3 ; ixyz++ ){
    for ( node = 0 ; node < maxnode ; node++ ){
      temp_xyz[o2n[node]] = gm->displacement[ixyz+3*node];
    }
    for ( node = 0 ; node < maxnode ; node++ ){
      gm->displacement[ixyz+3*node] = temp_xyz[node];
    }
  }
  free(temp_xyz);
}

void gridmoveReallocator(void *voidGridMove, int reallocType, 
			 int lastSize, int newSize)
{
  GridMove *gm = (GridMove *)voidGridMove;
  int i;
  if (gridREALLOC_NODE == reallocType) {
    gm->displacement = realloc(gm->displacement, 3*newSize*sizeof(double));
    for (i=3*lastSize;i<3*newSize;i++) gm->displacement[i] = 0.0;
  }
}

GridMove *gridmoveDisplace(GridMove *gm, int node, double *displace)
{
  if (node < 0 || node >= gridMaxNode(gm->grid)) return NULL;
  gm->displacement[0+3*node] = displace[0];
  gm->displacement[1+3*node] = displace[1];
  gm->displacement[2+3*node] = displace[2];
  return gm;
}

GridMove *gridmoveDisplacement(GridMove *gm, int node, double *displacement)
{
  if (node < 0 || node >= gridMaxNode(gm->grid)) return NULL;
  displacement[0] = gm->displacement[0+3*node];
  displacement[1] = gm->displacement[1+3*node];
  displacement[2] = gm->displacement[2+3*node];
  return gm;
}
