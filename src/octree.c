
/* Octree data structure for fast geometric searches
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "octree.h"

typedef struct Octant Octant;
struct Octant {
  double data[9];
};

struct Octree {
  double xMin, xMax, yMin, yMax, zMin, zMax;
  int nOctant;
  Octant *octant;
};


Octree* octreeCreate( double xMin, double xMax, 
		      double yMin, double yMax, 
		      double zMin, double zMax )
{
  Octree *octree;

  octree = malloc( sizeof(Octree) );
  octree->xMin = xMin;
  octree->xMax = xMax;
  octree->yMin = yMin;
  octree->yMax = yMax;
  octree->zMin = zMin;
  octree->zMax = zMax;
  octree->nOctant=1;
  octree->octant = malloc(sizeof(Octant));
  
  octree->octant->data[0]=1.0;
  octree->octant->data[1]=0.0;
  octree->octant->data[2]=0.0;

  octree->octant->data[3]=0.0;
  octree->octant->data[4]=1.0;
  octree->octant->data[5]=0.0;

  octree->octant->data[6]=0.0;
  octree->octant->data[7]=0.0;
  octree->octant->data[8]=1.0;

  return octree;
}

void octreeFree( Octree *octree )
{
  free( octree );
}

Octree *octreeBoundingBox( Octree *octree, double *boundingBox )
{
  boundingBox[0] = octree->xMin;
  boundingBox[1] = octree->xMax;
  boundingBox[2] = octree->yMin;
  boundingBox[3] = octree->yMax;
  boundingBox[4] = octree->zMin;
  boundingBox[5] = octree->zMax;
  return octree;
}

int octreeNOctant( Octree *octree )
{
  return octree->nOctant;
}

Octree *octreeInsert( Octree *octree, double *location, double *data )
{
  int i;

  for (i=0;i<9;i++) octree->octant->data[i]=data[i];

  return octree;
}

Octree *octreeQuery( Octree *octree, double *location, double *data )
{
  int i;

  for (i=0;i<9;i++) data[i] = octree->octant->data[i];

  return octree;
}

