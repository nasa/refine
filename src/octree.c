
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

struct Octree {
  double xMin, xMax, yMin, yMax, zMin, zMax;
  int nOctant;
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
  octree->nOctant=0;
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
