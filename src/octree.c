
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
  bool hasChildren;
  double data[9];
  Octant *o000;
  Octant *o001;
  Octant *o010;
  Octant *o011;
  Octant *o100;
  Octant *o101;
  Octant *o110;
  Octant *o111;
};

struct Octree {
  double xMin, xMax, yMin, yMax, zMin, zMax;
  int nOctant;
  int maxOctant;
  Octant *octant;
};

Octant *octantInit(Octant *octant)
{
  octant->hasChildren = FALSE;

  octant->data[0]=1.0;
  octant->data[1]=0.0;
  octant->data[2]=0.0;

  octant->data[3]=0.0;
  octant->data[4]=1.0;
  octant->data[5]=0.0;

  octant->data[6]=0.0;
  octant->data[7]=0.0;
  octant->data[8]=1.0;

  octant->o000=NULL;
  octant->o001=NULL;
  octant->o010=NULL;
  octant->o011=NULL;
  octant->o100=NULL;
  octant->o101=NULL;
  octant->o110=NULL;
  octant->o111=NULL;
}

Octree* octreeCreate( double xMin, double xMax, 
		      double yMin, double yMax, 
		      double zMin, double zMax )
{
  int i;
  Octree *octree;

  octree = malloc( sizeof(Octree) );
  octree->xMin = xMin;
  octree->xMax = xMax;
  octree->yMin = yMin;
  octree->yMax = yMax;
  octree->zMin = zMin;
  octree->zMax = zMax;

  octree->nOctant=1;
  octree->maxOctant=5000;

  octree->octant = malloc(octree->maxOctant*sizeof(Octant));
  
  for (i=0;i<octree->maxOctant;i++) octantInit(&octree->octant[i]);

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

