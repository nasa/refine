
/* Octree data structure for fast geometric searches
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef OCTREE_H
#define OCTREE_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Octree Octree;

Octree *octreeCreate( double xMin, double xMax, 
		      double yMin, double yMax, 
		      double zMin, double zMax );
void octreeFree( Octree * );
Octree *octreeBoundingBox( Octree *, double *boundingBox );
int octreeNOctant( Octree * );

END_C_DECLORATION

#endif /* OCTREE_H */
