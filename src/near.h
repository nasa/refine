
/* Node in a near tree with cached child distance
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef NEAR_H
#define NEAR_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Near Near;

struct Near {
  int index;
  double x, y, z;
  double radius;
  Near *leftChild, *rightChild;
  double leftRadius, rightRadius;
};

Near *nearCreate(int index, double x, double y, double z, double radius );
Near *nearInit(Near *, int index, double x, double y, double z, double radius );
void nearFree( Near * );
int nearIndex( Near * );
int nearLeftIndex( Near * );
int nearRightIndex( Near * );
double nearDistance( Near *, Near *other);
double nearClearance( Near *, Near *other);
Near *nearInsert( Near *, Near *child );
double nearLeftRadius( Near * );
double nearRightRadius( Near * );

Near *nearVisualize( Near * );

int nearCollisions(Near *, Near *target);
Near *nearTouched(Near *, Near *target, int *found, int maxfound, int *list);

int nearNearestIndex(Near *root, Near *key);
Near *nearNearestIndexAndDistance(Near *root, Near *key, 
				  int *index, double *distance);

END_C_DECLORATION

#endif /* NEAR_H */
