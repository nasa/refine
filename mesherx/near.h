
/* Node in a near tree with cached child distance
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#ifndef NEAR_H
#define NEAR_H

#include "refine_defs.h"
#include <math.h>

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
#define nearDistance(near,other) \
  ( sqrt( (near->x - other->x)*(near->x - other->x) + \
          (near->y - other->y)*(near->y - other->y) + \
          (near->z - other->z)*(near->z - other->z) ) )
double nearClearance( Near *, Near *other);
Near *nearInsert( Near *, Near *child );
#define nearLeftRadius(near) (near->leftRadius)
#define nearRightRadius(near) (near->rightRadius)

Near *nearVisualize( Near * );

int nearCollisions(Near *, Near *target);
Near *nearTouched(Near *, Near *target, int *found, int maxfound, int *list);

int nearNearestIndex(Near *root, Near *key);
Near *nearNearestIndexAndDistance(Near *root, Near *key, 
				  int *index, double *distance);

END_C_DECLORATION

#endif /* NEAR_H */
