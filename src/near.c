
/* Node in a near tree with cached child distance
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <math.h>
#include "near.h"

Near* nearCreate( int index, double x, double y, double z, double radius )
{
  Near *near;

  near = malloc( sizeof(Near) );

  return nearInit( near, index, x, y, z, radius );
}

Near* nearInit( Near *near, 
		int index, double x, double y, double z, double radius )
{
  near->index = index;
  near->x = x;
  near->y = y;
  near->z = z;
  near->radius = radius;

  near->rightChild = NULL;
  near->leftChild = NULL;

  near->farChild = 0;

  return near;
}

void nearFree( Near *near )
{
  free( near );
}

int nearIndex( Near *near )
{
  return (NULL==near?EMPTY:near->index);
}

int nearRightIndex( Near *near )
{
  return nearIndex(near->rightChild);
}

int nearLeftIndex( Near *near )
{
  return nearIndex(near->leftChild);
}

Near *nearInsert( Near *near, Near *child )
{
  double distance;

  distance = nearDistance(near,child) + child->radius;

  near->farChild = MAX(near->farChild,distance);

  if (NULL==near->leftChild){ 
    near->leftChild = child;
    child->farChild = distance;
    return near;
  }
  if (NULL==near->rightChild){ 
    near->rightChild = child;
    child->farChild = distance;
    return near;
  }
  if ( nearDistance(near->leftChild,child)
     < nearDistance(near->rightChild,child) ) {
    return (near->leftChild==nearInsert(near->leftChild,child)?near:NULL);
  }else{
    return (near->rightChild==nearInsert(near->rightChild,child)?near:NULL);
  }
}

double nearDistance( Near *near, Near *other)
{
  double dx, dy, dz, distance;

  dx = near->x - other->x;
  dy = near->y - other->y;
  dz = near->z - other->z;

  distance  = sqrt( dx*dx + dy*dy + dz*dz );
  return distance;
}

double nearClearance( Near *near, Near *other)
{
  double clearance;

  clearance = nearDistance(near,other) - near->radius - other->radius;

  return clearance;
}

double nearFarChild( Near *near )
{
  return (NULL==near?0:near->farChild);
}

double nearRightDistance( Near *near )
{
  return nearFarChild(near->rightChild);
}

double nearLeftDistance( Near *near )
{
  return nearFarChild(near->leftChild);
}

int nearCollisions(Near *near, Near *target)
{
  int collisions = 0;

  if (NULL==near || NULL == target) return 0;

  collisions += nearCollisions(near->rightChild,target);
  collisions += nearCollisions(near->leftChild,target);

  if (nearClearance(near,target) <= 0) collisions++;

  return collisions;
}

Near *nearTouched(Near *near, Near *target, int *found, int maxfound, int *list)
{
  double distance, safeZone;

  if (NULL==near || NULL==target) return NULL;

  distance = nearDistance( near, target);
  safeZone = distance - target->radius;

  if (safeZone <= nearRightDistance(near) ) 
    nearTouched(near->rightChild, target, found, maxfound, list);
  if (safeZone <= nearLeftDistance(near) ) 
    nearTouched(near->leftChild, target, found, maxfound, list);

  if ( (safeZone - near->radius) <= 0) {
    if (*found >= maxfound) return NULL;
    list[*found] = near->index;
    (*found)++;
  }

  return near;
}
