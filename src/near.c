
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

struct Near {
  int index;
  double x, y, z;
  double radius;
  Near *rightChild, *leftChild;
  double farChild, rightDistance, leftDistance;
};

Near* nearCreate( int index, double x, double y, double z, double radius )
{
  Near *near;

  near = malloc( sizeof(Near) );

  near->index = index;
  near->x = x;
  near->y = y;
  near->z = z;
  near->radius = radius;

  near->rightChild = NULL;
  near->leftChild = NULL;

  near->farChild = 0;
  near->rightDistance = 0;
  near->leftDistance = 0;

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
  if (NULL==near->leftChild){ 
    near->leftChild = child;
    return near;
  }
  if (NULL==near->rightChild){ 
    near->rightChild = child;
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
  clearance = MAX( 0, clearance );

  return clearance;
}

double nearFarChild( Near *near )
{
  return near->farChild;
}

double nearRightDistance( Near *near )
{
  return near->rightDistance;
}

double nearLeftDistance( Near *near )
{
  return near->leftDistance;
}
