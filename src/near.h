
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

Near *nearCreate(int index, double x, double y, double z, double radius );
void adjFree( Near *near );
int nearIndex( Near *near );
int nearRightIndex( Near *near );
int nearLeftIndex( Near *near );
double nearDistance( Near *near, Near *other);
double nearClearance( Near *near, Near *other);
Near *nearInsert( Near *near, Near *child );
double nearFarChild( Near *near );
double nearRightDistance( Near *near );
double nearLeftDistance( Near *near );

int nearCollisions(Near *near, Near *target);
Near *nearTouched(Near *near, Near *target, int *found, int maxfound, int *list);

END_C_DECLORATION

#endif /* NEAR_H */
