
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
Near *nearInit(Near *, int index, double x, double y, double z, double radius );
void nearFree( Near * );
int nearIndex( Near * );
int nearRightIndex( Near * );
int nearLeftIndex( Near * );
double nearDistance( Near *, Near *other);
double nearClearance( Near *, Near *other);
Near *nearInsert( Near *, Near *child );
double nearFarChild( Near * );
double nearRightDistance( Near * );
double nearLeftDistance( Near * );

int nearCollisions(Near *, Near *target);
Near *nearTouched(Near *, Near *target, int *found, int maxfound, int *list);

END_C_DECLORATION

#endif /* NEAR_H */
