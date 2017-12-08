
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
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
