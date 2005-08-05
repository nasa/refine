
/* Ring, continuous loop of linear segements on a geometry face that
 *       is meshed by filling with triangles.
 *
 * Michael A. Park
 * Computational AeroSciences Branch
 * NASA Langley Research Center
 * Phone: (757) 864-6604
 * Email: Mike.Park@NASA.Gov
 */

/* $Id$ */

#ifndef RING_H
#define RING_H

#include "refine_defs.h"
#include <stdio.h>

BEGIN_C_DECLORATION

typedef struct Ring Ring;

struct Ring {
  int segments;
  int triangles;
};

Ring *ringCreate( void );
void ringFree( Ring * );

int ringSegments( Ring * );
int ringTriangles( Ring * );

END_C_DECLORATION

#endif /* RING_H */
