
/* Plan, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#ifndef PLAN_H
#define PLAN_H

#include "refine_defs.h"
#include <stdio.h>

BEGIN_C_DECLORATION

typedef struct Plan Plan;

struct Plan {
  int size;
  int max_size;
  int chunk_size;
};

Plan *planCreate( int initial_size, int chunk_size );
void planFree( Plan * );

int planSize( Plan * );
int planMaxSize( Plan * );
int planChunkSize( Plan * );

END_C_DECLORATION

#endif /* PLAN_H */
