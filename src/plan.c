
/* Plan, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#include <stdlib.h>
#include "plan.h"


Plan* planCreate( int max_size, int chunk_size )
{
  Plan *plan;

  plan = malloc( sizeof(Plan) );

  plan->size       = 0;
  plan->max_size   = MAX(max_size,1);
  plan->chunk_size = MAX(chunk_size,1);

  return plan;
}

void planFree( Plan *plan )
{
  free( plan );
}

int planSize( Plan *plan )
{
  return plan->size;
}

int planMaxSize( Plan *plan )
{
  return plan->max_size;
}

int planChunkSize( Plan *plan )
{
  return plan->chunk_size;
}

