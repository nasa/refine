
/* Plan, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#include <stdlib.h>
#ifndef __APPLE__       /* Not needed on Mac OS X */
#include <malloc.h>
#endif
#include "plan.h"
#include "sort.h"

Plan* planCreate( int max_size, int chunk_size )
{
  Plan *plan;

  plan = (Plan *)malloc( sizeof(Plan) );

  plan->size       = 0;
  plan->max_size   = MAX(max_size,1);
  plan->chunk_size = MAX(chunk_size,1);

  plan->item     = (int *)   malloc( planMaxSize(plan)*sizeof(int) );
  plan->priority = (double *)malloc( planMaxSize(plan)*sizeof(double) );
  plan->ranking  = NULL;

  return plan;
}

void planFree( Plan *plan )
{
  if (NULL!=plan->ranking)  free(plan->ranking);
  if (NULL!=plan->priority) free(plan->priority);
  if (NULL!=plan->item)     free(plan->item);
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

Plan *planAddItemWithPriority( Plan *plan, int item, double priority )
{
  int newitem;
  if (planSize(plan)>=planMaxSize(plan)){
    plan->max_size += plan->chunk_size;
    plan->item     = (int *)   realloc(plan->item, 
				       planMaxSize(plan)*sizeof(int) );
    plan->priority = (double *)realloc(plan->priority, 
				       planMaxSize(plan)*sizeof(double) );
  }
  if (NULL!=plan->ranking) {
    free(plan->ranking);
    plan->ranking = NULL;
  }

  newitem = planSize(plan);
  plan->size++;

  plan->item[newitem] = item;
  plan->priority[newitem] = priority;

  return plan;
}

Plan *planDeriveRankingsFromPriorities( Plan *plan )
{
  int i;
  if (0==planSize(plan)) return NULL;
  plan->ranking = (int *)malloc( planSize(plan)*sizeof(int) );
  for(i=0;i<planSize(plan);i++) plan->ranking[i] = EMPTY;
  sortDoubleHeap(planSize(plan),plan->priority,plan->ranking);
  return plan;
}

int planItemWithThisRanking( Plan *plan, int ranking )
{
  if (NULL==plan->ranking) return EMPTY;
  if (ranking<0||ranking>=planSize(plan)) return EMPTY;
  return plan->item[plan->ranking[ranking]];
}

