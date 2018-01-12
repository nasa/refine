
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

#include <stdlib.h>

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

double planPriorityWithThisRanking( Plan *plan, int ranking )
{
  if (NULL==plan->ranking) return -999.0;
  if (ranking<0||ranking>=planSize(plan)) return -999.0;
  return plan->priority[plan->ranking[ranking]];
}

