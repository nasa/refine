
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
  int *item;
  double *priority;
  int *ranking;
};

Plan *planCreate( int initial_size, int chunk_size );
void planFree( Plan * );

int planSize( Plan * );
int planMaxSize( Plan * );
int planChunkSize( Plan * );

Plan *planAddItemWithPriority( Plan *, int item, double priority );
Plan *planDeriveRankingsFromPriorities( Plan * );
int planItemWithThisRanking( Plan *, int ranking );
double planPriorityWithThisRanking( Plan *, int ranking );

END_C_DECLORATION

#endif /* PLAN_H */
