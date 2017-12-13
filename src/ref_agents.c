
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
#include <stdio.h>
#include <math.h>

#include "ref_agents.h"

#include "ref_malloc.h"

/* REF_EMPTY is terminatior, next avalable is shifted by 2*/
#define next2index(next) (-(next)-2)
#define index2next(index) (-2-(index))

REF_STATUS ref_agents_create( REF_AGENTS *ref_agents_ptr, REF_MPI ref_mpi )
{
  REF_AGENTS ref_agents;
  REF_INT id;
  
  ref_malloc( *ref_agents_ptr, 1, REF_AGENTS_STRUCT );
  ref_agents = ( *ref_agents_ptr );

  ref_agents->ref_mpi = ref_mpi;

  ref_agents->n = 0;
  ref_agents->max = 10;

  ref_malloc( ref_agents->agent, ref_agents->max, REF_AGENT_STRUCT );

  for (id=0;id<ref_agents->max;id++)
    ref_agents->agent[id].node = index2next(id+1);
  ref_agents->agent[(ref_agents->max)-1].node = REF_EMPTY;
  ref_agents->blank = index2next(0);

  return REF_SUCCESS;
}

REF_STATUS ref_agents_free( REF_AGENTS ref_agents )
{
  if ( NULL == (void *)ref_agents )
    return REF_NULL;
  ref_free( ref_agents->agent );
  ref_free( ref_agents );
  return REF_SUCCESS;
}

