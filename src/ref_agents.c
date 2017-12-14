
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

#define ref_agent_previous(ref_agents,id) ((ref_agents)->agent[(id)].previous)
#define ref_agent_next(ref_agents,id) ((ref_agents)->agent[(id)].next)

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
    {
      ref_agent_mode(ref_agents,id) = REF_AGENT_UNUSED;
      ref_agent_previous(ref_agents,id) = REF_EMPTY;
      ref_agent_next(ref_agents,id) = id+1;
    }
  ref_agent_next(ref_agents,ref_agents_max(ref_agents)-1) = REF_EMPTY;
  ref_agents->blank = 0;
  ref_agents->last = REF_EMPTY;

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

REF_STATUS ref_agents_inspect( REF_AGENTS ref_agents )
{
  REF_INT id;
  for ( id = 0 ; id < ref_agents->max ; id++ )
    {
      if ( REF_AGENT_UNUSED != ref_agent_mode(ref_agents,id))
	printf("%2d: %2d mode %2d prev %2d next\n",
	       id,
	       ref_agent_mode(ref_agents,id),
	       ref_agent_previous(ref_agents,id),
	       ref_agent_next(ref_agents,id));
    }
  printf("%d n %d max %d\nblank %d last\n",
	 ref_agents->n, ref_agents->max, ref_agents->blank, ref_agents->last );
  printf("%d REF_AGENT_UNUSED\n",REF_AGENT_UNUSED );
  printf("%d REF_AGENT_WALKING\n",REF_AGENT_WALKING );
  printf("%d REF_AGENT_AT_BOUNDARY\n",REF_AGENT_AT_BOUNDARY );
  printf("%d REF_AGENT_HOP_PART\n",REF_AGENT_HOP_PART );
  printf("%d REF_AGENT_ENCLOSING\n",REF_AGENT_ENCLOSING );
  printf("%d REF_AGENT_MODE_LAST\n",REF_AGENT_MODE_LAST );

  return REF_SUCCESS;
}

REF_STATUS ref_agents_push( REF_AGENTS ref_agents, 
			    REF_INT node, REF_INT part,
			    REF_INT seed, REF_DBL *xyz )
{
  REF_INT i, id;

  if ( REF_EMPTY == ref_agents->blank )
    {
      REF_INT orig, chunk, extra;
      orig = ref_agents->max;
      chunk = MAX(5000,(REF_INT)(1.5*(REF_DBL)orig));
      ref_agents->max = orig + chunk;
      ref_realloc( ref_agents->agent, ref_agents->max, REF_AGENT_STRUCT);
      
      for (extra=orig;extra < ref_agents->max; extra++ ) 
	{
	  ref_agent_mode(ref_agents,extra) = REF_AGENT_UNUSED;
	  ref_agent_previous(ref_agents,extra) = REF_EMPTY;
	  ref_agent_next(ref_agents,extra) = extra+1;
	}
      ref_agent_next(ref_agents,ref_agents_max(ref_agents)-1) = REF_EMPTY;
      ref_agents->blank = orig;
    }

  id = ref_agents->blank;
  ref_agents->blank = ref_agent_next(ref_agents,id);

  if ( REF_EMPTY != ref_agents->last )
    ref_agent_next(ref_agents,ref_agents->last) = id;

  ref_agent_mode(ref_agents,id) = REF_AGENT_WALKING;
  ref_agent_previous(ref_agents,id) = ref_agents->last;
  ref_agent_next(ref_agents,id) = REF_EMPTY;

  ref_agent_home(ref_agents,id) = ref_mpi_rank(ref_agents->ref_mpi);
  ref_agent_node(ref_agents,id) = node;
  ref_agent_part(ref_agents,id) = part;
  ref_agent_seed(ref_agents,id) = seed;
  ref_agent_step(ref_agents,id) = 0;
  for (i=0;i<3;i++)
    ref_agent_xyz(ref_agents,i,id) = xyz[i];
 
  ref_agents->last = id;

  ref_agents_n(ref_agents)++;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_remove( REF_AGENTS ref_agents, REF_INT id )
{

  if ( id < 0 || ref_agents_max(ref_agents) <= id )
    return REF_INVALID;
  if ( !ref_agent_valid(ref_agents,id) )
    return REF_INVALID;

  if ( ref_agents->last == id )
    ref_agents->last = ref_agent_previous(ref_agents,id);

  if ( REF_EMPTY != ref_agent_previous(ref_agents,id) )
    ref_agent_next(ref_agents,ref_agent_previous(ref_agents,id)) = 
      ref_agent_next(ref_agents,id);

  if ( REF_EMPTY != ref_agent_next(ref_agents,id) )
    ref_agent_previous(ref_agents,ref_agent_next(ref_agents,id)) =
      ref_agent_previous(ref_agents,id);

  ref_agent_mode(ref_agents,id) = REF_AGENT_UNUSED;
  ref_agent_previous(ref_agents,id) = REF_EMPTY;
  ref_agent_next(ref_agents,id) = ref_agents->blank;
  ref_agents->blank = id;

  ref_agents_n(ref_agents)--;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_pop( REF_AGENTS ref_agents, 
			   REF_INT *node, REF_INT *part,
			   REF_INT *seed, REF_DBL *xyz )
{
  REF_INT i, id;

  if ( REF_EMPTY == ref_agents->last )
    return REF_NOT_FOUND;

  id = ref_agents->last;

  *node = ref_agent_node(ref_agents,id);
  *part = ref_agent_part(ref_agents,id);
  *seed = ref_agent_seed(ref_agents,id);
  for (i=0;i<3;i++)
    xyz[i] = ref_agent_xyz(ref_agents,i,id);

  RSS( ref_agents_remove( ref_agents, id ), "rm" );

  return REF_SUCCESS;
}

REF_STATUS ref_agents_delete( REF_AGENTS ref_agents, REF_INT node )
{
  REF_INT id;

  if ( node < 0 )
    return REF_INVALID;

  each_active_ref_agent( ref_agents, id )
    {
      if ( node == ref_agent_node(ref_agents,id) )
	RSS( ref_agents_remove( ref_agents, id ), "rm" );
    }

  return REF_SUCCESS;
}

