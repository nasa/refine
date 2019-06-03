

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

#ifndef REF_AGENTS_H
#define REF_AGENTS_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_AGENTS_STRUCT REF_AGENTS_STRUCT;
typedef REF_AGENTS_STRUCT *REF_AGENTS;
typedef struct REF_AGENT_STRUCT REF_AGENT_STRUCT;
typedef REF_AGENT_STRUCT *REF_AGENT;
typedef enum REF_AGENT_MODES { /* 0 */ REF_AGENT_UNUSED,
                               /* 1 */ REF_AGENT_WALKING,
                               /* 2 */ REF_AGENT_AT_BOUNDARY,
                               /* 3 */ REF_AGENT_HOP_PART,
                               /* 4 */ REF_AGENT_SUGGESTION,
                               /* 5 */ REF_AGENT_ENCLOSING,
                               /* 6 */ REF_AGENT_TERMINATED,
                               /* 7 */ REF_AGENT_NEW,
                               /* 8 */ REF_AGENT_MODE_LAST } REF_AGENT_MODE;
END_C_DECLORATION

#include "ref_mpi.h"

BEGIN_C_DECLORATION
struct REF_AGENT_STRUCT {
  REF_AGENT_MODE mode;
  REF_INT previous; /* agent list navigation */
  REF_INT next;

  REF_INT home;    /* mpi rank of the to node that needs an interpolant */
  REF_INT node;    /* the local node id that needs an interpolant
                    * empty, map global when SUGGESTION */
  REF_INT part;    /* mpi rank of the seed */
  REF_INT seed;    /* cell guess when WALKING or BOUNDARY
                    * empty when HOP_PART
                    * from cell when ENCLOSE or SUGGESTION */
  REF_INT global;  /* global node guess when HOP_PART
                    * global node that needs an interpolant when SUGGESTION
                    * empty otherwise */
  REF_INT step;    /* number of cells visited */
  REF_DBL xyz[3];  /* the to xyz that needs an interpolant */
  REF_DBL bary[4]; /* the from bary of the from cell when ENCLOSE */
};
struct REF_AGENTS_STRUCT {
  REF_INT n, max;
  REF_INT blank;
  REF_INT last;
  REF_AGENT_STRUCT *agent;
  REF_MPI ref_mpi;
};

#define ref_agents_n(ref_agents) ((ref_agents)->n)
#define ref_agents_max(ref_agents) ((ref_agents)->max)

#define ref_agent_mode(ref_agents, id) ((ref_agents)->agent[(id)].mode)
#define ref_agent_valid(ref_agents, id) \
  (REF_AGENT_UNUSED != ref_agent_mode(ref_agents, id))

#define ref_agent_home(ref_agents, id) ((ref_agents)->agent[(id)].home)
#define ref_agent_node(ref_agents, id) ((ref_agents)->agent[(id)].node)
#define ref_agent_part(ref_agents, id) ((ref_agents)->agent[(id)].part)
#define ref_agent_seed(ref_agents, id) ((ref_agents)->agent[(id)].seed)
#define ref_agent_global(ref_agents, id) ((ref_agents)->agent[(id)].global)
#define ref_agent_step(ref_agents, id) ((ref_agents)->agent[(id)].step)

#define ref_agent_xyz(ref_agents, j, id) ((ref_agents)->agent[(id)].xyz[j])
#define ref_agent_xyz_ptr(ref_agents, id) ((ref_agents)->agent[(id)].xyz)
#define ref_agent_bary(ref_agents, j, id) ((ref_agents)->agent[(id)].bary[j])

#define each_active_ref_agent(ref_agents, id)               \
  for ((id) = 0; (id) < ref_agents_max(ref_agents); (id)++) \
    if (ref_agent_valid(ref_agents, id))

#define each_ref_agent_step(ref_agents, id, limit) \
  for (; ref_agent_step(ref_agents, id) < (limit); \
       ref_agent_step(ref_agents, id)++)

REF_STATUS ref_agents_create(REF_AGENTS *ref_agents, REF_MPI ref_mpi);

REF_STATUS ref_agents_free(REF_AGENTS ref_agents);

REF_STATUS ref_agents_inspect(REF_AGENTS ref_agents);
REF_STATUS ref_agents_tattle(REF_AGENTS ref_agents, REF_INT id,
                             const char *context);
REF_STATUS ref_agents_population(REF_AGENTS ref_agents, const char *context);

REF_STATUS ref_agents_restart(REF_AGENTS ref_agents, REF_INT part, REF_INT seed,
                              REF_INT id);
REF_STATUS ref_agents_push(REF_AGENTS ref_agents, REF_INT node, REF_INT part,
                           REF_INT seed, REF_DBL *xyz, REF_INT *id);
REF_STATUS ref_agents_remove(REF_AGENTS ref_agents, REF_INT id);
REF_STATUS ref_agents_pop(REF_AGENTS ref_agents, REF_INT *node, REF_INT *part,
                          REF_INT *seed, REF_DBL *xyz);
REF_STATUS ref_agents_delete(REF_AGENTS ref_agents, REF_INT node);

REF_STATUS ref_agents_migrate(REF_AGENTS ref_agents);

END_C_DECLORATION

#endif /* REF_AGENTS_H */
