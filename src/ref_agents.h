

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
typedef REF_AGENTS_STRUCT * REF_AGENTS;
typedef struct REF_AGENT_STRUCT REF_AGENT_STRUCT;
typedef REF_AGENT_STRUCT * REF_AGENT;
END_C_DECLORATION

#include "ref_mpi.h"

BEGIN_C_DECLORATION
struct REF_AGENT_STRUCT {
  REF_INT node;
  REF_INT guess;
};
struct REF_AGENTS_STRUCT {
  REF_INT n, max;
  REF_INT blank;
  REF_AGENT_STRUCT *agent;
  REF_MPI ref_mpi;
};

REF_STATUS ref_agents_create( REF_AGENTS *ref_agents, REF_MPI ref_mpi );

REF_STATUS ref_agents_free( REF_AGENTS ref_agents );

#endif /* REF_AGENTS_H */
