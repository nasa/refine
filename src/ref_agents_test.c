
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

#include "ref_agents.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_mpi.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  { /* create */
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents, ref_mpi), "make agents");
    RSS(ref_agents_free(ref_agents), "agents free");
  }

  { /* push sets home */
    REF_INT node = 0, part = 0, seed = 0, id;
    REF_DBL xyz[] = {1.0, 2.0, 3.0};
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents, ref_mpi), "make agents");
    RSS(ref_agents_push(ref_agents, node, part, seed, xyz, &id), "add");
    REIS(ref_mpi_rank(ref_mpi), ref_agents->agent[0].home, "home not set");
    RSS(ref_agents_free(ref_agents), "agents free");
  }

  { /* remove middle, pop all */
    REF_INT node, part = 0, seed = 0, id;
    REF_DBL xyz[] = {1.0, 2.0, 3.0};
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents, ref_mpi), "create");

    RSS(ref_agents_push(ref_agents, 10, part, seed, xyz, &id), "add");
    RSS(ref_agents_push(ref_agents, 11, part, seed, xyz, &id), "add");
    RSS(ref_agents_push(ref_agents, 12, part, seed, xyz, &id), "add");

    RSS(ref_agents_remove(ref_agents, 1), "remove middle");

    RSS(ref_agents_pop(ref_agents, &node, &part, &seed, xyz), "pop");
    REIS(12, node, "wrong node");
    RSS(ref_agents_pop(ref_agents, &node, &part, &seed, xyz), "pop");
    REIS(10, node, "wrong node");

    RSS(ref_agents_free(ref_agents), "free");
  }

  { /* delete first, last, middle */
    REF_AGENTS ref_agents;
    REF_INT part = 0, seed = 0, id;
    REF_DBL xyz[] = {1.0, 2.0, 3.0};
    RSS(ref_agents_create(&ref_agents, ref_mpi), "create");

    RSS(ref_agents_push(ref_agents, 10, part, seed, xyz, &id), "add");
    RSS(ref_agents_push(ref_agents, 11, part, seed, xyz, &id), "add");
    RSS(ref_agents_push(ref_agents, 12, part, seed, xyz, &id), "add");
    REIS(3, ref_agents_n(ref_agents), "wrong count");

    RSS(ref_agents_delete(ref_agents, 10), "del first");
    RSS(ref_agents_delete(ref_agents, 12), "del last");
    RSS(ref_agents_delete(ref_agents, 11), "del middle");
    REIS(0, ref_agents_n(ref_agents), "wrong count");

    RSS(ref_agents_free(ref_agents), "free");
  }

  { /* remove max id */
    REF_INT node, max, part = 0, seed = 0, id;
    REF_DBL xyz[] = {1.0, 2.0, 3.0};
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents, ref_mpi), "create");

    max = ref_agents->max;
    for (node = 0; node < max; node++)
      RSS(ref_agents_push(ref_agents, node, part, seed, xyz, &id), "add");

    RSS(ref_agents_remove(ref_agents, max - 1), "remove last agents");

    RSS(ref_agents_free(ref_agents), "free");
  }

  { /* add bunch testing realloc */
    REF_INT node, max, part = 0, seed = 0, id;
    REF_DBL xyz[] = {1.0, 2.0, 3.0};
    REF_AGENTS ref_agents;
    RSS(ref_agents_create(&ref_agents, ref_mpi), "create");

    max = ref_agents->max;
    for (node = 0; node < 30000; node++)
      RSS(ref_agents_push(ref_agents, node, part, seed, xyz, &id), "add");

    RAS(max < ref_agents->max, "grow max");

    RSS(ref_agents_free(ref_agents), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
