
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_agents.h"

#include "ref_malloc.h"

#define ref_agent_previous(ref_agents, id) ((ref_agents)->agent[(id)].previous)
#define ref_agent_next(ref_agents, id) ((ref_agents)->agent[(id)].next)

REF_STATUS ref_agents_create(REF_AGENTS *ref_agents_ptr, REF_MPI ref_mpi) {
  REF_AGENTS ref_agents;
  REF_INT id;

  ref_malloc(*ref_agents_ptr, 1, REF_AGENTS_STRUCT);
  ref_agents = (*ref_agents_ptr);

  ref_agents->ref_mpi = ref_mpi;

  ref_agents->n = 0;
  ref_agents->max = 10;

  ref_malloc(ref_agents->agent, ref_agents->max, REF_AGENT_STRUCT);

  for (id = 0; id < ref_agents->max; id++) {
    ref_agent_mode(ref_agents, id) = REF_AGENT_UNUSED;
    ref_agent_previous(ref_agents, id) = REF_EMPTY;
    ref_agent_next(ref_agents, id) = id + 1;
  }
  ref_agent_next(ref_agents, ref_agents_max(ref_agents) - 1) = REF_EMPTY;
  ref_agents->blank = 0;
  ref_agents->last = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_free(REF_AGENTS ref_agents) {
  if (NULL == (void *)ref_agents) return REF_NULL;
  ref_free(ref_agents->agent);
  ref_free(ref_agents);
  return REF_SUCCESS;
}

REF_STATUS ref_agents_inspect(REF_AGENTS ref_agents) {
  REF_INT id;
  for (id = 0; id < ref_agents->max; id++) {
    if (REF_AGENT_UNUSED != ref_agent_mode(ref_agents, id))
      printf("%2d: %2d mode %2d prev %2d next\n", id,
             (int)ref_agent_mode(ref_agents, id),
             ref_agent_previous(ref_agents, id),
             ref_agent_next(ref_agents, id));
  }
  printf("%d n %d max %d\nblank %d last\n", ref_agents->n, ref_agents->max,
         ref_agents->blank, ref_agents->last);

  return REF_SUCCESS;
}

REF_STATUS ref_agents_tattle(REF_AGENTS ref_agents, REF_INT id,
                             const char *context) {
  printf(
      "%d: %d id %d mode %d home %d node %d part %d seed %d global %f %f %f "
      "%s\n",
      ref_mpi_rank(ref_agents->ref_mpi), id,
      (int)ref_agent_mode(ref_agents, id), ref_agent_home(ref_agents, id),
      ref_agent_node(ref_agents, id), ref_agent_part(ref_agents, id),
      ref_agent_seed(ref_agents, id), ref_agent_global(ref_agents, id),
      ref_agent_xyz(ref_agents, 1, id), ref_agent_xyz(ref_agents, 1, id),
      ref_agent_xyz(ref_agents, 2, id), context);
  return REF_SUCCESS;
}

REF_STATUS ref_agents_population(REF_AGENTS ref_agents, const char *context) {
  REF_MPI ref_mpi = ref_agents->ref_mpi;
  REF_INT id, *counts;
  ref_malloc_init(counts, REF_AGENT_MODE_LAST, REF_INT, 0);
  for (id = 0; id < ref_agents->max; id++) {
    RAS(0 <= ref_agent_mode(ref_agents, id), "last");
    RAS(REF_AGENT_MODE_LAST > ref_agent_mode(ref_agents, id), "last");
    counts[ref_agent_mode(ref_agents, id)]++;
  }
  RSS(ref_mpi_allsum(ref_mpi, counts, REF_AGENT_MODE_LAST, REF_INT_TYPE), "as");
  if (ref_mpi_once(ref_mpi)) {
    REF_INT total = 0;
    for (id = 1; id < REF_AGENT_MODE_LAST; id++) total += counts[id];
    printf(" %5d of", total);
    for (id = 1; id < REF_AGENT_MODE_LAST; id++)
      printf(" %d:%4d", id, counts[id]);
    printf(" %s\n", context);
  }
  ref_free(counts);
  return REF_SUCCESS;
}

REF_STATUS ref_agents_new(REF_AGENTS ref_agents, REF_INT *id) {
  if (REF_EMPTY == ref_agents->blank) {
    REF_INT orig, chunk, extra;
    orig = ref_agents->max;
    chunk = MAX(5000, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_agents->max = orig + chunk;
    ref_realloc(ref_agents->agent, ref_agents->max, REF_AGENT_STRUCT);

    for (extra = orig; extra < ref_agents->max; extra++) {
      ref_agent_mode(ref_agents, extra) = REF_AGENT_UNUSED;
      ref_agent_previous(ref_agents, extra) = REF_EMPTY;
      ref_agent_next(ref_agents, extra) = extra + 1;
    }
    ref_agent_next(ref_agents, ref_agents_max(ref_agents) - 1) = REF_EMPTY;
    ref_agents->blank = orig;
  }

  *id = ref_agents->blank;
  ref_agents->blank = ref_agent_next(ref_agents, *id);

  if (REF_EMPTY != ref_agents->last)
    ref_agent_next(ref_agents, ref_agents->last) = *id;

  ref_agent_mode(ref_agents, *id) = REF_AGENT_NEW;
  ref_agent_previous(ref_agents, *id) = ref_agents->last;
  ref_agent_next(ref_agents, *id) = REF_EMPTY;

  ref_agents->last = *id;

  ref_agents_n(ref_agents)++;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_restart(REF_AGENTS ref_agents, REF_INT part, REF_INT seed,
                              REF_INT id) {
  ref_agent_mode(ref_agents, id) = REF_AGENT_WALKING;
  ref_agent_part(ref_agents, id) = part;
  ref_agent_seed(ref_agents, id) = seed;
  ref_agent_global(ref_agents, id) = REF_EMPTY;
  ref_agent_step(ref_agents, id) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_push(REF_AGENTS ref_agents, REF_INT node, REF_INT part,
                           REF_INT seed, REF_DBL *xyz, REF_INT *id_ptr) {
  REF_INT i, id;

  RSS(ref_agents_new(ref_agents, &id), "new");
  *id_ptr = id;

  ref_agent_mode(ref_agents, id) = REF_AGENT_WALKING;
  ref_agent_home(ref_agents, id) = ref_mpi_rank(ref_agents->ref_mpi);
  ref_agent_node(ref_agents, id) = node;
  ref_agent_part(ref_agents, id) = part;
  ref_agent_seed(ref_agents, id) = seed;
  ref_agent_global(ref_agents, id) = REF_EMPTY;
  ref_agent_step(ref_agents, id) = 0;
  for (i = 0; i < 3; i++) ref_agent_xyz(ref_agents, i, id) = xyz[i];

  return REF_SUCCESS;
}

REF_STATUS ref_agents_remove(REF_AGENTS ref_agents, REF_INT id) {
  if (id < 0 || ref_agents_max(ref_agents) <= id) return REF_INVALID;
  if (!ref_agent_valid(ref_agents, id)) return REF_INVALID;

  if (ref_agents->last == id)
    ref_agents->last = ref_agent_previous(ref_agents, id);

  if (REF_EMPTY != ref_agent_previous(ref_agents, id))
    ref_agent_next(ref_agents, ref_agent_previous(ref_agents, id)) =
        ref_agent_next(ref_agents, id);

  if (REF_EMPTY != ref_agent_next(ref_agents, id))
    ref_agent_previous(ref_agents, ref_agent_next(ref_agents, id)) =
        ref_agent_previous(ref_agents, id);

  ref_agent_mode(ref_agents, id) = REF_AGENT_UNUSED;
  ref_agent_previous(ref_agents, id) = REF_EMPTY;
  ref_agent_next(ref_agents, id) = ref_agents->blank;
  ref_agents->blank = id;

  ref_agents_n(ref_agents)--;

  return REF_SUCCESS;
}

REF_STATUS ref_agents_pop(REF_AGENTS ref_agents, REF_INT *node, REF_INT *part,
                          REF_INT *seed, REF_DBL *xyz) {
  REF_INT i, id;

  if (REF_EMPTY == ref_agents->last) return REF_NOT_FOUND;

  id = ref_agents->last;

  *node = ref_agent_node(ref_agents, id);
  *part = ref_agent_part(ref_agents, id);
  *seed = ref_agent_seed(ref_agents, id);
  for (i = 0; i < 3; i++) xyz[i] = ref_agent_xyz(ref_agents, i, id);

  RSS(ref_agents_remove(ref_agents, id), "rm");

  return REF_SUCCESS;
}

REF_STATUS ref_agents_delete(REF_AGENTS ref_agents, REF_INT node) {
  REF_INT id;

  if (node < 0) return REF_INVALID;

  each_active_ref_agent(ref_agents, id) {
    if (node == ref_agent_node(ref_agents, id))
      RSS(ref_agents_remove(ref_agents, id), "rm");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_agents_dest(REF_AGENTS ref_agents, REF_INT id, REF_INT *dest) {
  *dest = REF_EMPTY;
  switch (ref_agent_mode(ref_agents, id)) {
    case REF_AGENT_WALKING:
      *dest = ref_agent_part(ref_agents, id);
      break;
    case REF_AGENT_AT_BOUNDARY:
      *dest = ref_agent_home(ref_agents, id);
      break;
    case REF_AGENT_HOP_PART:
      *dest = ref_agent_part(ref_agents, id);
      break;
    case REF_AGENT_SUGGESTION:
      *dest = ref_agent_home(ref_agents, id);
      break;
    case REF_AGENT_ENCLOSING:
      *dest = ref_agent_home(ref_agents, id);
      break;
    case REF_AGENT_TERMINATED:
      *dest = ref_agent_home(ref_agents, id);
      break;
    default:
      THROW("has no dest");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_agents_migrate(REF_AGENTS ref_agents) {
  REF_MPI ref_mpi = ref_agents->ref_mpi;
  REF_INT i, id, nsend, nrecv, dest, rec;
  REF_INT n_ints, n_dbls;
  REF_INT *destination, *send_int, *recv_int;
  REF_DBL *send_dbl, *recv_dbl;
  nsend = 0;
  each_active_ref_agent(ref_agents, id) {
    RSS(ref_agents_dest(ref_agents, id, &dest), "dest");
    if (ref_mpi_rank(ref_mpi) != dest) {
      nsend++;
    }
  }
  n_ints = 7;
  n_dbls = 7;
  ref_malloc_init(destination, nsend, REF_INT, REF_EMPTY);
  ref_malloc_init(send_int, nsend * n_ints, REF_INT, REF_EMPTY);
  ref_malloc_init(send_dbl, nsend * n_dbls, REF_DBL, 0.0);
  nsend = 0;
  each_active_ref_agent(ref_agents, id) {
    RSS(ref_agents_dest(ref_agents, id, &dest), "dest");
    if (ref_mpi_rank(ref_mpi) != dest) {
      destination[nsend] = dest;
      send_int[0 + nsend * n_ints] = (REF_INT)ref_agent_mode(ref_agents, id);
      send_int[1 + nsend * n_ints] = ref_agent_home(ref_agents, id);
      send_int[2 + nsend * n_ints] = ref_agent_node(ref_agents, id);
      send_int[3 + nsend * n_ints] = ref_agent_part(ref_agents, id);
      send_int[4 + nsend * n_ints] = ref_agent_seed(ref_agents, id);
      send_int[5 + nsend * n_ints] = ref_agent_global(ref_agents, id);
      send_int[6 + nsend * n_ints] = ref_agent_step(ref_agents, id);

      for (i = 0; i < 3; i++)
        send_dbl[i + nsend * n_dbls] = ref_agent_xyz(ref_agents, i, id);
      for (i = 0; i < 4; i++)
        send_dbl[3 + i + nsend * n_dbls] = ref_agent_bary(ref_agents, i, id);

      RSS(ref_agents_remove(ref_agents, id), "poof");

      nsend++;
    }
  }

  RSS(ref_mpi_blindsend(ref_mpi, destination, (void *)send_int, n_ints, nsend,
                        (void **)(&recv_int), &nrecv, REF_INT_TYPE),
      "is");
  RSS(ref_mpi_blindsend(ref_mpi, destination, (void *)send_dbl, n_dbls, nsend,
                        (void **)(&recv_dbl), &nrecv, REF_DBL_TYPE),
      "ds");
  ref_free(send_dbl);
  ref_free(send_int);
  ref_free(destination);

  for (rec = 0; rec < nrecv; rec++) {
    RSS(ref_agents_new(ref_agents, &id), "new");
    ref_agent_mode(ref_agents, id) = (REF_AGENT_MODE)recv_int[0 + rec * n_ints];
    ref_agent_home(ref_agents, id) = recv_int[1 + rec * n_ints];
    ref_agent_node(ref_agents, id) = recv_int[2 + rec * n_ints];
    ref_agent_part(ref_agents, id) = recv_int[3 + rec * n_ints];
    ref_agent_seed(ref_agents, id) = recv_int[4 + rec * n_ints];
    ref_agent_global(ref_agents, id) = recv_int[5 + rec * n_ints];
    ref_agent_step(ref_agents, id) = recv_int[6 + rec * n_ints];

    for (i = 0; i < 3; i++)
      ref_agent_xyz(ref_agents, i, id) = recv_dbl[i + rec * n_dbls];
    for (i = 0; i < 4; i++)
      ref_agent_bary(ref_agents, i, id) = recv_dbl[3 + i + rec * n_dbls];
  }

  ref_free(recv_dbl);
  ref_free(recv_int);

  return REF_SUCCESS;
}
