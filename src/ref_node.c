
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

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

/* REF_EMPTY is terminatior, next avalable is shifted by 2*/
#define next2index(next) (-(next)-2)
#define index2next(index) (-2 - (index))

REF_STATUS ref_node_create(REF_NODE *ref_node_ptr, REF_MPI ref_mpi) {
  REF_INT max, node;
  REF_NODE ref_node;

  ref_malloc(*ref_node_ptr, 1, REF_NODE_STRUCT);

  ref_node = *ref_node_ptr;

  max = 20;

  ref_node_n(ref_node) = 0;
  ref_node_max(ref_node) = max;

  ref_malloc(ref_node->global, max, REF_INT);

  for (node = 0; node < ref_node_max(ref_node); node++)
    ref_node->global[node] = index2next(node + 1);
  ref_node->global[(ref_node->max) - 1] = REF_EMPTY;
  ref_node->blank = index2next(0);

  ref_malloc(ref_node->sorted_global, max, REF_INT);
  ref_malloc(ref_node->sorted_local, max, REF_INT);

  ref_malloc(ref_node->part, max, REF_INT);
  ref_malloc(ref_node->age, max, REF_INT);

  ref_malloc(ref_node->real, REF_NODE_REAL_PER * max, REF_DBL);

  ref_node_naux(ref_node) = 0;
  ref_node->aux = NULL;

  ref_node_mpi(ref_node) = ref_mpi; /* reference only */

  ref_node_n_unused(ref_node) = 0;
  ref_node_max_unused(ref_node) = 10;
  ref_malloc(ref_node->unused_global, ref_node_max_unused(ref_node), REF_INT);

  ref_node->old_n_global = REF_EMPTY;
  ref_node->new_n_global = REF_EMPTY;

  ref_node->twod_mid_plane = 0.5;
  ref_node->min_volume = 1.0e-15;
  ref_node->min_uv_area = 1.0e-12;

  ref_node->tet_quality = REF_NODE_JAC_QUALITY;
  ref_node->tri_quality = REF_NODE_JAC_QUALITY;

  return REF_SUCCESS;
}

REF_STATUS ref_node_free(REF_NODE ref_node) {
  if (NULL == (void *)ref_node) return REF_NULL;
  ref_free(ref_node->unused_global);
  /* ref_mpi reference only */
  ref_free(ref_node->aux);
  ref_free(ref_node->real);
  ref_free(ref_node->age);
  ref_free(ref_node->part);
  ref_free(ref_node->sorted_local);
  ref_free(ref_node->sorted_global);
  ref_free(ref_node->global);
  ref_free(ref_node);
  return REF_SUCCESS;
}

REF_STATUS ref_node_deep_copy(REF_NODE *ref_node_ptr, REF_NODE original) {
  REF_INT max, node, i;
  REF_NODE ref_node;

  ref_malloc(*ref_node_ptr, 1, REF_NODE_STRUCT);
  ref_node = *ref_node_ptr;

  max = ref_node_max(original);

  ref_node_n(ref_node) = ref_node_n(original);
  ref_node_max(ref_node) = max;

  ref_malloc(ref_node->global, max, REF_INT);
  ref_node->blank = original->blank;
  for (node = 0; node < max; node++)
    ref_node->global[node] = original->global[node];

  ref_malloc(ref_node->sorted_global, max, REF_INT);
  ref_malloc(ref_node->sorted_local, max, REF_INT);
  for (node = 0; node < max; node++)
    ref_node->sorted_global[node] = original->sorted_global[node];
  for (node = 0; node < max; node++)
    ref_node->sorted_local[node] = original->sorted_local[node];

  ref_malloc(ref_node->part, max, REF_INT);
  for (node = 0; node < max; node++)
    ref_node_part(ref_node, node) = ref_node_part(original, node);

  ref_malloc(ref_node->age, max, REF_INT);
  for (node = 0; node < max; node++)
    ref_node_age(ref_node, node) = ref_node_age(original, node);

  ref_malloc(ref_node->real, REF_NODE_REAL_PER * max, REF_DBL);
  for (node = 0; node < max; node++)
    for (i = 0; i < REF_NODE_REAL_PER; i++)
      ref_node_real(ref_node, i, node) = ref_node_real(original, i, node);

  ref_node_naux(ref_node) = ref_node_naux(original);
  ref_node->aux = NULL;
  if (ref_node_naux(original) > 0) {
    ref_realloc(ref_node->aux, ref_node_naux(ref_node) * ref_node_max(ref_node),
                REF_DBL);
    for (node = 0; node < max; node++)
      for (i = 0; i < ref_node_naux(ref_node); i++)
        ref_node_aux(ref_node, i, node) = ref_node_aux(original, i, node);
  }

  ref_node_mpi(ref_node) = ref_node_mpi(original); /* reference only */

  ref_node->n_unused = original->n_unused;
  ref_node->max_unused = original->max_unused;
  ref_malloc(ref_node->unused_global, ref_node_max_unused(ref_node), REF_INT);
  for (i = 0; i < ref_node_n_unused(ref_node); i++)
    ref_node->unused_global[i] = original->unused_global[i];

  ref_node->old_n_global = original->old_n_global;
  ref_node->new_n_global = original->new_n_global;

  ref_node->twod_mid_plane = original->twod_mid_plane;
  ref_node->min_volume = original->min_volume;
  ref_node->min_uv_area = original->min_uv_area;

  ref_node->tet_quality = original->tet_quality;
  ref_node->tri_quality = original->tri_quality;

  return REF_SUCCESS;
}

REF_STATUS ref_node_pack(REF_NODE ref_node, REF_INT *o2n, REF_INT *n2o) {
  REF_INT i, node;
  REF_NODE copy;
  RSS(ref_node_deep_copy(&copy, ref_node), "make a copy first");

  for (node = 0; node < ref_node_n(ref_node); node++)
    ref_node->global[node] = copy->global[n2o[node]];
  if (ref_node_n(ref_node) < ref_node_max(ref_node)) {
    for (node = ref_node_n(ref_node); node < ref_node_max(ref_node); node++)
      ref_node->global[node] = index2next(node + 1);
    ref_node->global[(ref_node->max) - 1] = REF_EMPTY;
    ref_node->blank = index2next(ref_node_n(ref_node));
  } else {
    ref_node->blank = REF_EMPTY;
  }

  for (node = 0; node < ref_node_n(ref_node); node++)
    ref_node->sorted_local[node] = o2n[copy->sorted_local[node]];

  for (node = 0; node < ref_node_n(ref_node); node++)
    ref_node->part[node] = copy->part[n2o[node]];

  for (node = 0; node < ref_node_n(ref_node); node++)
    ref_node->age[node] = copy->age[n2o[node]];

  for (node = 0; node < ref_node_n(ref_node); node++)
    for (i = 0; i < REF_NODE_REAL_PER; i++)
      ref_node_real(ref_node, i, node) = ref_node_real(copy, i, n2o[node]);

  if (ref_node_naux(ref_node) > 0) {
    for (node = 0; node < ref_node_n(ref_node); node++)
      for (i = 0; i < ref_node_naux(ref_node); i++)
        ref_node_aux(ref_node, i, node) = ref_node_aux(copy, i, n2o[node]);
  }

  RSS(ref_node_free(copy), "release copy");

  return REF_SUCCESS;
}

REF_STATUS ref_node_inspect(REF_NODE ref_node) {
  REF_INT node;
  printf("ref_node = %p\n", (void *)ref_node);
  printf(" n = %d\n", ref_node_n(ref_node));
  printf(" max = %d\n", ref_node_max(ref_node));
  printf(" blank = %d\n", ref_node->blank);
  for (node = 0; node < ref_node_max(ref_node); node++)
    if (0 <= ref_node->global[node])
      printf(" global[%d] = %3d; part[%d] = %3d;\n", node,
             ref_node->global[node], node, ref_node->part[node]);
  for (node = 0; node < ref_node_n(ref_node); node++)
    printf(" sorted_global[%d] = %d sorted_local[%d] = %d\n", node,
           ref_node->sorted_global[node], node, ref_node->sorted_local[node]);
  printf(" old_n_global = %d\n", ref_node->old_n_global);
  printf(" new_n_global = %d\n", ref_node->new_n_global);
  return REF_SUCCESS;
}

REF_STATUS ref_node_location(REF_NODE ref_node, REF_INT node) {
  printf("ref_node %d", node);
  if (ref_node_valid(ref_node, node)) {
    printf(" part %d mine %d (%.15e,%.15e,%.15e)\n",
           ref_node_part(ref_node, node), ref_mpi_rank(ref_node_mpi(ref_node)),
           ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
           ref_node_xyz(ref_node, 2, node));

  } else {
    printf(" invalid\n");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tattle_global(REF_NODE ref_node, REF_INT global) {
  REF_INT local_from_sorted;
  REF_INT local_from_exhaustive;
  REF_BOOL found_from_sorted;
  REF_BOOL found_from_exhaustive;
  REF_INT node;

  REF_STATUS ref_status;

  ref_status = ref_node_local(ref_node, global, &local_from_sorted);
  if (REF_NOT_FOUND == ref_status) {
    found_from_sorted = REF_FALSE;
  } else {
    RSS(ref_status, "local search");
    found_from_sorted = REF_TRUE;
  }

  found_from_exhaustive = REF_FALSE;
  local_from_exhaustive = REF_EMPTY;
  each_ref_node_valid_node(ref_node, node) {
    if (global == ref_node_global(ref_node, node)) {
      if (found_from_exhaustive) RSS(REF_FAILURE, "twice");
      local_from_exhaustive = node;
      found_from_exhaustive = REF_TRUE;
    }
  }

  printf("%d: global %d: search%d %d exhast%d %d\n",
         ref_mpi_rank(ref_node_mpi(ref_node)), global, found_from_sorted,
         local_from_sorted, found_from_exhaustive, local_from_exhaustive);

  return REF_SUCCESS;
}

static REF_STATUS ref_node_add_core(REF_NODE ref_node, REF_INT global,
                                    REF_INT *node) {
  REF_INT orig, chunk, extra;

  if (global < 0) RSS(REF_INVALID, "invalid global node");

  if (REF_EMPTY == ref_node->blank) {
    orig = ref_node_max(ref_node);
    chunk = MAX(5000, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_node->max = orig + chunk;
    ref_realloc(ref_node->global, ref_node_max(ref_node), REF_INT);
    for (extra = orig; extra < ref_node_max(ref_node); extra++)
      ref_node->global[extra] = index2next(extra + 1);
    ref_node->global[ref_node_max(ref_node) - 1] = REF_EMPTY;
    ref_node->blank = index2next(orig);

    ref_realloc(ref_node->sorted_global, ref_node_max(ref_node), REF_INT);
    ref_realloc(ref_node->sorted_local, ref_node_max(ref_node), REF_INT);

    ref_realloc(ref_node->part, ref_node_max(ref_node), REF_INT);
    ref_realloc(ref_node->age, ref_node_max(ref_node), REF_INT);

    ref_realloc(ref_node->real,
                ((unsigned long)REF_NODE_REAL_PER *
                 (unsigned long)ref_node_max(ref_node)),
                REF_DBL);

    if (ref_node_naux(ref_node) > 0)
      ref_realloc(ref_node->aux,
                  ((unsigned long)ref_node_naux(ref_node) *
                   (unsigned long)ref_node_max(ref_node)),
                  REF_DBL);
  }

  *node = next2index(ref_node->blank);
  ref_node->blank = ref_node->global[*node];

  ref_node->global[*node] = global;
  ref_node->part[*node] =
      ref_mpi_rank(ref_node_mpi(ref_node)); /*local default*/
  ref_node->age[*node] = 0;                 /* default new born */

  RSS(ref_node_metric_form(ref_node, *node, 1, 0, 0, 1, 0, 1), "set ident");

  (ref_node->n)++;
  return REF_SUCCESS;
}

REF_STATUS ref_node_add(REF_NODE ref_node, REF_INT global, REF_INT *node) {
  REF_INT location, insert_point;
  REF_STATUS status;

  if (global < 0) RSS(REF_INVALID, "invalid global node");

  status = ref_node_local(ref_node, global, node);
  if (REF_SUCCESS == status) return REF_SUCCESS;

  RSS(ref_node_add_core(ref_node, global, node), "core");

  /* general case of non-ascending global node, requires:
     search and shift (but looks to see if bigger than last early) */
  insert_point = 0;
  for (location = ref_node_n(ref_node) - 2; location >= 0; location--) {
    if (ref_node->sorted_global[location] < global) {
      insert_point = location + 1;
      break;
    }
  }

  /* shift down to clear insert_point */
  for (location = ref_node_n(ref_node) - 1; location > insert_point; location--)
    ref_node->sorted_global[location] = ref_node->sorted_global[location - 1];
  for (location = ref_node_n(ref_node) - 1; location > insert_point; location--)
    ref_node->sorted_local[location] = ref_node->sorted_local[location - 1];

  /* insert in empty location */
  ref_node->sorted_global[insert_point] = global;
  ref_node->sorted_local[insert_point] = *node;

  return REF_SUCCESS;
}

REF_STATUS ref_node_add_many(REF_NODE ref_node, REF_INT n,
                             REF_INT *global_orig) {
  REF_STATUS status;
  REF_INT i, j, local, new;

  REF_INT *global;
  REF_INT *sorted;

  /* copy, removing existing nodes from list */

  ref_malloc(global, n, REF_INT);

  new = 0;
  for (i = 0; i < n; i++) {
    status = ref_node_local(ref_node, global_orig[i], &local);
    if (REF_NOT_FOUND == status) {
      global[new] = global_orig[i];
      new ++;
    }
  }

  /* remove duplicates from list so core add can be used with existing check */

  ref_malloc(sorted, new, REF_INT);

  RSS(ref_sort_heap_int(new, global, sorted), "heap");

  j = 0;
  for (i = 1; i < new; i++) {
    if (global[sorted[i]] != global[sorted[j]]) {
      j = i;
      continue;
    }
    global[sorted[i]] = REF_EMPTY;
  }

  /* add remaining via core */

  for (i = 0; i < new; i++)
    if (REF_EMPTY != global[i]) {
      RSS(ref_node_add_core(ref_node, global[i], &local), "add core");
    }

  RSS(ref_node_rebuild_sorted_global(ref_node), "rebuild globals");

  ref_free(sorted);
  ref_free(global);

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove(REF_NODE ref_node, REF_INT node) {
  REF_INT location, sorted_node;
  if (!ref_node_valid(ref_node, node)) return REF_INVALID;

  RSS(ref_sort_search(ref_node_n(ref_node), ref_node->sorted_global,
                      ref_node->global[node], &location),
      "find global in sort list");

  for (sorted_node = location; sorted_node < ref_node_n(ref_node) - 1;
       sorted_node++)
    ref_node->sorted_global[sorted_node] =
        ref_node->sorted_global[sorted_node + 1];
  for (sorted_node = location; sorted_node < ref_node_n(ref_node) - 1;
       sorted_node++)
    ref_node->sorted_local[sorted_node] =
        ref_node->sorted_local[sorted_node + 1];

  RSS(ref_node_push_unused(ref_node, ref_node->global[node]),
      "store unused global");

  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove_invalidates_sorted(REF_NODE ref_node, REF_INT node) {
  if (!ref_node_valid(ref_node, node)) return REF_INVALID;

  RSS(ref_node_push_unused(ref_node, ref_node->global[node]),
      "store unused global");

  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove_without_global(REF_NODE ref_node, REF_INT node) {
  REF_INT location, sorted_node;
  if (!ref_node_valid(ref_node, node)) return REF_INVALID;

  RSS(ref_sort_search(ref_node_n(ref_node), ref_node->sorted_global,
                      ref_node->global[node], &location),
      "find global in sort list");

  for (sorted_node = location; sorted_node < ref_node_n(ref_node) - 1;
       sorted_node++)
    ref_node->sorted_global[sorted_node] =
        ref_node->sorted_global[sorted_node + 1];
  for (sorted_node = location; sorted_node < ref_node_n(ref_node) - 1;
       sorted_node++)
    ref_node->sorted_local[sorted_node] =
        ref_node->sorted_local[sorted_node + 1];

  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove_without_global_invalidates_sorted(REF_NODE ref_node,
                                                             REF_INT node) {
  if (!ref_node_valid(ref_node, node)) return REF_INVALID;

  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_rebuild_sorted_global(REF_NODE ref_node) {
  REF_INT node, nnode, *pack;

  ref_malloc(pack, ref_node_n(ref_node), REF_INT);

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    ref_node->sorted_global[nnode] = ref_node->global[node];
    pack[nnode] = node;
    nnode++;
  }

  RSS(ref_sort_heap_int(ref_node_n(ref_node), ref_node->sorted_global,
                        ref_node->sorted_local),
      "heap");

  for (node = 0; node < ref_node_n(ref_node); node++) {
    ref_node->sorted_local[node] = pack[ref_node->sorted_local[node]];
    ref_node->sorted_global[node] =
        ref_node->global[ref_node->sorted_local[node]];
  }

  ref_free(pack);
  return REF_SUCCESS;
}

REF_STATUS ref_node_initialize_n_global(REF_NODE ref_node, REF_INT n_global) {
  ref_node->old_n_global = n_global;
  ref_node->new_n_global = n_global;

  return REF_SUCCESS;
}

REF_STATUS ref_node_next_global(REF_NODE ref_node, REF_INT *global) {
  if (0 < ref_node_n_unused(ref_node)) {
    RSS(ref_node_pop_unused(ref_node, global),
        "grab an unused global from list");
  } else {
    if (REF_EMPTY == ref_node->new_n_global)
      RSS(ref_node_initialize_n_global(ref_node, ref_node_n(ref_node)),
          "init with n");
    (*global) = ref_node->new_n_global;
    (ref_node->new_n_global)++;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_node_synchronize_globals(REF_NODE ref_node) {
  RSS(ref_node_shift_new_globals(ref_node), "shift");
  RSS(ref_node_eliminate_unused_globals(ref_node), "shift");

  return REF_SUCCESS;
}

REF_STATUS ref_node_shift_new_globals(REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT new_nodes;
  REF_INT *everyones_new_nodes;
  REF_INT offset, proc, total_new_nodes, node;

  ref_malloc(everyones_new_nodes, ref_mpi_n(ref_mpi), REF_INT);

  new_nodes = ref_node->new_n_global - ref_node->old_n_global;

  RSS(ref_mpi_allgather(ref_mpi, &new_nodes, everyones_new_nodes, REF_INT_TYPE),
      "allgather");

  offset = 0;
  for (proc = 0; proc < ref_mpi_rank(ref_mpi); proc++)
    offset += everyones_new_nodes[proc];

  total_new_nodes = 0;
  each_ref_mpi_part(ref_mpi, proc) total_new_nodes += everyones_new_nodes[proc];

  ref_free(everyones_new_nodes);

  if (0 != offset) {
    each_ref_node_valid_node(ref_node, node) {
      if (ref_node_global(ref_node, node) >= ref_node->old_n_global) {
        (ref_node->global[node]) += offset;
      }
    }
    for (node = ref_node_n(ref_node) - 1;
         node >= 0 && ref_node->sorted_global[node] >= ref_node->old_n_global;
         node--)
      ref_node->sorted_global[node] += offset;

    RSS(ref_node_shift_unused(ref_node, ref_node->old_n_global, offset),
        "shift");
  }

  RSS(ref_node_initialize_n_global(ref_node,
                                   total_new_nodes + ref_node->old_n_global),
      "re-init");

  return REF_SUCCESS;
}

REF_STATUS ref_node_implicit_global_from_local(REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT node, nnode, *global;
  REF_INT *everyones_nnode, offset, proc;

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      nnode++;
    }
  }

  ref_malloc(everyones_nnode, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &nnode, everyones_nnode, REF_INT_TYPE),
      "allgather");

  offset = 0;
  for (proc = 0; proc < ref_mpi_rank(ref_mpi); proc++)
    offset += everyones_nnode[proc];

  ref_free(everyones_nnode);
  ref_malloc(global, ref_node_max(ref_node), REF_INT);

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      global[node] = offset + nnode;
      nnode++;
    }
  }

  RSS(ref_node_ghost_int(ref_node, global, 1), "ghost int");

  each_ref_node_valid_node(ref_node, node) {
    ref_node->global[node] = global[node];
  }

  ref_free(global);

  RSS(ref_node_rebuild_sorted_global(ref_node), "rebuild globals");

  return REF_SUCCESS;
}

REF_STATUS ref_node_eliminate_unused_globals(REF_NODE ref_node) {
  REF_INT sort, offset, local;

  RSS(ref_node_allgather_unused(ref_node), "gather unused global");
  RSS(ref_node_sort_unused(ref_node), "sort unused global");

  offset = 0;
  for (sort = 0; sort < ref_node_n(ref_node); sort++) {
    while ((offset < ref_node_n_unused(ref_node)) &&
           (ref_node->unused_global[offset] < ref_node->sorted_global[sort])) {
      offset++;
    }
    local = ref_node->sorted_local[sort];
    ref_node->global[local] -= offset; /* move to separate loop for cashe? */
    ref_node->sorted_global[sort] -= offset;
  }

  RSS(ref_node_initialize_n_global(
          ref_node, ref_node->old_n_global - ref_node_n_unused(ref_node)),
      "re-init");

  RSS(ref_node_erase_unused(ref_node), "erase unused list");

  return REF_SUCCESS;
}

REF_STATUS ref_node_collect_ghost_age(REF_NODE ref_node) {
  RSS(ref_node_localize_ghost_int(ref_node, (ref_node->age)),
      "localize ghost age");
  return REF_SUCCESS;
}

REF_STATUS ref_node_local(REF_NODE ref_node, REF_GLOB global, REF_INT *local) {
  REF_INT location;

  (*local) = REF_EMPTY;

  RAISE(ref_sort_search(ref_node_n(ref_node), ref_node->sorted_global, global,
                        &location));

  if ((location) == REF_EMPTY) return REF_NOT_FOUND;

  (*local) = ref_node->sorted_local[location];

  return REF_SUCCESS;
}

REF_STATUS ref_node_compact(REF_NODE ref_node, REF_INT **o2n_ptr,
                            REF_INT **n2o_ptr) {
  REF_INT node;
  REF_INT nnode;
  REF_INT *o2n, *n2o;

  ref_malloc_init(*o2n_ptr, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  o2n = *o2n_ptr;
  ref_malloc(*n2o_ptr, ref_node_n(ref_node), REF_INT);
  n2o = *n2o_ptr;

  nnode = 0;

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_node_mpi(ref_node)) == ref_node_part(ref_node, node)) {
      o2n[node] = nnode;
      nnode++;
    }
  }

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_node_mpi(ref_node)) != ref_node_part(ref_node, node)) {
      o2n[node] = nnode;
      nnode++;
    }
  }

  RES(nnode, ref_node_n(ref_node), "nnode miscount");

  each_ref_node_valid_node(ref_node, node) n2o[o2n[node]] = node;

  return REF_SUCCESS;
}

REF_STATUS ref_node_ghost_real(REF_NODE ref_node) {
  RSS(ref_node_ghost_dbl(ref_node, ref_node->real, REF_NODE_REAL_PER),
      "ghost dbl");
  return REF_SUCCESS;
}

REF_STATUS ref_node_ghost_int(REF_NODE ref_node, REF_INT *vector,
                              REF_INT ldim) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_INT *a_vector, *b_vector;
  REF_INT i, local;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      a_size[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) { a_total += a_size[part]; }
  ref_malloc(a_global, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) { b_total += b_size[part]; }
  ref_malloc(b_global, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      a_next[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_size, b_global, b_size, 1,
                        REF_INT_TYPE),
      "alltoallv global");

  if (a_total < REF_INT_MAX / ldim && b_total < REF_INT_MAX / ldim) {
    ref_malloc(a_vector, ldim * a_total, REF_INT);
    ref_malloc(b_vector, ldim * b_total, REF_INT);
    for (node = 0; node < b_total; node++) {
      RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        b_vector[i + ldim * node] = vector[i + ldim * local];
    }

    RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, ldim,
                          REF_INT_TYPE),
        "alltoallv global");

    for (node = 0; node < a_total; node++) {
      RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        vector[i + ldim * local] = a_vector[i + ldim * node];
    }
    free(b_vector);
    free(a_vector);
  } else {
    ref_malloc(a_vector, a_total, REF_INT);
    ref_malloc(b_vector, b_total, REF_INT);
    for (i = 0; i < ldim; i++) {
      for (node = 0; node < b_total; node++) {
        RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
        b_vector[node] = vector[i + ldim * local];
      }

      RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, 1,
                            REF_INT_TYPE),
          "alltoallv global");

      for (node = 0; node < a_total; node++) {
        RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
        vector[i + ldim * local] = a_vector[node];
      }
    }
    free(b_vector);
    free(a_vector);
  }

  free(a_next);
  free(b_global);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_node_ghost_glob(REF_NODE ref_node, REF_GLOB *vector,
                               REF_INT ldim) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_GLOB *a_vector, *b_vector;
  REF_INT i, local;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      a_size[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) { a_total += a_size[part]; }
  ref_malloc(a_global, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) { b_total += b_size[part]; }
  ref_malloc(b_global, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      a_next[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_size, b_global, b_size, 1,
                        REF_INT_TYPE),
      "alltoallv global");

  if (a_total < REF_INT_MAX / ldim && b_total < REF_INT_MAX / ldim) {
    ref_malloc(a_vector, ldim * a_total, REF_GLOB);
    ref_malloc(b_vector, ldim * b_total, REF_GLOB);
    for (node = 0; node < b_total; node++) {
      RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        b_vector[i + ldim * node] = vector[i + ldim * local];
    }

    RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, ldim,
                          REF_GLOB_TYPE),
        "alltoallv global");

    for (node = 0; node < a_total; node++) {
      RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        vector[i + ldim * local] = a_vector[i + ldim * node];
    }
    free(b_vector);
    free(a_vector);
  } else {
    ref_malloc(a_vector, a_total, REF_GLOB);
    ref_malloc(b_vector, b_total, REF_GLOB);
    for (i = 0; i < ldim; i++) {
      for (node = 0; node < b_total; node++) {
        RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
        b_vector[node] = vector[i + ldim * local];
      }

      RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, 1,
                            REF_GLOB_TYPE),
          "alltoallv global");

      for (node = 0; node < a_total; node++) {
        RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
        vector[i + ldim * local] = a_vector[node];
      }
    }
    free(b_vector);
    free(a_vector);
  }

  free(a_next);
  free(b_global);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_node_ghost_dbl(REF_NODE ref_node, REF_DBL *vector,
                              REF_INT ldim) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_vector, *b_vector;
  REF_INT i, local;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      a_size[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_global, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_global, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      a_next[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_size, b_global, b_size, 1,
                        REF_INT_TYPE),
      "alltoallv global");

  if (a_total < REF_INT_MAX / ldim && b_total < REF_INT_MAX / ldim) {
    ref_malloc(a_vector, ldim * a_total, REF_DBL);
    ref_malloc(b_vector, ldim * b_total, REF_DBL);
    for (node = 0; node < b_total; node++) {
      RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        b_vector[i + ldim * node] = vector[i + ldim * local];
    }

    RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, ldim,
                          REF_DBL_TYPE),
        "alltoallv global");

    for (node = 0; node < a_total; node++) {
      RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
      for (i = 0; i < ldim; i++)
        vector[i + ldim * local] = a_vector[i + ldim * node];
    }
    free(b_vector);
    free(a_vector);
  } else {
    ref_malloc(a_vector, a_total, REF_DBL);
    ref_malloc(b_vector, b_total, REF_DBL);
    for (i = 0; i < ldim; i++) {
      for (node = 0; node < b_total; node++) {
        RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
        b_vector[node] = vector[i + ldim * local];
      }

      RSS(ref_mpi_alltoallv(ref_mpi, b_vector, b_size, a_vector, a_size, 1,
                            REF_DBL_TYPE),
          "alltoallv global");

      for (node = 0; node < a_total; node++) {
        RSS(ref_node_local(ref_node, a_global[node], &local), "g2l");
        vector[i + ldim * local] = a_vector[node];
      }
    }
    free(b_vector);
    free(a_vector);
  }
  free(a_next);
  free(b_global);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_node_localize_ghost_int(REF_NODE ref_node, REF_INT *scalar) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_GLOB *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_INT *a_scalar, *b_scalar;
  REF_INT local;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      a_size[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) { a_total += a_size[part]; }
  ref_malloc(a_global, a_total, REF_GLOB);
  ref_malloc(a_scalar, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) { b_total += b_size[part]; }
  ref_malloc(b_global, b_total, REF_GLOB);
  ref_malloc(b_scalar, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      a_scalar[a_next[part]] = scalar[node];
      a_next[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_size, b_global, b_size, 1,
                        REF_GLOB_TYPE),
      "alltoallv global");

  RSS(ref_mpi_alltoallv(ref_mpi, a_scalar, a_size, b_scalar, b_size, 1,
                        REF_INT_TYPE),
      "alltoallv scalar");

  for (node = 0; node < b_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
    scalar[local] += b_scalar[node];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node)) {
      scalar[node] = 0;
    }
  }

  free(a_next);
  free(b_scalar);
  free(b_global);
  free(a_scalar);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_node_edge_twod(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                              REF_BOOL *twod) {
  REF_DBL mid_plane = ref_node_twod_mid_plane(ref_node);

  *twod = ((ref_node_xyz(ref_node, 1, node0) < mid_plane) &&
           (ref_node_xyz(ref_node, 1, node1) < mid_plane));

  return REF_SUCCESS;
}

REF_STATUS ref_node_node_twod(REF_NODE ref_node, REF_INT node, REF_BOOL *twod) {
  REF_DBL mid_plane = ref_node_twod_mid_plane(ref_node);

  *twod = (ref_node_xyz(ref_node, 1, node) < mid_plane);

  return REF_SUCCESS;
}

REF_STATUS ref_node_metric_form(REF_NODE ref_node, REF_INT node, REF_DBL m11,
                                REF_DBL m12, REF_DBL m13, REF_DBL m22,
                                REF_DBL m23, REF_DBL m33) {
  REF_DBL m[6];
  m[0] = m11;
  m[1] = m12;
  m[2] = m13;
  m[3] = m22;
  m[4] = m23;
  m[5] = m33;
  RSS(ref_node_metric_set(ref_node, node, m), "set");
  return REF_SUCCESS;
}

REF_STATUS ref_node_metric_set(REF_NODE ref_node, REF_INT node, REF_DBL *m) {
  REF_INT i;
  REF_DBL log_m[6];
  for (i = 0; i < 6; i++) {
    ((ref_node)->real[(i + 3) + REF_NODE_REAL_PER * (node)]) = m[i];
  }
  RSS(ref_matrix_log_m(m, log_m), "exp");
  for (i = 0; i < 6; i++) {
    ((ref_node)->real[(i + 9) + REF_NODE_REAL_PER * (node)]) = log_m[i];
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_metric_get(REF_NODE ref_node, REF_INT node, REF_DBL *m) {
  REF_INT i;
  for (i = 0; i < 6; i++) {
    m[i] = ((ref_node)->real[(i + 3) + REF_NODE_REAL_PER * (node)]);
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_metric_set_log(REF_NODE ref_node, REF_INT node,
                                   REF_DBL *log_m) {
  REF_INT i;
  REF_DBL m[6];
  for (i = 0; i < 6; i++) {
    ((ref_node)->real[(i + 9) + REF_NODE_REAL_PER * (node)]) = log_m[i];
  }
  RSS(ref_matrix_exp_m(log_m, m), "exp");
  for (i = 0; i < 6; i++) {
    ((ref_node)->real[(i + 3) + REF_NODE_REAL_PER * (node)]) = m[i];
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_metric_get_log(REF_NODE ref_node, REF_INT node,
                                   REF_DBL *log_m) {
  REF_INT i;
  for (i = 0; i < 6; i++) {
    log_m[i] = ((ref_node)->real[(i + 9) + REF_NODE_REAL_PER * (node)]);
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_ratio(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                          REF_DBL *ratio) {
  REF_DBL direction[3], length;
  REF_DBL ratio0, ratio1;
  REF_DBL r, r_min, r_max;
  REF_DBL m[6];

  if (!ref_node_valid(ref_node, node0) || !ref_node_valid(ref_node, node1))
    RSS(REF_INVALID, "node invalid");

  direction[0] =
      (ref_node_xyz(ref_node, 0, node1) - ref_node_xyz(ref_node, 0, node0));
  direction[1] =
      (ref_node_xyz(ref_node, 1, node1) - ref_node_xyz(ref_node, 1, node0));
  direction[2] =
      (ref_node_xyz(ref_node, 2, node1) - ref_node_xyz(ref_node, 2, node0));

  length = ref_math_dot(direction, direction);
  length = sqrt(length);

  if (!ref_math_divisible(direction[0], length) ||
      !ref_math_divisible(direction[1], length) ||
      !ref_math_divisible(direction[2], length)) {
    *ratio = 0.0;
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get(ref_node, node0, m), "node0 m");
  ratio0 = ref_matrix_sqrt_vt_m_v(m, direction);
  RSS(ref_node_metric_get(ref_node, node1, m), "node1 m");
  ratio1 = ref_matrix_sqrt_vt_m_v(m, direction);

  /* Loseille Lohner IMR 18 (2009) pg 613 */
  /* Alauzet Finite Elements in Analysis and Design 46 (2010) pg 185 */

  if (ratio0 < 1.0e-12 || ratio1 < 1.0e-12) {
    *ratio = MIN(ratio0, ratio1);
    return REF_SUCCESS;
  }

  r_min = MIN(ratio0, ratio1);
  r_max = MAX(ratio0, ratio1);

  r = r_min / r_max;

  if (ABS(r - 1.0) < 1.0e-12) {
    *ratio = 0.5 * (ratio0 + ratio1);
    return REF_SUCCESS;
  }

  *ratio = r_min * (r - 1.0) / (r * log(r));

  return REF_SUCCESS;
}

REF_STATUS ref_node_dratio_dnode0(REF_NODE ref_node, REF_INT node0,
                                  REF_INT node1, REF_DBL *ratio,
                                  REF_DBL *d_ratio) {
  REF_DBL direction[3], length;
  REF_DBL ratio0, d_ratio0[3], ratio1, d_ratio1[3];
  REF_DBL r, d_r[3], r_min, d_r_min[3], r_max, d_r_max[3];
  REF_INT i;
  REF_DBL r_log_r;
  REF_DBL m[6];

  if (!ref_node_valid(ref_node, node0) || !ref_node_valid(ref_node, node1))
    RSS(REF_INVALID, "node invalid");

  direction[0] =
      (ref_node_xyz(ref_node, 0, node1) - ref_node_xyz(ref_node, 0, node0));
  direction[1] =
      (ref_node_xyz(ref_node, 1, node1) - ref_node_xyz(ref_node, 1, node0));
  direction[2] =
      (ref_node_xyz(ref_node, 2, node1) - ref_node_xyz(ref_node, 2, node0));

  length = ref_math_dot(direction, direction);
  length = sqrt(length);

  if (!ref_math_divisible(direction[0], length) ||
      !ref_math_divisible(direction[1], length) ||
      !ref_math_divisible(direction[2], length)) {
    *ratio = 0.0;
    d_ratio[0] = 0.0;
    d_ratio[1] = 0.0;
    d_ratio[2] = 0.0;
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get(ref_node, node0, m), "node0 m");
  RSS(ref_matrix_sqrt_vt_m_v_deriv(m, direction, &ratio0, d_ratio0), "vt m v0");
  for (i = 0; i < 3; i++) d_ratio0[i] = -d_ratio0[i]; /* node 0 is neg */
  RSS(ref_node_metric_get(ref_node, node1, m), "node1 m");
  RSS(ref_matrix_sqrt_vt_m_v_deriv(m, direction, &ratio1, d_ratio1), "vt m v0");
  for (i = 0; i < 3; i++) d_ratio1[i] = -d_ratio1[i]; /* node 0 is neg */

  /* Loseille Lohner IMR 18 (2009) pg 613 */
  /* Alauzet Finite Elements in Analysis and Design 46 (2010) pg 185 */

  if (ratio0 < 1.0e-12 || ratio1 < 1.0e-12) {
    if (ratio0 < ratio1) {
      *ratio = ratio0;
      d_ratio[0] = d_ratio0[0];
      d_ratio[1] = d_ratio0[1];
      d_ratio[2] = d_ratio0[2];
    } else {
      *ratio = ratio1;
      d_ratio[0] = d_ratio1[0];
      d_ratio[1] = d_ratio1[1];
      d_ratio[2] = d_ratio1[2];
    }
    return REF_SUCCESS;
  }

  if (ratio0 < ratio1) {
    r_min = ratio0;
    d_r_min[0] = d_ratio0[0];
    d_r_min[1] = d_ratio0[1];
    d_r_min[2] = d_ratio0[2];
    r_max = ratio1;
    d_r_max[0] = d_ratio1[0];
    d_r_max[1] = d_ratio1[1];
    d_r_max[2] = d_ratio1[2];
  } else {
    r_min = ratio1;
    d_r_min[0] = d_ratio1[0];
    d_r_min[1] = d_ratio1[1];
    d_r_min[2] = d_ratio1[2];
    r_max = ratio0;
    d_r_max[0] = d_ratio0[0];
    d_r_max[1] = d_ratio0[1];
    d_r_max[2] = d_ratio0[2];
  }

  r = r_min / r_max;
  for (i = 0; i < 3; i++)
    d_r[i] = (d_r_min[i] * r_max - r_min * d_r_max[i]) / r_max / r_max;

  if (ABS(r - 1.0) < 1.0e-12) {
    *ratio = 0.5 * (r_min + r_max);
    for (i = 0; i < 3; i++) d_ratio[i] = 0.5 * (d_r_min[i] + d_r_max[i]);
    return REF_SUCCESS;
  }

  r_log_r = r * log(r);
  *ratio = r_min * (r - 1.0) / r_log_r;

  for (i = 0; i < 3; i++)
    d_ratio[i] = ((r_min * d_r[i] + d_r_min[i] * (r - 1.0)) * r_log_r -
                  r_min * (r - 1.0) * (r * 1 / r * d_r[i] + d_r[i] * log(r))) /
                 r_log_r / r_log_r;

  return REF_SUCCESS;
}

REF_STATUS ref_node_ratio_node0(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                                REF_DBL *ratio_node0) {
  REF_DBL direction[3], length;
  REF_DBL m[6];

  if (!ref_node_valid(ref_node, node0) || !ref_node_valid(ref_node, node1))
    RSS(REF_INVALID, "node invalid");

  direction[0] =
      (ref_node_xyz(ref_node, 0, node1) - ref_node_xyz(ref_node, 0, node0));
  direction[1] =
      (ref_node_xyz(ref_node, 1, node1) - ref_node_xyz(ref_node, 1, node0));
  direction[2] =
      (ref_node_xyz(ref_node, 2, node1) - ref_node_xyz(ref_node, 2, node0));

  length = ref_math_dot(direction, direction);
  length = sqrt(length);

  if (!ref_math_divisible(direction[0], length) ||
      !ref_math_divisible(direction[1], length) ||
      !ref_math_divisible(direction[2], length)) {
    *ratio_node0 = 0.0;
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get(ref_node, node0, m), "node0 m");
  *ratio_node0 = ref_matrix_sqrt_vt_m_v(m, direction);

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_epic_quality(REF_NODE ref_node, REF_INT *nodes,
                                     REF_DBL *quality) {
  REF_DBL l0, l1, l2, l3, l4, l5;

  REF_DBL det, min_det, volume;
  REF_DBL volume_in_metric;
  REF_DBL num, denom;
  REF_DBL m[6];

  RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");

  if (volume <= ref_node_min_volume(ref_node)) {
    *quality = volume - ref_node_min_volume(ref_node);
    return REF_SUCCESS;
  }

  RSS(ref_node_ratio(ref_node, nodes[0], nodes[1], &l0), "l0");
  RSS(ref_node_ratio(ref_node, nodes[0], nodes[2], &l1), "l1");
  RSS(ref_node_ratio(ref_node, nodes[0], nodes[3], &l2), "l2");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[2], &l3), "l3");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[3], &l4), "l4");
  RSS(ref_node_ratio(ref_node, nodes[2], nodes[3], &l5), "l5");

  RSS(ref_node_metric_get(ref_node, nodes[0], m), "nodes[0] m");
  RSS(ref_matrix_det_m(m, &det), "n0");
  min_det = det;

  RSS(ref_node_metric_get(ref_node, nodes[1], m), "nodes[1] m");
  RSS(ref_matrix_det_m(m, &det), "n1");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[2], m), "nodes[2] m");
  RSS(ref_matrix_det_m(m, &det), "n2");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[3], m), "nodes[3] m");
  RSS(ref_matrix_det_m(m, &det), "n3");
  min_det = MIN(min_det, det);

  volume_in_metric = sqrt(min_det) * volume;

  num = pow(volume_in_metric, 2.0 / 3.0);
  denom = l0 * l0 + l1 * l1 + l2 * l2 + l3 * l3 + l4 * l4 + l5 * l5;

  if (ref_math_divisible(num, denom)) {
    /* 36/3^(1/3) */
    *quality = 24.9610058766228 * num / denom;
  } else {
    /* printf("%s: %d: %s: div zero vol %.18e min_det %.18e (%.18e / %.18e)\n",
       __FILE__, __LINE__, __func__, volume, min_det, num, denom); */
    *quality = -1.0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_epic_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                             REF_DBL *quality,
                                             REF_DBL *d_quality) {
  REF_DBL l0, l1, l2, l3, l4, l5;
  REF_DBL d_l0[3], d_l1[3], d_l2[3];

  REF_DBL det, min_det, volume, d_volume[3];
  REF_DBL volume_in_metric, d_volume_in_metric[3];
  REF_DBL num, denom;
  REF_DBL d_num[3], d_denom[3];
  REF_INT i;
  REF_DBL m[6];

  RSS(ref_node_dratio_dnode0(ref_node, nodes[0], nodes[1], &l0, d_l0), "l0");
  RSS(ref_node_dratio_dnode0(ref_node, nodes[0], nodes[2], &l1, d_l1), "l1");
  RSS(ref_node_dratio_dnode0(ref_node, nodes[0], nodes[3], &l2, d_l2), "l2");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[2], &l3), "l3");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[3], &l4), "l4");
  RSS(ref_node_ratio(ref_node, nodes[2], nodes[3], &l5), "l5");

  RSS(ref_node_tet_dvol_dnode0(ref_node, nodes, &volume, d_volume), "vol");

  if (volume <= ref_node_min_volume(ref_node)) {
    *quality = volume - ref_node_min_volume(ref_node);
    for (i = 0; i < 3; i++) d_quality[i] = d_volume[i];
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get(ref_node, nodes[0], m), "nodes[0] m");
  RSS(ref_matrix_det_m(m, &det), "n0");
  min_det = det;

  RSS(ref_node_metric_get(ref_node, nodes[1], m), "nodes[1] m");
  RSS(ref_matrix_det_m(m, &det), "n1");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[2], m), "nodes[2] m");
  RSS(ref_matrix_det_m(m, &det), "n2");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[3], m), "nodes[3] m");
  RSS(ref_matrix_det_m(m, &det), "n3");
  min_det = MIN(min_det, det);

  volume_in_metric = sqrt(min_det) * volume;
  for (i = 0; i < 3; i++) d_volume_in_metric[i] = sqrt(min_det) * d_volume[i];

  num = pow(volume_in_metric, 2.0 / 3.0);
  for (i = 0; i < 3; i++)
    d_num[i] =
        2.0 / 3.0 * pow(volume_in_metric, -1.0 / 3.0) * d_volume_in_metric[i];
  denom = l0 * l0 + l1 * l1 + l2 * l2 + l3 * l3 + l4 * l4 + l5 * l5;
  for (i = 0; i < 3; i++)
    d_denom[i] = 2.0 * l0 * d_l0[i] + 2.0 * l1 * d_l1[i] + 2.0 * l2 * d_l2[i];

  if (ref_math_divisible(num, denom)) {
    /* 36/3^(1/3) */
    *quality = 24.9610058766228 * num / denom;
    for (i = 0; i < 3; i++)
      d_quality[i] = 24.9610058766228 * (d_num[i] * denom - num * d_denom[i]) /
                     denom / denom;
  } else {
    /* printf("%s: %d: %s: div zero vol %.18e min_det %.18e (%.18e / %.18e)\n",
       __FILE__, __LINE__, __func__, volume, min_det, num, denom); */
    *quality = -1.0;
    for (i = 0; i < 3; i++) d_quality[i] = 0.0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_jac_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                            REF_DBL *quality,
                                            REF_DBL *d_quality) {
  REF_DBL mlog0[6], mlog1[6], mlog2[6], mlog3[6];
  REF_DBL mlog[6], m[6], jac[9];
  REF_DBL e0[3], e1[3], e2[3], e3[3], e4[3], e5[3];
  REF_INT i;

  REF_DBL l2, det, volume, volume_in_metric, num;
  REF_DBL d_volume[3];
  REF_DBL temp, d_e0[3], d_e1[3], d_e2[3], d_l2[3], d_num[3];
  REF_DBL sqrt_det, pow_vim;

  RSS(ref_node_tet_dvol_dnode0(ref_node, nodes, &volume, d_volume), "vol");
  if (volume <= ref_node_min_volume(ref_node)) {
    *quality = volume - ref_node_min_volume(ref_node);
    for (i = 0; i < 3; i++) d_quality[i] = d_volume[i];
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get_log(ref_node, nodes[0], mlog0), "log0");
  RSS(ref_node_metric_get_log(ref_node, nodes[1], mlog1), "log1");
  RSS(ref_node_metric_get_log(ref_node, nodes[2], mlog2), "log2");
  RSS(ref_node_metric_get_log(ref_node, nodes[3], mlog3), "log3");
  for (i = 0; i < 6; i++)
    mlog[i] = (mlog0[i] + mlog1[i] + mlog2[i] + mlog3[i]) / 4.0;
  RSS(ref_matrix_exp_m(mlog, m), "exp");
  RSS(ref_matrix_jacob_m(m, jac), "jac");

  for (i = 0; i < 3; i++)
    e0[i] = ref_node_xyz(ref_node, i, nodes[1]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e1[i] = ref_node_xyz(ref_node, i, nodes[2]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e2[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e3[i] = ref_node_xyz(ref_node, i, nodes[2]) -
            ref_node_xyz(ref_node, i, nodes[1]);
  for (i = 0; i < 3; i++)
    e4[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[1]);
  for (i = 0; i < 3; i++)
    e5[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[2]);

  l2 = ref_matrix_vt_m_v(m, e0) + ref_matrix_vt_m_v(m, e1) +
       ref_matrix_vt_m_v(m, e2) + ref_matrix_vt_m_v(m, e3) +
       ref_matrix_vt_m_v(m, e4) + ref_matrix_vt_m_v(m, e5);

  RSS(ref_matrix_vt_m_v_deriv(m, e0, &temp, d_e0), "d_e0");
  RSS(ref_matrix_vt_m_v_deriv(m, e1, &temp, d_e1), "d_e1");
  RSS(ref_matrix_vt_m_v_deriv(m, e2, &temp, d_e2), "d_e2");
  d_l2[0] = -d_e0[0] - d_e1[0] - d_e2[0];
  d_l2[1] = -d_e0[1] - d_e1[1] - d_e2[1];
  d_l2[2] = -d_e0[2] - d_e1[2] - d_e2[2];

  RSS(ref_matrix_det_m(m, &det), "det(mavg)");
  sqrt_det = sqrt(det);
  volume_in_metric = sqrt_det * volume;

  num = pow(volume_in_metric, 2.0 / 3.0);

  pow_vim = pow(volume_in_metric, -1.0 / 3.0);
  d_num[0] = (2.0 / 3.0) * pow_vim * sqrt_det * d_volume[0];
  d_num[1] = (2.0 / 3.0) * pow_vim * sqrt_det * d_volume[1];
  d_num[2] = (2.0 / 3.0) * pow_vim * sqrt_det * d_volume[2];

  if (ref_math_divisible(num, l2)) {
    /* 36/3^(1/3) */
    *quality = 24.9610058766228 * num / l2;
    d_quality[0] =
        24.9610058766228 * (d_num[0] * l2 - num * d_l2[0]) / (l2 * l2);
    d_quality[1] =
        24.9610058766228 * (d_num[1] * l2 - num * d_l2[1]) / (l2 * l2);
    d_quality[2] =
        24.9610058766228 * (d_num[2] * l2 - num * d_l2[2]) / (l2 * l2);
  } else {
    /* printf("%s: %d: %s: div zero vol %.18e (%.18e / %.18e)\n", __FILE__,
       __LINE__, __func__, volume, num, l2); */
    *quality = -1.0;
    d_quality[0] = 0.0;
    d_quality[1] = 0.0;
    d_quality[2] = 0.0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                        REF_DBL *quality, REF_DBL *d_quality) {
  switch (ref_node->tet_quality) {
    case REF_NODE_EPIC_QUALITY:
      RSS(ref_node_tet_epic_dquality_dnode0(ref_node, nodes, quality,
                                            d_quality),
          "epic");
      break;
    case REF_NODE_JAC_QUALITY:
      RSS(ref_node_tet_jac_dquality_dnode0(ref_node, nodes, quality, d_quality),
          "jac");
      break;
    default:
      THROW("case not recognized");
      break;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_jac_quality(REF_NODE ref_node, REF_INT *nodes,
                                    REF_DBL *quality) {
  REF_DBL mlog0[6], mlog1[6], mlog2[6], mlog3[6];
  REF_DBL mlog[6], m[6], jac[9];
  REF_DBL e0[3], e1[3], e2[3], e3[3], e4[3], e5[3];
  REF_INT i;

  REF_DBL l2, det, volume, volume_in_metric, num;

  RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");
  if (volume <= ref_node_min_volume(ref_node)) {
    *quality = volume - ref_node_min_volume(ref_node);
    return REF_SUCCESS;
  }

  RSS(ref_node_metric_get_log(ref_node, nodes[0], mlog0), "log0");
  RSS(ref_node_metric_get_log(ref_node, nodes[1], mlog1), "log1");
  RSS(ref_node_metric_get_log(ref_node, nodes[2], mlog2), "log2");
  RSS(ref_node_metric_get_log(ref_node, nodes[3], mlog3), "log3");
  for (i = 0; i < 6; i++)
    mlog[i] = (mlog0[i] + mlog1[i] + mlog2[i] + mlog3[i]) / 4.0;
  RSS(ref_matrix_exp_m(mlog, m), "exp");
  RSS(ref_matrix_jacob_m(m, jac), "jac");

  for (i = 0; i < 3; i++)
    e0[i] = ref_node_xyz(ref_node, i, nodes[1]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e1[i] = ref_node_xyz(ref_node, i, nodes[2]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e2[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    e3[i] = ref_node_xyz(ref_node, i, nodes[2]) -
            ref_node_xyz(ref_node, i, nodes[1]);
  for (i = 0; i < 3; i++)
    e4[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[1]);
  for (i = 0; i < 3; i++)
    e5[i] = ref_node_xyz(ref_node, i, nodes[3]) -
            ref_node_xyz(ref_node, i, nodes[2]);

  l2 = ref_matrix_vt_m_v(m, e0) + ref_matrix_vt_m_v(m, e1) +
       ref_matrix_vt_m_v(m, e2) + ref_matrix_vt_m_v(m, e3) +
       ref_matrix_vt_m_v(m, e4) + ref_matrix_vt_m_v(m, e5);

  RSS(ref_matrix_det_m(m, &det), "det(mavg)");
  volume_in_metric = sqrt(det) * volume;

  num = pow(volume_in_metric, 2.0 / 3.0);

  if (ref_math_divisible(num, l2)) {
    /* 36/3^(1/3) */
    *quality = 24.9610058766228 * num / l2;
  } else {
    /* printf("%s: %d: %s: div zero vol %.18e (%.18e / %.18e)\n", __FILE__,
       __LINE__, __func__, volume, num, l2); */
    *quality = -1.0;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_node_tet_quality(REF_NODE ref_node, REF_INT *nodes,
                                REF_DBL *quality) {
  if (REF_FALSE) {
    REF_DBL epic, jac;
    RSS(ref_node_tet_epic_quality(ref_node, nodes, &epic), "epic");
    RSS(ref_node_tet_jac_quality(ref_node, nodes, &jac), "epic");
    printf("tet epic %11.8f jac %11.8f\n", epic, jac);
  }

  switch (ref_node->tet_quality) {
    case REF_NODE_EPIC_QUALITY:
      RSS(ref_node_tet_epic_quality(ref_node, nodes, quality), "epic");
      break;
    case REF_NODE_JAC_QUALITY:
      RSS(ref_node_tet_jac_quality(ref_node, nodes, quality), "jac");
      break;
    default:
      THROW("case not recognized");
      break;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_node_tri_epic_quality(REF_NODE ref_node, REF_INT *nodes,
                                            REF_DBL *quality) {
  REF_DBL l0, l1, l2;

  REF_DBL det, min_det, area;
  REF_DBL area_in_metric;
  REF_DBL num, denom;
  REF_DBL m[6];

  RSS(ref_node_ratio(ref_node, nodes[0], nodes[1], &l0), "l0");
  RSS(ref_node_ratio(ref_node, nodes[0], nodes[2], &l1), "l1");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[2], &l2), "l2");

  RSS(ref_node_tri_area(ref_node, nodes, &area), "area");

  RSS(ref_node_metric_get(ref_node, nodes[0], m), "nodes[0] m");
  RSS(ref_matrix_det_m(m, &det), "n0");
  min_det = det;

  RSS(ref_node_metric_get(ref_node, nodes[1], m), "nodes[1] m");
  RSS(ref_matrix_det_m(m, &det), "n1");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[2], m), "nodes[2] m");
  RSS(ref_matrix_det_m(m, &det), "n2");
  min_det = MIN(min_det, det);

  area_in_metric = pow(min_det, 1.0 / 3.0) * area;

  num = area_in_metric;
  denom = l0 * l0 + l1 * l1 + l2 * l2;

  if (ref_math_divisible(num, denom)) {
    *quality = 4.0 / sqrt(3.0) * 3 * num / denom;
  } else {
    /* printf("%s: %d: %s: div zero area %.18e min_det %.18e (%.18e / %.18e)\n",
       __FILE__, __LINE__, __func__, area, min_det, num, denom); */
    *quality = -1.0;
  }

  return REF_SUCCESS;
}
static REF_STATUS ref_node_tri_jac_quality(REF_NODE ref_node, REF_INT *nodes,
                                           REF_DBL *quality) {
  REF_DBL mlog0[6], mlog1[6], mlog2[6];
  REF_DBL mlog[6], m[6], jac[9];
  REF_DBL xyz0[3], xyz1[3], xyz2[3];
  REF_DBL e0[3], e1[3], e2[3], n[3];
  REF_DBL a, l2;
  REF_INT i;

  RSS(ref_node_metric_get_log(ref_node, nodes[0], mlog0), "log0");
  RSS(ref_node_metric_get_log(ref_node, nodes[1], mlog1), "log1");
  RSS(ref_node_metric_get_log(ref_node, nodes[2], mlog2), "log2");
  for (i = 0; i < 6; i++) mlog[i] = (mlog0[i] + mlog1[i] + mlog2[i]) / 3.0;
  RSS(ref_matrix_exp_m(mlog, m), "exp");
  RSS(ref_matrix_jacob_m(m, jac), "jac");

  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[0]), xyz0),
      "xyz0");
  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[1]), xyz1),
      "xyz1");
  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[2]), xyz2),
      "xyz2");

  for (i = 0; i < 3; i++) e0[i] = xyz2[i] - xyz1[i];
  for (i = 0; i < 3; i++) e1[i] = xyz0[i] - xyz2[i];
  for (i = 0; i < 3; i++) e2[i] = xyz1[i] - xyz0[i];

  ref_math_cross_product(e2, e0, n);
  l2 = ref_math_dot(e0, e0) + ref_math_dot(e1, e1) + ref_math_dot(e2, e2);

  a = 0.5 * sqrt(ref_math_dot(n, n));

  if (ref_math_divisible(a, l2)) {
    *quality = 4.0 * sqrt(3.0) * (a / l2);
  } else {
    /* printf("%s: %d: %s: div zero area %.18e l2 %.18e\n", __FILE__, __LINE__,
       __func__, a, l2); */
    *quality = -1.0;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_node_tri_quality(REF_NODE ref_node, REF_INT *nodes,
                                REF_DBL *quality) {
  if (REF_FALSE) {
    REF_DBL epic, jac;
    RSS(ref_node_tri_epic_quality(ref_node, nodes, &epic), "epic");
    RSS(ref_node_tri_jac_quality(ref_node, nodes, &jac), "epic");
    printf("tri epic %11.8f jac %11.8f\n", epic, jac);
  }

  switch (ref_node->tri_quality) {
    case REF_NODE_EPIC_QUALITY:
      RSS(ref_node_tri_epic_quality(ref_node, nodes, quality), "epic");
      break;
    case REF_NODE_JAC_QUALITY:
      RSS(ref_node_tri_jac_quality(ref_node, nodes, quality), "epic");
      break;
    default:
      THROW("case not recognized");
      break;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_node_tri_epic_dquality_dnode0(REF_NODE ref_node,
                                                    REF_INT *nodes,
                                                    REF_DBL *quality,
                                                    REF_DBL *d_quality) {
  REF_DBL l0, l1, l2;

  REF_DBL det, min_det, area, d_area[3];
  REF_DBL area_in_metric, d_area_in_metric[3];
  REF_DBL num, d_num[3], denom, d_denom[3];
  REF_DBL d_l0[3], d_l1[3];
  REF_INT i;
  REF_DBL m[6];

  RSS(ref_node_dratio_dnode0(ref_node, nodes[0], nodes[1], &l0, d_l0), "l0");
  RSS(ref_node_dratio_dnode0(ref_node, nodes[0], nodes[2], &l1, d_l1), "l1");
  RSS(ref_node_ratio(ref_node, nodes[1], nodes[2], &l2), "l2");

  RSS(ref_node_tri_darea_dnode0(ref_node, nodes, &area, d_area), "area");

  RSS(ref_node_metric_get(ref_node, nodes[0], m), "nodes[0] m");
  RSS(ref_matrix_det_m(m, &det), "n0");
  min_det = det;

  RSS(ref_node_metric_get(ref_node, nodes[1], m), "nodes[1] m");
  RSS(ref_matrix_det_m(m, &det), "n1");
  min_det = MIN(min_det, det);

  RSS(ref_node_metric_get(ref_node, nodes[2], m), "nodes[2] m");
  RSS(ref_matrix_det_m(m, &det), "n2");
  min_det = MIN(min_det, det);

  area_in_metric = pow(min_det, 1.0 / 3.0) * area;
  for (i = 0; i < 3; i++)
    d_area_in_metric[i] = pow(min_det, 1.0 / 3.0) * d_area[i];

  num = area_in_metric;
  for (i = 0; i < 3; i++) d_num[i] = d_area_in_metric[i];
  denom = l0 * l0 + l1 * l1 + l2 * l2;
  for (i = 0; i < 3; i++) d_denom[i] = 2.0 * l0 * d_l0[i] + 2.0 * l1 * d_l1[i];

  if (ref_math_divisible(num, denom)) {
    *quality = 4.0 / sqrt(3.0) * 3 * num / denom;
    for (i = 0; i < 3; i++)
      d_quality[i] = 4.0 / sqrt(3.0) * 3 *
                     (d_num[i] * denom - num * d_denom[i]) / denom / denom;
  } else {
    /* printf("%s: %d: %s: div zero area %.18e min_det %.18e (%.18e / %.18e)\n",
       __FILE__, __LINE__, __func__, area, min_det, num, denom); */
    *quality = -1.0;
    for (i = 0; i < 3; i++) d_quality[i] = 0.0;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_node_tri_jac_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                            REF_DBL *quality,
                                            REF_DBL *d_quality) {
  REF_DBL mlog0[6], mlog1[6], mlog2[6];
  REF_DBL mlog[6], m[6], jac[9];
  REF_DBL xyz0[3], xyz1[3], xyz2[3];
  REF_DBL dxyz0[3][3];
  REF_DBL e0[3], e1[3], e2[3], n[3];
  REF_DBL de1[3][3], de2[3][3], dn[3][3];
  REF_DBL a, l2, da[3], dl2[3];
  REF_INT i, j;

  RSS(ref_node_metric_get_log(ref_node, nodes[0], mlog0), "log0");
  RSS(ref_node_metric_get_log(ref_node, nodes[1], mlog1), "log1");
  RSS(ref_node_metric_get_log(ref_node, nodes[2], mlog2), "log2");
  for (i = 0; i < 6; i++) mlog[i] = (mlog0[i] + mlog1[i] + mlog2[i]) / 3.0;
  RSS(ref_matrix_exp_m(mlog, m), "exp");
  RSS(ref_matrix_jacob_m(m, jac), "jac");

  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[0]), xyz0),
      "xyz0");
  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[1]), xyz1),
      "xyz1");
  RSS(ref_matrix_vect_mult(jac, ref_node_xyz_ptr(ref_node, nodes[2]), xyz2),
      "xyz2");

  dxyz0[0][0] = jac[0];
  dxyz0[0][1] = jac[1];
  dxyz0[0][2] = jac[2];

  dxyz0[1][0] = jac[3];
  dxyz0[1][1] = jac[4];
  dxyz0[1][2] = jac[5];

  dxyz0[2][0] = jac[6];
  dxyz0[2][1] = jac[7];
  dxyz0[2][2] = jac[8];

  for (i = 0; i < 3; i++) e0[i] = xyz2[i] - xyz1[i];
  for (i = 0; i < 3; i++) e1[i] = xyz0[i] - xyz2[i];
  for (i = 0; i < 3; i++) e2[i] = xyz1[i] - xyz0[i];

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) de1[i][j] = dxyz0[i][j];
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) de2[i][j] = -dxyz0[i][j];

  ref_math_cross_product(e2, e0, n);
  for (j = 0; j < 3; j++) {
    dn[0][j] = de2[1][j] * e0[2] - de2[2][j] * e0[1];
    dn[1][j] = de2[2][j] * e0[0] - de2[0][j] * e0[2];
    dn[2][j] = de2[0][j] * e0[1] - de2[1][j] * e0[0];
  }

  l2 = ref_math_dot(e0, e0) + ref_math_dot(e1, e1) + ref_math_dot(e2, e2);

  for (j = 0; j < 3; j++)
    dl2[j] = 2.0 * e1[0] * de1[0][j] + 2.0 * e1[1] * de1[1][j] +
             2.0 * e1[2] * de1[2][j] + 2.0 * e2[0] * de2[0][j] +
             2.0 * e2[1] * de2[1][j] + 2.0 * e2[2] * de2[2][j];

  a = 0.5 * sqrt(ref_math_dot(n, n));

  for (j = 0; j < 3; j++)
    da[j] =
        0.5 * 0.5 / sqrt(ref_math_dot(n, n)) *
        (2.0 * n[0] * dn[0][j] + 2.0 * n[1] * dn[1][j] + 2.0 * n[2] * dn[2][j]);

  if (ref_math_divisible(a, l2)) {
    *quality = 4.0 * sqrt(3.0) * (a / l2);
    for (j = 0; j < 3; j++)
      d_quality[j] = 4.0 * sqrt(3.0) * (da[j] * l2 - a * dl2[j]) / l2 / l2;
  } else {
    /* printf("%s: %d: %s: div zero area %.18e l2 %.18e\n", __FILE__, __LINE__,
       __func__, a, l2); */
    *quality = -1.0;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_node_tri_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                        REF_DBL *quality, REF_DBL *d_quality) {
  switch (ref_node->tri_quality) {
    case REF_NODE_EPIC_QUALITY:
      RSS(ref_node_tri_epic_dquality_dnode0(ref_node, nodes, quality,
                                            d_quality),
          "epic");
      break;
    case REF_NODE_JAC_QUALITY:
      RSS(ref_node_tri_jac_dquality_dnode0(ref_node, nodes, quality, d_quality),
          "epic");
      break;
    default:
      THROW("case not recognized");
      break;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_xyz_normal(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz2,
                               REF_DBL *normal) {
  REF_DBL edge10[3], edge20[3];

  edge10[0] = xyz1[0] - xyz0[0];
  edge10[1] = xyz1[1] - xyz0[1];
  edge10[2] = xyz1[2] - xyz0[2];

  edge20[0] = xyz2[0] - xyz0[0];
  edge20[1] = xyz2[1] - xyz0[1];
  edge20[2] = xyz2[2] - xyz0[2];

  ref_math_cross_product(edge10, edge20, normal);

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_normal(REF_NODE ref_node, REF_INT *nodes,
                               REF_DBL *normal) {
  REF_DBL *xyz0, *xyz1, *xyz2;

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]))
    RSS(REF_INVALID, "node invalid");

  xyz0 = ref_node_xyz_ptr(ref_node, nodes[0]);
  xyz1 = ref_node_xyz_ptr(ref_node, nodes[1]);
  xyz2 = ref_node_xyz_ptr(ref_node, nodes[2]);

  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz2, normal), "xyz norm");

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_centroid(REF_NODE ref_node, REF_INT *nodes,
                                 REF_DBL *centroid) {
  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]))
    RSS(REF_INVALID, "node invalid");

  centroid[0] = (ref_node_xyz(ref_node, 0, nodes[0]) +
                 ref_node_xyz(ref_node, 0, nodes[1]) +
                 ref_node_xyz(ref_node, 0, nodes[2])) /
                3.0;
  centroid[1] = (ref_node_xyz(ref_node, 1, nodes[0]) +
                 ref_node_xyz(ref_node, 1, nodes[1]) +
                 ref_node_xyz(ref_node, 1, nodes[2])) /
                3.0;
  centroid[2] = (ref_node_xyz(ref_node, 2, nodes[0]) +
                 ref_node_xyz(ref_node, 2, nodes[1]) +
                 ref_node_xyz(ref_node, 2, nodes[2])) /
                3.0;

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_y_projection(REF_NODE ref_node, REF_INT *nodes,
                                     REF_DBL *y_projection) {
  REF_DBL normal[3];

  RSS(ref_node_tri_normal(ref_node, nodes, normal), "norm inside of area");

  *y_projection = 0.5 * normal[1];

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_twod_orientation(REF_NODE ref_node, REF_INT *nodes,
                                         REF_BOOL *valid) {
  REF_DBL mid_plane = ref_node_twod_mid_plane(ref_node);
  REF_DBL normal[3];

  *valid = REF_FALSE;

  RSS(ref_node_tri_normal(ref_node, nodes, normal), "norm inside of area");

  if ((ref_node_xyz(ref_node, 1, nodes[0]) > mid_plane && normal[1] < 0.0) ||
      (ref_node_xyz(ref_node, 1, nodes[0]) < mid_plane && normal[1] > 0.0))
    *valid = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_node_angle(REF_NODE ref_node, REF_INT *nodes,
                                   REF_INT node, REF_DBL *angle) {
  REF_INT node1, node2, i;
  REF_DBL edge1[3], edge2[3];

  *angle = -9.0;

  node1 = REF_EMPTY;
  node2 = REF_EMPTY;
  if (node == nodes[0]) {
    node1 = nodes[1];
    node2 = nodes[2];
  }
  if (node == nodes[1]) {
    node1 = nodes[2];
    node2 = nodes[0];
  }
  if (node == nodes[2]) {
    node1 = nodes[0];
    node2 = nodes[1];
  }
  if (REF_EMPTY == node1 || REF_EMPTY == node2) return REF_NOT_FOUND;

  for (i = 0; i < 3; i++) {
    edge1[i] =
        ref_node_xyz(ref_node, i, node1) - ref_node_xyz(ref_node, i, node);
    edge2[i] =
        ref_node_xyz(ref_node, i, node2) - ref_node_xyz(ref_node, i, node);
  }

  RSS(ref_math_normalize(edge1), "normalize zero length edge1");
  RSS(ref_math_normalize(edge2), "normalize zero length edge2");
  *angle = acos(ref_math_dot(edge1, edge2));

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_area(REF_NODE ref_node, REF_INT *nodes, REF_DBL *area) {
  REF_DBL normal[3];

  RSS(ref_node_tri_normal(ref_node, nodes, normal), "norm inside of area");

  *area = 0.5 * sqrt(ref_math_dot(normal, normal));

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_darea_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                     REF_DBL *area, REF_DBL *d_area) {
  REF_DBL *xyz0, *xyz1, *xyz2;
  REF_DBL v0[3], v1[3];
  REF_DBL normx, normy, normz;
  REF_DBL d_normx[3], d_normy[3], d_normz[3];
  REF_INT i;

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]))
    RSS(REF_INVALID, "node invalid");

  xyz0 = ref_node_xyz_ptr(ref_node, nodes[0]);
  xyz1 = ref_node_xyz_ptr(ref_node, nodes[1]);
  xyz2 = ref_node_xyz_ptr(ref_node, nodes[2]);

  v0[0] = xyz1[0] - xyz0[0];
  v0[1] = xyz1[1] - xyz0[1];
  v0[2] = xyz1[2] - xyz0[2];

  v1[0] = xyz2[0] - xyz0[0];
  v1[1] = xyz2[1] - xyz0[1];
  v1[2] = xyz2[2] - xyz0[2];

  normx = (v0)[1] * (v1)[2] - (v0)[2] * (v1)[1];
  d_normx[0] = 0.0;
  d_normx[1] = -(v1)[2] + (v0)[2];
  d_normx[2] = -(v0)[1] + (v1)[1];

  normy = (v0)[2] * (v1)[0] - (v0)[0] * (v1)[2];
  d_normy[0] = -(v0)[2] + (v1)[2];
  d_normy[1] = 0.0;
  d_normy[2] = -(v1)[0] + (v0)[0];

  normz = (v0)[0] * (v1)[1] - (v0)[1] * (v1)[0];
  d_normz[0] = -(v1)[1] + (v0)[1];
  d_normz[1] = -(v0)[0] + (v1)[0];
  d_normz[2] = 0.0;

  *area = 0.5 * sqrt(normx * normx + normy * normy + normz * normz);
  for (i = 0; i < 3; i++)
    d_area[i] = 0.5 * 0.5 /
                sqrt(normx * normx + normy * normy + normz * normz) *
                (2.0 * normx * d_normx[i] + 2.0 * normy * d_normy[i] +
                 2.0 * normz * d_normz[i]);

  return REF_SUCCESS;
}

REF_STATUS ref_node_xyz_vol(REF_DBL *xyzs[4], REF_DBL *volume) {
  REF_DBL *a, *b, *c, *d;
  REF_DBL m11, m12, m13;
  REF_DBL det;

  a = xyzs[0];
  b = xyzs[1];
  c = xyzs[2];
  d = xyzs[3];

  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  det = (m11 - m12 + m13);

  *volume = -det / 6.0;

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_vol(REF_NODE ref_node, REF_INT *nodes,
                            REF_DBL *volume) {
  REF_DBL *a, *b, *c, *d;
  REF_DBL m11, m12, m13;
  REF_DBL det;

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]) ||
      !ref_node_valid(ref_node, nodes[3]))
    RSS(REF_INVALID, "node invalid");

  a = ref_node_xyz_ptr(ref_node, nodes[0]);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);
  c = ref_node_xyz_ptr(ref_node, nodes[2]);
  d = ref_node_xyz_ptr(ref_node, nodes[3]);

  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  det = (m11 - m12 + m13);

  *volume = -det / 6.0;

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_dvol_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                    REF_DBL *vol, REF_DBL *d_vol) {
  REF_DBL *a, *b, *c, *d;
  REF_DBL m11, m12, m13;
  REF_DBL det;

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]) ||
      !ref_node_valid(ref_node, nodes[3]))
    RSS(REF_INVALID, "node invalid");

  a = ref_node_xyz_ptr(ref_node, nodes[0]);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);
  c = ref_node_xyz_ptr(ref_node, nodes[2]);
  d = ref_node_xyz_ptr(ref_node, nodes[3]);

  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  det = (m11 - m12 + m13);

  *vol = -det / 6.0;
  d_vol[0] =
      -((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2])) / 6.0;
  d_vol[1] =
      ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2])) / 6.0;
  d_vol[2] =
      -((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1])) / 6.0;

  return REF_SUCCESS;
}

REF_STATUS ref_node_twod_clone(REF_NODE ref_node, REF_INT original,
                               REF_INT *clone_ptr) {
  REF_DBL mid_plane = ref_node_twod_mid_plane(ref_node);
  REF_INT global, clone;
  REF_INT i;
  REF_DBL m[6];

  RSS(ref_node_next_global(ref_node, &global), "next global");
  RSS(ref_node_add(ref_node, global, clone_ptr), "new node");
  clone = *clone_ptr;

  ref_node_xyz(ref_node, 0, clone) = ref_node_xyz(ref_node, 0, original);
  ref_node_xyz(ref_node, 1, clone) =
      2 * mid_plane - ref_node_xyz(ref_node, 1, original);
  ref_node_xyz(ref_node, 2, clone) = ref_node_xyz(ref_node, 2, original);

  for (i = 0; i < ref_node_naux(ref_node); i++)
    ref_node_aux(ref_node, i, clone) = ref_node_aux(ref_node, i, original);

  RSS(ref_node_metric_get(ref_node, original, m), "get original m");
  RSS(ref_node_metric_set(ref_node, clone, m), "set clone m");

  return REF_SUCCESS;
}

REF_STATUS ref_node_interpolate_edge(REF_NODE ref_node, REF_INT node0,
                                     REF_INT node1, REF_DBL node1_weight,
                                     REF_INT new_node) {
  REF_DBL log_m0[6], log_m1[6], log_m[6];
  REF_INT i;
  REF_DBL node0_weight = 1.0 - node1_weight;

  if (!ref_node_valid(ref_node, node0) || !ref_node_valid(ref_node, node1))
    RSS(REF_INVALID, "node invalid");

  for (i = 0; i < 3; i++)
    ref_node_xyz(ref_node, i, new_node) =
        node0_weight * ref_node_xyz(ref_node, i, node0) +
        node1_weight * ref_node_xyz(ref_node, i, node1);

  for (i = 0; i < ref_node_naux(ref_node); i++)
    ref_node_aux(ref_node, i, new_node) =
        node0_weight * ref_node_aux(ref_node, i, node0) +
        node1_weight * ref_node_aux(ref_node, i, node1);

  RSS(ref_node_metric_get_log(ref_node, node0, log_m0), "log 0");
  RSS(ref_node_metric_get_log(ref_node, node1, log_m1), "log 1");

  RSS(ref_matrix_weight_m(log_m0, log_m1, node1_weight, log_m), "log weight");

  RSS(ref_node_metric_set_log(ref_node, new_node, log_m), "log new");

  return REF_SUCCESS;
}

REF_STATUS ref_node_interpolate_face(REF_NODE ref_node, REF_INT node0,
                                     REF_INT node1, REF_INT node2,
                                     REF_INT new_node) {
  REF_DBL log_m0[6], log_m1[6], log_m2[6], log_m[6];
  REF_INT i;

  if (!ref_node_valid(ref_node, node0) || !ref_node_valid(ref_node, node1) ||
      !ref_node_valid(ref_node, node2))
    RSS(REF_INVALID, "node invalid");

  for (i = 0; i < 3; i++)
    ref_node_xyz(ref_node, i, new_node) =
        (1.0 / 3.0) *
        (ref_node_xyz(ref_node, i, node0) + ref_node_xyz(ref_node, i, node1) +
         ref_node_xyz(ref_node, i, node2));

  for (i = 0; i < ref_node_naux(ref_node); i++)
    ref_node_aux(ref_node, i, new_node) =
        (1.0 / 3.0) *
        (ref_node_aux(ref_node, i, node0) + ref_node_aux(ref_node, i, node1) +
         ref_node_aux(ref_node, i, node2));

  RSS(ref_node_metric_get_log(ref_node, node0, log_m0), "log 0");
  RSS(ref_node_metric_get_log(ref_node, node1, log_m1), "log 1");
  RSS(ref_node_metric_get_log(ref_node, node2, log_m2), "log 2");

  for (i = 0; i < 6; i++)
    log_m[i] = (1.0 / 3.0) * (log_m0[i] + log_m1[i] + log_m2[i]);

  RSS(ref_node_metric_set_log(ref_node, new_node, log_m), "log new");

  return REF_SUCCESS;
}

REF_STATUS ref_node_resize_aux(REF_NODE ref_node) {
  if (NULL == ref_node->aux) {
    ref_malloc(ref_node->aux, ref_node_naux(ref_node) * ref_node_max(ref_node),
               REF_DBL);
  } else {
    ref_realloc(ref_node->aux, ref_node_naux(ref_node) * ref_node_max(ref_node),
                REF_DBL);
  }
  return REF_SUCCESS;
}

REF_STATUS ref_node_bary3(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                          REF_DBL *bary) {
  REF_DBL *xyz0, *xyz1, *xyz2;
  REF_DBL total, normal[3];

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]))
    RSS(REF_INVALID, "node invalid");

  xyz0 = ref_node_xyz_ptr(ref_node, nodes[0]);
  xyz1 = ref_node_xyz_ptr(ref_node, nodes[1]);
  xyz2 = ref_node_xyz_ptr(ref_node, nodes[2]);

  RSS(ref_node_xyz_normal(xyz, xyz1, xyz2, normal), "n0");
  bary[0] = normal[1];
  RSS(ref_node_xyz_normal(xyz0, xyz, xyz2, normal), "n1");
  bary[1] = normal[1];
  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz, normal), "n2");
  bary[2] = normal[1];

  total = bary[0] + bary[1] + bary[2];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
  } else {
    REF_INT i;
    printf("%s: %d: %s: div zero total %.18e norms %.18e %.18e %.18e\n",
           __FILE__, __LINE__, __func__, total, bary[0], bary[1], bary[2]);
    for (i = 0; i < 3; i++) bary[i] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_bary3d(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                           REF_DBL *bary) {
  REF_DBL *xyz0, *xyz1, *xyz2;
  REF_DBL total, normal[3], total_normal[3];

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]))
    RSS(REF_INVALID, "node invalid");

  xyz0 = ref_node_xyz_ptr(ref_node, nodes[0]);
  xyz1 = ref_node_xyz_ptr(ref_node, nodes[1]);
  xyz2 = ref_node_xyz_ptr(ref_node, nodes[2]);

  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz2, total_normal), "n0");

  RSS(ref_node_xyz_normal(xyz, xyz1, xyz2, normal), "n0");
  bary[0] = ref_math_dot(normal, total_normal);
  RSS(ref_node_xyz_normal(xyz0, xyz, xyz2, normal), "n1");
  bary[1] = ref_math_dot(normal, total_normal);
  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz, normal), "n2");
  bary[2] = ref_math_dot(normal, total_normal);

  total = bary[0] + bary[1] + bary[2];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
  } else {
    REF_INT i;
    printf("%s: %d: %s: div zero total %.18e norms %.18e %.18e %.18e\n",
           __FILE__, __LINE__, __func__, total, bary[0], bary[1], bary[2]);
    for (i = 0; i < 3; i++) bary[i] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_bary4(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                          REF_DBL *bary) {
  REF_DBL *a, *b, *c, *d;
  REF_DBL total, m11, m12, m13;

  if (!ref_node_valid(ref_node, nodes[0]) ||
      !ref_node_valid(ref_node, nodes[1]) ||
      !ref_node_valid(ref_node, nodes[2]) ||
      !ref_node_valid(ref_node, nodes[3]))
    RSS(REF_INVALID, "node invalid");

  a = ref_node_xyz_ptr(ref_node, nodes[0]);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);
  c = ref_node_xyz_ptr(ref_node, nodes[2]);
  d = ref_node_xyz_ptr(ref_node, nodes[3]);

  a = xyz;
  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  bary[0] = (m11 - m12 + m13);
  a = ref_node_xyz_ptr(ref_node, nodes[0]);

  b = xyz;
  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  bary[1] = (m11 - m12 + m13);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);

  c = xyz;
  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  bary[2] = (m11 - m12 + m13);
  c = ref_node_xyz_ptr(ref_node, nodes[2]);

  d = xyz;
  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  bary[3] = (m11 - m12 + m13);
  d = ref_node_xyz_ptr(ref_node, nodes[3]);

  total = bary[0] + bary[1] + bary[2] + bary[3];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total) &&
      ref_math_divisible(bary[3], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
    bary[3] /= total;
  } else {
    REF_DBL volume;
    REF_INT i, smallest;
    RSS(ref_node_tet_vol(ref_node, nodes, &volume), "bary vol chk");
    printf("%s: %d: %s: div zero\ntot %.18e\nbary %.18e %.18e\n%.18e %.18e\n",
           __FILE__, __LINE__, __func__, total, bary[0], bary[1], bary[2],
           bary[3]);
    /* for walking set the smallest as the direction */
    smallest = 0;
    for (i = 1; i < 4; i++)
      if (bary[i] < bary[smallest]) smallest = i;
    for (i = 0; i < 4; i++) bary[i] = 0.0;
    bary[smallest] = -1.0;
    RSS(ref_node_tet_vol(ref_node, nodes, &volume), "bary vol chk");
    printf("vol %.18e modified bary\n%.18e %.18e\n%.18e %.18e\n", volume,
           bary[0], bary[1], bary[2], bary[3]);
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_clip_bary4(REF_DBL *orig_bary, REF_DBL *bary) {
  REF_DBL total;

  bary[0] = MAX(0.0, orig_bary[0]);
  bary[1] = MAX(0.0, orig_bary[1]);
  bary[2] = MAX(0.0, orig_bary[2]);
  bary[3] = MAX(0.0, orig_bary[3]);

  total = bary[0] + bary[1] + bary[2] + bary[3];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total) &&
      ref_math_divisible(bary[3], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
    bary[3] /= total;
  } else {
    REF_INT i, largest;
    printf("%s: %d: %s: div zero\ntot %.18e\nbary %.18e %.18e\n%.18e %.18e\n",
           __FILE__, __LINE__, __func__, total, orig_bary[0], orig_bary[1],
           orig_bary[2], orig_bary[3]);
    printf("clipped bary\n%.18e %.18e\n%.18e %.18e\n", bary[0], bary[1],
           bary[2], bary[3]);
    /* chose one node */
    largest = 0;
    for (i = 1; i < 4; i++)
      if (bary[i] > bary[largest]) largest = i;
    for (i = 0; i < 4; i++) bary[i] = 0.0;
    bary[largest] = 1.0;
    printf("modified bary\n%.18e %.18e\n%.18e %.18e\n", bary[0], bary[1],
           bary[2], bary[3]);
    return REF_DIV_ZERO;
  }

  RAS(bary[0] >= 0.0, "bary[0] not positve");
  RAS(bary[1] >= 0.0, "bary[1] not positve");
  RAS(bary[2] >= 0.0, "bary[2] not positve");
  RAS(bary[3] >= 0.0, "bary[3] not positve");

  return REF_SUCCESS;
}

REF_STATUS ref_node_tri_projection(REF_NODE ref_node, REF_INT *nodes,
                                   REF_DBL *xyz, REF_DBL *projection) {
  REF_DBL area;
  REF_DBL *a, *b, *c, *d;
  REF_DBL m11, m12, m13;
  REF_DBL vol;
  RSS(ref_node_tri_area(ref_node, nodes, &area), "area");

  a = ref_node_xyz_ptr(ref_node, nodes[0]);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);
  c = ref_node_xyz_ptr(ref_node, nodes[2]);
  d = xyz;

  m11 = (a[0] - d[0]) *
        ((b[1] - d[1]) * (c[2] - d[2]) - (c[1] - d[1]) * (b[2] - d[2]));
  m12 = (a[1] - d[1]) *
        ((b[0] - d[0]) * (c[2] - d[2]) - (c[0] - d[0]) * (b[2] - d[2]));
  m13 = (a[2] - d[2]) *
        ((b[0] - d[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - d[1]));
  vol = -(m11 - m12 + m13) / 6.0;

  if (ref_math_divisible(vol, area)) {
    *projection = 3.0 * (vol / area);
  } else {
    *projection = 0.0;
    printf("%s: %d: %s: div zero vol %.18e area %.18e\n", __FILE__, __LINE__,
           __func__, vol, area);
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_dist_to_edge(REF_NODE ref_node, REF_INT *nodes,
                                 REF_DBL *xyz, REF_DBL *distance) {
  REF_DBL *a, *b;
  REF_DBL direction[3], diff[3], normal[3];
  REF_DBL len, dis, t;

  a = ref_node_xyz_ptr(ref_node, nodes[0]);
  b = ref_node_xyz_ptr(ref_node, nodes[1]);
  direction[0] = b[0] - a[0];
  direction[1] = b[1] - a[1];
  direction[2] = b[2] - a[2];
  diff[0] = xyz[0] - a[0];
  diff[1] = xyz[1] - a[1];
  diff[2] = xyz[2] - a[2];

  dis = ref_math_dot(direction, diff);
  len = ref_math_dot(direction, direction);

  if (ref_math_divisible(dis, len)) {
    t = dis / len;
  } else {
    *distance = 0.0;
    printf("%s: %d: %s: div zero dis %.18e len %.18e\n", __FILE__, __LINE__,
           __func__, dis, len);
    return REF_DIV_ZERO;
  }

  normal[0] = xyz[0] - (t * b[0] + (1.0 - t) * a[0]);
  normal[1] = xyz[1] - (t * b[1] + (1.0 - t) * a[1]);
  normal[2] = xyz[2] - (t * b[2] + (1.0 - t) * a[2]);

  *distance = sqrt(ref_math_dot(normal, normal));

  return REF_SUCCESS;
}
REF_STATUS ref_node_dist_to_tri(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                                REF_DBL *distance) {
  REF_DBL projection;
  REF_DBL edge_dist0, edge_dist1, edge_dist2;
  REF_DBL node_dist0, node_dist1, node_dist2;
  REF_INT edge_nodes[2];
  RSS(ref_node_tri_projection(ref_node, nodes, xyz, &projection), "proj");
  projection = ABS(projection);
  edge_nodes[0] = nodes[1];
  edge_nodes[1] = nodes[2];
  RSS(ref_node_dist_to_edge(ref_node, edge_nodes, xyz, &edge_dist0), "e0");
  edge_nodes[0] = nodes[2];
  edge_nodes[1] = nodes[0];
  RSS(ref_node_dist_to_edge(ref_node, edge_nodes, xyz, &edge_dist1), "e1");
  edge_nodes[0] = nodes[0];
  edge_nodes[1] = nodes[1];
  RSS(ref_node_dist_to_edge(ref_node, edge_nodes, xyz, &edge_dist2), "e2");
  node_dist0 = pow(xyz[0] - ref_node_xyz(ref_node, 0, nodes[0]), 2) +
               pow(xyz[1] - ref_node_xyz(ref_node, 1, nodes[0]), 2) +
               pow(xyz[2] - ref_node_xyz(ref_node, 2, nodes[0]), 2);
  node_dist1 = pow(xyz[0] - ref_node_xyz(ref_node, 0, nodes[1]), 2) +
               pow(xyz[1] - ref_node_xyz(ref_node, 1, nodes[1]), 2) +
               pow(xyz[2] - ref_node_xyz(ref_node, 2, nodes[1]), 2);
  node_dist2 = pow(xyz[0] - ref_node_xyz(ref_node, 0, nodes[2]), 2) +
               pow(xyz[1] - ref_node_xyz(ref_node, 1, nodes[2]), 2) +
               pow(xyz[2] - ref_node_xyz(ref_node, 2, nodes[2]), 2);
  *distance = MIN(MIN(MIN(projection, edge_dist0), MIN(edge_dist1, edge_dist2)),
                  MIN(node_dist0, MIN(node_dist1, node_dist2)));

  return REF_SUCCESS;
}

REF_STATUS ref_node_xyz_grad(REF_DBL *xyzs[4], REF_DBL *scalar,
                             REF_DBL *gradient) {
  REF_DBL vol, norm1[3], norm2[3], norm3[3];
  REF_DBL *xyz0, *xyz1, *xyz2;

  gradient[0] = 0.0;
  gradient[1] = 0.0;
  gradient[2] = 0.0;

  RSS(ref_node_xyz_vol(xyzs, &vol), "vol");
  vol *= -6.0;

  xyz0 = xyzs[0];
  xyz1 = xyzs[3];
  xyz2 = xyzs[2];
  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz2, norm1), "vol");

  xyz0 = xyzs[0];
  xyz1 = xyzs[1];
  xyz2 = xyzs[3];
  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz2, norm2), "vol");

  xyz0 = xyzs[0];
  xyz1 = xyzs[2];
  xyz2 = xyzs[1];
  RSS(ref_node_xyz_normal(xyz0, xyz1, xyz2, norm3), "vol");

  gradient[0] = (scalar[1] - scalar[0]) * norm1[0] +
                (scalar[2] - scalar[0]) * norm2[0] +
                (scalar[3] - scalar[0]) * norm3[0];
  gradient[1] = (scalar[1] - scalar[0]) * norm1[1] +
                (scalar[2] - scalar[0]) * norm2[1] +
                (scalar[3] - scalar[0]) * norm3[1];
  gradient[2] = (scalar[1] - scalar[0]) * norm1[2] +
                (scalar[2] - scalar[0]) * norm2[2] +
                (scalar[3] - scalar[0]) * norm3[2];

  if (ref_math_divisible(gradient[0], vol) &&
      ref_math_divisible(gradient[1], vol) &&
      ref_math_divisible(gradient[2], vol)) {
    gradient[0] /= vol;
    gradient[1] /= vol;
    gradient[2] /= vol;
  } else {
    gradient[0] = 0.0;
    gradient[1] = 0.0;
    gradient[2] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tet_grad_nodes(REF_NODE ref_node, REF_INT *nodes,
                                   REF_DBL *scalar, REF_DBL *gradient) {
  REF_DBL vol, norm1[3], norm2[3], norm3[3];
  REF_INT face[3];

  gradient[0] = 0.0;
  gradient[1] = 0.0;
  gradient[2] = 0.0;

  RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
  vol *= -6.0;

  face[0] = nodes[0];
  face[1] = nodes[3];
  face[2] = nodes[2];
  RSS(ref_node_tri_normal(ref_node, face, norm1), "vol");

  face[0] = nodes[0];
  face[1] = nodes[1];
  face[2] = nodes[3];
  RSS(ref_node_tri_normal(ref_node, face, norm2), "vol");

  face[0] = nodes[0];
  face[1] = nodes[2];
  face[2] = nodes[1];
  RSS(ref_node_tri_normal(ref_node, face, norm3), "vol");

  gradient[0] = (scalar[nodes[1]] - scalar[nodes[0]]) * norm1[0] +
                (scalar[nodes[2]] - scalar[nodes[0]]) * norm2[0] +
                (scalar[nodes[3]] - scalar[nodes[0]]) * norm3[0];
  gradient[1] = (scalar[nodes[1]] - scalar[nodes[0]]) * norm1[1] +
                (scalar[nodes[2]] - scalar[nodes[0]]) * norm2[1] +
                (scalar[nodes[3]] - scalar[nodes[0]]) * norm3[1];
  gradient[2] = (scalar[nodes[1]] - scalar[nodes[0]]) * norm1[2] +
                (scalar[nodes[2]] - scalar[nodes[0]]) * norm2[2] +
                (scalar[nodes[3]] - scalar[nodes[0]]) * norm3[2];

  if (ref_math_divisible(gradient[0], vol) &&
      ref_math_divisible(gradient[1], vol) &&
      ref_math_divisible(gradient[2], vol)) {
    gradient[0] /= vol;
    gradient[1] /= vol;
    gradient[2] /= vol;
  } else {
    gradient[0] = 0.0;
    gradient[1] = 0.0;
    gradient[2] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_nearest_xyz(REF_NODE ref_node, REF_DBL *xyz,
                                REF_INT *closest_node, REF_DBL *distance) {
  REF_INT node;
  REF_DBL dist;
  *closest_node = REF_EMPTY;
  *distance = 1.0e100;
  each_ref_node_valid_node(ref_node, node) {
    dist = pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
           pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
           pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2);
    if (dist < *distance) {
      *closest_node = node;
      *distance = dist;
    }
  }
  *distance = sqrt(*distance);
  return REF_SUCCESS;
}

REF_STATUS ref_node_selection(REF_NODE ref_node, REF_DBL *elements,
                              REF_INT position, REF_DBL *value) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_DBL *pack, *sorted;
  REF_INT *order;
  REF_DBL *all_elements;
  REF_INT *all_source, *all_order;
  REF_INT node, nnode, all_nnode;

  ref_malloc(pack, ref_node_n(ref_node), REF_DBL);
  ref_malloc(sorted, ref_node_n(ref_node), REF_DBL);
  ref_malloc(order, ref_node_n(ref_node), REF_INT);
  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      pack[nnode] = elements[node];
      nnode++;
    }
  }
  RSS(ref_sort_heap_dbl(nnode, pack, order), "heap");
  for (node = 0; node < nnode; node++) {
    sorted[node] = pack[order[node]];
  }

  RSS(ref_mpi_allconcat(ref_mpi, 1, nnode, (void *)(sorted), &all_nnode,
                        &all_source, (void **)(&all_elements), REF_DBL_TYPE),
      "concat");
  ref_malloc(all_order, all_nnode, REF_INT);

  RSS(ref_sort_heap_dbl(all_nnode, all_elements, all_order), "heap");

  *value = all_elements[all_order[position]];

  ref_free(all_order);
  ref_free(all_source);
  ref_free(all_elements);
  ref_free(order);
  ref_free(sorted);
  ref_free(pack);

  return REF_SUCCESS;
}

REF_STATUS ref_node_push_unused(REF_NODE ref_node, REF_GLOB unused_global) {
  if (ref_node_max_unused(ref_node) == ref_node_n_unused(ref_node)) {
    ref_node_max_unused(ref_node) += 1000;
    ref_realloc(ref_node->unused_global, ref_node_max_unused(ref_node),
                REF_GLOB);
  }

  ref_node->unused_global[ref_node_n_unused(ref_node)] = unused_global;

  ref_node_n_unused(ref_node)++;

  return REF_SUCCESS;
}

REF_STATUS ref_node_pop_unused(REF_NODE ref_node, REF_GLOB *new_global) {
  if (0 == ref_node_n_unused(ref_node)) {
    *new_global = REF_EMPTY;
    return REF_FAILURE;
  }

  ref_node_n_unused(ref_node)--;
  *new_global = ref_node->unused_global[ref_node_n_unused(ref_node)];

  return REF_SUCCESS;
}

REF_STATUS ref_node_shift_unused(REF_NODE ref_node, REF_GLOB equal_and_above,
                                 REF_GLOB shift) {
  REF_INT i;

  for (i = 0; i < ref_node_n_unused(ref_node); i++) {
    if (ref_node->unused_global[i] >= equal_and_above) {
      ref_node->unused_global[i] += shift;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_node_sort_unused(REF_NODE ref_node) {
  REF_INT *order;
  REF_GLOB *new_order;
  REF_INT i;

  /* see if it is too short to require sorting */
  if (2 > ref_node_n_unused(ref_node)) return REF_SUCCESS;

  ref_malloc(order, ref_node_n_unused(ref_node), REF_INT);
  ref_malloc(new_order, ref_node_n_unused(ref_node), REF_GLOB);

  RSS(ref_sort_heap_glob(ref_node_n_unused(ref_node), ref_node->unused_global,
                         order),
      "heap");

  for (i = 0; i < ref_node_n_unused(ref_node); i++) {
    new_order[i] = ref_node->unused_global[order[i]];
  }

  for (i = 0; i < ref_node_n_unused(ref_node); i++) {
    ref_node->unused_global[i] = new_order[i];
  }

  ref_free(new_order);
  ref_free(order);

  return REF_SUCCESS;
}

REF_STATUS ref_node_erase_unused(REF_NODE ref_node) {
  ref_node_n_unused(ref_node) = 0;
  return REF_SUCCESS;
}

REF_STATUS ref_node_allgather_unused(REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT i;
  REF_GLOB *local_copy;
  REF_INT proc;
  REF_INT *counts;
  REF_INT total_count;

  ref_malloc(counts, ref_mpi_n(ref_mpi), REF_INT);

  RSS(ref_mpi_allgather(ref_mpi, &(ref_node_n_unused(ref_node)), counts,
                        REF_INT_TYPE),
      "gather size");

  total_count = 0;
  each_ref_mpi_part(ref_mpi, proc) total_count += counts[proc];

  ref_malloc(local_copy, ref_node_n_unused(ref_node), REF_GLOB);
  for (i = 0; i < ref_node_n_unused(ref_node); i++) {
    local_copy[i] = ref_node->unused_global[i];
  }

  if (total_count > ref_node_max_unused(ref_node)) {
    ref_node_max_unused(ref_node) = total_count;
    ref_free(ref_node->unused_global);
    ref_malloc(ref_node->unused_global, ref_node_max_unused(ref_node),
               REF_GLOB);
  }

  RSS(ref_mpi_allgatherv(ref_mpi, local_copy, counts, ref_node->unused_global,
                         REF_GLOB_TYPE),
      "gather values");

  ref_node_n_unused(ref_node) = total_count;

  ref_free(local_copy);
  ref_free(counts);

  return REF_SUCCESS;
}
