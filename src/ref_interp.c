
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

#include "ref_interp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_search.h"

#define MAX_NODE_LIST (200)

#define ref_interp_bary_inside(ref_interp, bary)                             \
  ((bary)[0] >= (ref_interp)->inside && (bary)[1] >= (ref_interp)->inside && \
   (bary)[2] >= (ref_interp)->inside && (bary)[3] >= (ref_interp)->inside)

static REF_STATUS ref_interp_exhaustive_tet_around_node(REF_GRID ref_grid,
                                                        REF_INT node,
                                                        REF_DBL *xyz,
                                                        REF_INT *cell,
                                                        REF_DBL *bary) {
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, candidate, best_candidate;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
  REF_STATUS status;

  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_having_node(ref_cell, node, item, candidate) {
    RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
    status = ref_node_bary4(ref_node, nodes, xyz, current_bary);
    RXS(status, REF_DIV_ZERO, "bary");
    if (REF_SUCCESS == status) { /* exclude REF_DIV_ZERO */
      min_bary = MIN(MIN(current_bary[0], current_bary[1]),
                     MIN(current_bary[2], current_bary[3]));
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = candidate;
        best_bary = min_bary;
      }
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find cell");

  *cell = best_candidate;
  RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
  RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_exhaustive_tri_around_node(REF_GRID ref_grid,
                                                        REF_INT node,
                                                        REF_DBL *xyz,
                                                        REF_INT *cell,
                                                        REF_DBL *bary) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, candidate, best_candidate;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
  REF_STATUS status;

  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_having_node(ref_cell, node, item, candidate) {
    RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
    status = ref_node_bary3(ref_node, nodes, xyz, current_bary);
    current_bary[3] = 0.0;
    RXS(status, REF_DIV_ZERO, "bary");
    if (REF_SUCCESS == status) { /* exclude REF_DIV_ZERO */
      min_bary = MIN(MIN(current_bary[0], current_bary[1]),
                     MIN(current_bary[2], current_bary[3]));
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = candidate;
        best_bary = min_bary;
      }
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find cell");

  *cell = best_candidate;
  RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
  RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
  bary[3] = 0.0;

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_bounding_sphere4(REF_NODE ref_node, REF_INT *nodes,
                                              REF_DBL *center,
                                              REF_DBL *radius) {
  REF_INT i;
  for (i = 0; i < 3; i++)
    center[i] = 0.25 * (ref_node_xyz(ref_node, i, nodes[0]) +
                        ref_node_xyz(ref_node, i, nodes[1]) +
                        ref_node_xyz(ref_node, i, nodes[2]) +
                        ref_node_xyz(ref_node, i, nodes[3]));
  *radius = 0.0;
  for (i = 0; i < 4; i++)
    *radius = MAX(
        *radius, sqrt(pow(ref_node_xyz(ref_node, 0, nodes[i]) - center[0], 2) +
                      pow(ref_node_xyz(ref_node, 1, nodes[i]) - center[1], 2) +
                      pow(ref_node_xyz(ref_node, 2, nodes[i]) - center[2], 2)));
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_bounding_sphere3(REF_NODE ref_node, REF_INT *nodes,
                                              REF_DBL *center,
                                              REF_DBL *radius) {
  REF_INT i;
  for (i = 0; i < 3; i++)
    center[i] = (1.0 / 3.0) * (ref_node_xyz(ref_node, i, nodes[0]) +
                               ref_node_xyz(ref_node, i, nodes[1]) +
                               ref_node_xyz(ref_node, i, nodes[2]));
  *radius = 0.0;
  for (i = 0; i < 3; i++)
    *radius = MAX(
        *radius, sqrt(pow(ref_node_xyz(ref_node, 0, nodes[i]) - center[0], 2) +
                      pow(ref_node_xyz(ref_node, 1, nodes[i]) - center[1], 2) +
                      pow(ref_node_xyz(ref_node, 2, nodes[i]) - center[2], 2)));
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_create_search(REF_INTERP ref_interp) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tet = ref_grid_tet(from_grid);
  REF_CELL from_tri = ref_grid_tri(from_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius;
  REF_SEARCH ref_search;

  if (ref_grid_twod(from_grid)) {
    RSS(ref_search_create(&ref_search, ref_cell_n(from_tri)), "create search");
    each_ref_cell_valid_cell_with_nodes(from_tri, cell, nodes) {
      RSS(ref_interp_bounding_sphere3(from_node, nodes, center, &radius), "b");
      RSS(ref_search_insert(ref_search, cell, center,
                            ref_interp_search_donor_scale(ref_interp) * radius),
          "ins");
    }
  } else {
    RSS(ref_search_create(&ref_search, ref_cell_n(from_tet)), "create search");
    each_ref_cell_valid_cell_with_nodes(from_tet, cell, nodes) {
      RSS(ref_interp_bounding_sphere4(from_node, nodes, center, &radius), "b");
      RSS(ref_search_insert(ref_search, cell, center,
                            ref_interp_search_donor_scale(ref_interp) * radius),
          "ins");
    }
  }
  ref_interp_search(ref_interp) = ref_search;

  return REF_SUCCESS;
}

REF_STATUS ref_interp_create(REF_INTERP *ref_interp_ptr, REF_GRID from_grid,
                             REF_GRID to_grid) {
  REF_INTERP ref_interp;
  REF_INT max = ref_node_max(ref_grid_node(to_grid));

  ref_malloc(*ref_interp_ptr, 1, REF_INTERP_STRUCT);
  ref_interp = (*ref_interp_ptr);

  ref_interp_from_grid(ref_interp) = from_grid;
  ref_interp_to_grid(ref_interp) = to_grid;

  ref_interp_mpi(ref_interp) = ref_grid_mpi(ref_interp_from_grid(ref_interp));

  ref_interp->instrument = REF_FALSE;
  ref_interp->continuously = REF_FALSE;
  ref_interp->n_walk = 0;
  ref_interp->n_terminated = 0;
  ref_interp->walk_steps = 0;
  ref_interp->n_geom = 0;
  ref_interp->n_geom_fail = 0;
  ref_interp->n_tree = 0;
  ref_interp->tree_cells = 0;
  ref_interp_max(ref_interp) = max;
  ref_malloc_init(ref_interp->agent_hired, max, REF_BOOL, REF_FALSE);
  ref_malloc_init(ref_interp->cell, max, REF_INT, REF_EMPTY);
  ref_malloc_init(ref_interp->part, max, REF_INT, REF_EMPTY);
  ref_malloc(ref_interp->bary, 4 * max, REF_DBL);
  ref_interp->inside = -1.0e-12; /* inside tolerence */
  ref_interp->bound = -0.1;      /* bound tolerence */

  RSS(ref_agents_create(&(ref_interp->ref_agents), ref_interp_mpi(ref_interp)),
      "add agents");
  RSS(ref_list_create(&(ref_interp->visualize)), "add list");
  ref_interp_search_fuzz(ref_interp) = 1.0e-12;
  ref_interp_search_donor_scale(ref_interp) = 2.0;
  RSS(ref_interp_create_search(ref_interp), "fill search");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_resize(REF_INTERP ref_interp, REF_INT max) {
  REF_INT old = ref_interp_max(ref_interp);

  ref_realloc_init(ref_interp->agent_hired, old, max, REF_BOOL, REF_FALSE);
  ref_realloc_init(ref_interp->cell, old, max, REF_INT, REF_EMPTY);
  ref_realloc_init(ref_interp->part, old, max, REF_INT, REF_EMPTY);
  ref_realloc(ref_interp->bary, 4 * max, REF_DBL);

  ref_interp_max(ref_interp) = max;

  return REF_SUCCESS;
}

REF_STATUS ref_interp_reset(REF_INTERP ref_interp) {
  REF_INT node_max =
      ref_node_max(ref_grid_node(ref_interp_to_grid(ref_interp)));
  REF_INT node;
  for (node = 0; node < ref_interp_max(ref_interp); node++) {
    RAS(!ref_interp->agent_hired[node],
        "agent should not be hired during reset");
  }
  if (node_max > ref_interp_max(ref_interp)) {
    RSS(ref_interp_resize(ref_interp, node_max), "protective resize");
  }
  for (node = 0; node < ref_interp_max(ref_interp); node++) {
    ref_interp->cell[node] = REF_EMPTY;
  }
  for (node = 0; node < ref_interp_max(ref_interp); node++) {
    ref_interp->part[node] = REF_EMPTY;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_create_identity(REF_INTERP *ref_interp_ptr,
                                      REF_GRID to_grid) {
  REF_INTERP ref_interp;
  REF_GRID from_grid;
  REF_NODE to_node;
  REF_INT node;
  REF_DBL max_error;

  RSS(ref_grid_deep_copy(&from_grid, to_grid), "import");
  RSS(ref_interp_create(ref_interp_ptr, from_grid, to_grid), "create");

  ref_interp = *ref_interp_ptr;
  to_node = ref_grid_node(to_grid);

  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      REIS(REF_EMPTY, ref_interp->cell[node], "identity already found?");
      if (ref_grid_twod(to_grid)) {
        if (ref_cell_node_empty(ref_grid_tri(to_grid), node)) continue;
        RSS(ref_interp_exhaustive_tri_around_node(
                from_grid, node, ref_node_xyz_ptr(to_node, node),
                &(ref_interp->cell[node]), &(ref_interp->bary[4 * node])),
            "tri around node");
      } else {
        RSS(ref_interp_exhaustive_tet_around_node(
                from_grid, node, ref_node_xyz_ptr(to_node, node),
                &(ref_interp->cell[node]), &(ref_interp->bary[4 * node])),
            "tet around node");
      }
      ref_interp->part[node] = ref_mpi_rank(ref_grid_mpi(from_grid));
      if (!ref_interp_bary_inside(ref_interp, &(ref_interp->bary[4 * node]))) {
        printf("info bary %e %e %e %e identity tol\n",
               ref_interp->bary[0 + 4 * node], ref_interp->bary[1 + 4 * node],
               ref_interp->bary[2 + 4 * node], ref_interp->bary[3 + 4 * node]);
      }
    }
  }

  RSS(ref_interp_max_error(ref_interp, &max_error), "max error");
  if (ref_mpi_once(ref_grid_mpi(to_grid)) && max_error > 1.0e-12) {
    printf("warning %e max error for identity background grid\n", max_error);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_free(REF_INTERP ref_interp) {
  if (NULL == (void *)ref_interp) return REF_NULL;
  ref_search_free(ref_interp->ref_search);
  ref_list_free(ref_interp->visualize);
  ref_agents_free(ref_interp->ref_agents);
  ref_free(ref_interp->bary);
  ref_free(ref_interp->part);
  ref_free(ref_interp->cell);
  ref_free(ref_interp->agent_hired);
  ref_free(ref_interp);
  return REF_SUCCESS;
}

REF_STATUS ref_interp_pack(REF_INTERP ref_interp, REF_INT *n2o) {
  REF_INT *int_copy;
  REF_DBL *dbl_copy;
  REF_INT i, node, n, max;

  if (NULL == ref_interp) return REF_SUCCESS;

  n = ref_node_n(ref_grid_node(ref_interp_to_grid(ref_interp)));
  max = ref_interp_max(ref_interp);
  if (n > max) {
    RSS(ref_interp_resize(ref_interp, max), "match node max");
  }
  REIS(0, ref_agents_n(ref_interp->ref_agents), "can't pack active agents");

  for (node = 0; node < max; node++) {
    REIS(REF_FALSE, ref_interp->agent_hired[node], "can't pack hired agents");
  }

  ref_malloc(int_copy, max, REF_INT);
  for (node = 0; node < max; node++) {
    int_copy[node] = ref_interp->cell[node];
  }
  for (node = 0; node < n; node++) {
    ref_interp->cell[node] = int_copy[n2o[node]];
  }
  for (node = n; node < max; node++) {
    ref_interp->cell[node] = REF_EMPTY;
  }
  for (node = 0; node < max; node++) {
    int_copy[node] = ref_interp->part[node];
  }
  for (node = 0; node < n; node++) {
    ref_interp->part[node] = int_copy[n2o[node]];
  }
  for (node = n; node < max; node++) {
    ref_interp->part[node] = REF_EMPTY;
  }
  ref_free(int_copy);

  ref_malloc(dbl_copy, 4 * max, REF_DBL);
  for (node = 0; node < max; node++) {
    for (i = 0; i < 4; i++) {
      dbl_copy[i + 4 * node] = ref_interp->bary[i + 4 * node];
    }
  }
  for (node = 0; node < n; node++) {
    for (i = 0; i < 4; i++) {
      ref_interp->bary[i + 4 * node] = dbl_copy[i + 4 * n2o[node]];
    }
  }
  ref_free(dbl_copy);

  return REF_SUCCESS;
}

REF_STATUS ref_interp_remove(REF_INTERP ref_interp, REF_INT node) {
  if (NULL == ref_interp) return REF_SUCCESS;
  if (!ref_interp_continuously(ref_interp)) return REF_SUCCESS;
  REIS(REF_FALSE, ref_interp->agent_hired[node], "remove node agent hired");
  RUS(REF_EMPTY, ref_interp_cell(ref_interp, node), "remove node no located");
  ref_interp->cell[node] = REF_EMPTY;
  return REF_SUCCESS;
}

REF_STATUS ref_interp_exhaustive_enclosing_tet(REF_GRID ref_grid, REF_DBL *xyz,
                                               REF_INT *cell, REF_DBL *bary) {
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT candidate, best_candidate;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
  REF_STATUS status;

  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_valid_cell(ref_cell, candidate) {
    RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
    status = ref_node_bary4(ref_node, nodes, xyz, current_bary);
    RXS(status, REF_DIV_ZERO, "bary");
    if (REF_SUCCESS == status) { /* exclude REF_DIV_ZERO */
      min_bary = MIN(MIN(current_bary[0], current_bary[1]),
                     MIN(current_bary[2], current_bary[3]));
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = candidate;
        best_bary = min_bary;
      }
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find cell");

  *cell = best_candidate;
  RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
  RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_enclosing_tet_in_list(REF_GRID ref_grid,
                                            REF_LIST ref_list, REF_DBL *xyz,
                                            REF_INT *cell, REF_DBL *bary) {
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, candidate, best_candidate;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
  REF_STATUS status;

  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_list_item(ref_list, item) {
    candidate = ref_list_value(ref_list, item);
    RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
    status = ref_node_bary4(ref_node, nodes, xyz, current_bary);
    RXS(status, REF_DIV_ZERO, "bary");
    if (REF_SUCCESS == status) { /* exclude REF_DIV_ZERO */

      min_bary = MIN(MIN(current_bary[0], current_bary[1]),
                     MIN(current_bary[2], current_bary[3]));
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = candidate;
        best_bary = min_bary;
      }
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find cell");

  *cell = best_candidate;
  RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
  RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_best_tri_in_list(REF_GRID ref_grid, REF_LIST ref_list,
                                       REF_DBL *xyz, REF_INT *cell,
                                       REF_DBL *bary) {
  REF_CELL ref_tri = ref_grid_tri(ref_grid);
  REF_CELL ref_tet = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, tri_candidate, best_candidate, tet0, tet1, i, j;
  REF_DBL current_bary3[3];
  REF_DBL best_bary, min_bary;
  REF_STATUS status;

  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_list_item(ref_list, item) {
    tri_candidate = ref_list_value(ref_list, item);
    RSS(ref_cell_nodes(ref_tri, tri_candidate, tri_nodes), "cell");
    status = ref_node_bary3d(ref_node, tri_nodes, xyz, current_bary3);
    RXS(status, REF_DIV_ZERO, "bary");
    if (REF_SUCCESS == status) { /* exclude REF_DIV_ZERO */

      min_bary = MIN(MIN(current_bary3[0], current_bary3[1]), current_bary3[2]);
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = tri_candidate;
        best_bary = min_bary;
      }
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find tri");

  RSS(ref_cell_nodes(ref_tri, best_candidate, tri_nodes), "cell");
  RSS(ref_node_bary3d(ref_node, tri_nodes, xyz, current_bary3), "bary3d");
  tri_nodes[3] = tri_nodes[0];
  RSS(ref_cell_with_face(ref_tet, tri_nodes, &tet0, &tet1), "find with face");
  RUS(REF_EMPTY, tet0, "no tet for tri");
  REIS(REF_EMPTY, tet1, "found two tets for boundary face");

  *cell = tet0;
  RSS(ref_cell_nodes(ref_tet, tet0, tet_nodes), "cell");
  for (i = 0; i < 4; i++) {
    bary[i] = 0.0;
    for (j = 0; j < 3; j++) {
      if (tri_nodes[j] == tet_nodes[i]) bary[i] = current_bary3[j];
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_update_agent_tet_seed(REF_INTERP ref_interp, REF_INT id,
                                            REF_INT node0, REF_INT node1,
                                            REF_INT node2) {
  REF_GRID ref_grid = ref_interp_from_grid(ref_interp);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tets = ref_grid_tet(ref_grid);
  REF_CELL tris = ref_grid_tri(ref_grid);
  REF_AGENTS ref_agents = ref_interp->ref_agents;
  REF_INT face_nodes[4], cell0, cell1;
  REF_INT tri, node;

  face_nodes[0] = node0;
  face_nodes[1] = node1;
  face_nodes[2] = node2;
  face_nodes[3] = node0;

  RSS(ref_cell_with_face(tets, face_nodes, &cell0, &cell1), "next");
  if (REF_EMPTY == cell0) THROW("bary update missing first");
  if (REF_EMPTY == cell1) {
    /* if it is off proc */
    if (!ref_node_owned(ref_node, node0) && !ref_node_owned(ref_node, node1) &&
        !ref_node_owned(ref_node, node2)) {
      /* pick at pseudo random */
      node = face_nodes[rand() % 3];
      ref_agent_part(ref_agents, id) = ref_node_part(ref_node, node);
      ref_agent_seed(ref_agents, id) = REF_EMPTY;
      ref_agent_global(ref_agents, id) = ref_node_global(ref_node, node);
      ref_agent_mode(ref_agents, id) = REF_AGENT_HOP_PART;
      return REF_SUCCESS;
    }
    /* hit boundary, but verifying */
    RSS(ref_cell_with(tris, face_nodes, &tri), "boundary tri expected");
    ref_agent_mode(ref_agents, id) = REF_AGENT_AT_BOUNDARY;
    return REF_SUCCESS;
  }

  if (ref_agent_seed(ref_agents, id) == cell0) {
    ref_agent_seed(ref_agents, id) = cell1;
    return REF_SUCCESS;
  }
  if (ref_agent_seed(ref_agents, id) == cell1) {
    ref_agent_seed(ref_agents, id) = cell0;
    return REF_SUCCESS;
  }

  return REF_NOT_FOUND;
}

static REF_STATUS ref_update_agent_tri_seed(REF_INTERP ref_interp, REF_INT id,
                                            REF_INT node0, REF_INT node1) {
  REF_GRID ref_grid = ref_interp_from_grid(ref_interp);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tris = ref_grid_tri(ref_grid);
  REF_AGENTS ref_agents = ref_interp->ref_agents;
  REF_INT ncell, cells[2];
  REF_INT node;

  RSS(ref_cell_list_with2(tris, node0, node1, 2, &ncell, cells),
      "more then two");
  if (0 == ncell) THROW("bary update missing first");
  if (1 == ncell) {
    /* if it is off proc */
    if (!ref_node_owned(ref_node, node0) && !ref_node_owned(ref_node, node1)) {
      /* pick at pseudo random */
      node = node0;
      if (1 == rand() % 2) node = node1;
      ref_agent_part(ref_agents, id) = ref_node_part(ref_node, node);
      ref_agent_seed(ref_agents, id) = REF_EMPTY;
      ref_agent_global(ref_agents, id) = ref_node_global(ref_node, node);
      ref_agent_mode(ref_agents, id) = REF_AGENT_HOP_PART;
      return REF_SUCCESS;
    }
    /* hit boundary, skip verify for no geom */
    ref_agent_mode(ref_agents, id) = REF_AGENT_AT_BOUNDARY;
    return REF_SUCCESS;
  }

  if (ref_agent_seed(ref_agents, id) == cells[0]) {
    ref_agent_seed(ref_agents, id) = cells[1];
    return REF_SUCCESS;
  }
  if (ref_agent_seed(ref_agents, id) == cells[1]) {
    ref_agent_seed(ref_agents, id) = cells[0];
    return REF_SUCCESS;
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_interp_tattle(REF_INTERP ref_interp, REF_INT node) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_INT i, j, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[3], error;
  if (NULL == ref_interp) {
    printf("NULL interp %d\n", node);
    return REF_SUCCESS;
  }
  if (node >= ref_interp_max(ref_interp)) {
    printf("max %d too samll for %d\n", ref_interp_max(ref_interp), node);
    return REF_SUCCESS;
  }
  if (REF_EMPTY == ref_interp->cell[node]) {
    printf("empty cell %d\n", node);
    return REF_SUCCESS;
  }
  printf("cell %d part %d node %d\n", ref_interp->cell[node],
         ref_interp->part[node], node);
  printf("bary %f %f %f %f\n", ref_interp->bary[0 + 4 * node],
         ref_interp->bary[1 + 4 * node], ref_interp->bary[2 + 4 * node],
         ref_interp->bary[3 + 4 * node]);
  printf("target %f %f %f\n",
         ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)), 0, node),
         ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)), 1, node),
         ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)), 2, node));
  if (ref_mpi_rank(ref_mpi) == ref_interp->part[node]) {
    RSS(ref_cell_nodes(ref_grid_tet(ref_interp_from_grid(ref_interp)),
                       ref_interp->cell[node], nodes),
        "node needs to be localized");

    for (i = 0; i < 3; i++) xyz[i] = 0.0;
    for (j = 0; j < 4; j++) {
      for (i = 0; i < 3; i++) {
        xyz[i] += ref_interp->bary[j + 4 * node] *
                  ref_node_xyz(ref_grid_node(ref_interp_from_grid(ref_interp)),
                               i, nodes[j]);
      }
    }

    error =
        pow(xyz[0] - ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)),
                                  0, node),
            2) +
        pow(xyz[1] - ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)),
                                  1, node),
            2) +
        pow(xyz[2] - ref_node_xyz(ref_grid_node(ref_interp_to_grid(ref_interp)),
                                  2, node),
            2);
    error = sqrt(error);
    printf("interp %f %f %f error %e\n", xyz[0], xyz[1], xyz[2], error);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_walk_agent(REF_INTERP ref_interp, REF_INT id) {
  REF_GRID ref_grid = ref_interp_from_grid(ref_interp);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT i, limit;
  REF_DBL bary[4];
  REF_AGENTS ref_agents = ref_interp->ref_agents;

  if (ref_grid_twod(ref_grid)) ref_cell = ref_grid_tri(ref_grid);

  limit = 215; /* 10e6^(1/3), required 108 for twod testcase  */

  each_ref_agent_step(ref_agents, id, limit) {
    /* return if no longer walking */
    if (REF_AGENT_WALKING != ref_agent_mode(ref_agents, id)) {
      return REF_SUCCESS;
    }

    RSB(ref_cell_nodes(ref_cell, ref_agent_seed(ref_agents, id), nodes), "cell",
        { ref_agents_tattle(ref_agents, id, "cell_nodes in walk"); });
    /* when REF_DIV_ZERO, min bary is preserved */
    if (ref_grid_twod(ref_grid)) {
      bary[3] = 0.0;
      RXS(ref_node_bary3(ref_node, nodes, ref_agent_xyz_ptr(ref_agents, id),
                         bary),
          REF_DIV_ZERO, "bary");
    } else {
      RXS(ref_node_bary4(ref_node, nodes, ref_agent_xyz_ptr(ref_agents, id),
                         bary),
          REF_DIV_ZERO, "bary");
    }

    if (ref_agent_step(ref_agents, id) > (limit)) {
      printf("bary %e %e %e %e inside %e\n", bary[0], bary[1], bary[2], bary[3],
             ref_interp->inside);
      RSS(ref_agents_tattle(ref_agents, id, "many steps"), "tat");
    }

    if (ref_interp_bary_inside(ref_interp, bary)) {
      ref_agent_mode(ref_agents, id) = REF_AGENT_ENCLOSING;
      for (i = 0; i < 4; i++) ref_agent_bary(ref_agents, i, id) = bary[i];
      return REF_SUCCESS;
    }

    if (ref_grid_twod(ref_grid)) {
      /* less than */
      if (bary[0] < bary[1] && bary[0] < bary[2]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[1], nodes[2]),
            "1 2");
        continue;
      }
      if (bary[1] < bary[0] && bary[1] < bary[2]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[2], nodes[0]),
            "2 0");
        continue;
      }
      if (bary[2] < bary[0] && bary[2] < bary[0]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[0], nodes[1]),
            "0 1");
        continue;
      }
      /* less than or equal */
      if (bary[0] <= bary[1] && bary[0] <= bary[2]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[1], nodes[2]),
            "1 2");
        continue;
      }
      if (bary[1] <= bary[0] && bary[1] <= bary[2]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[2], nodes[0]),
            "2 0");
        continue;
      }
      if (bary[2] <= bary[0] && bary[2] <= bary[0]) {
        RSS(ref_update_agent_tri_seed(ref_interp, id, nodes[0], nodes[1]),
            "0 1");
        continue;
      }
    } else {
      /* less than */
      if (bary[0] < bary[1] && bary[0] < bary[2] && bary[0] < bary[3]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[1], nodes[2],
                                      nodes[3]),
            "1 2 3");
        continue;
      }

      if (bary[1] < bary[0] && bary[1] < bary[3] && bary[1] < bary[2]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[3],
                                      nodes[2]),
            "0 3 2");
        continue;
      }

      if (bary[2] < bary[0] && bary[2] < bary[1] && bary[2] < bary[3]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[1],
                                      nodes[3]),
            "0 1 3");
        continue;
      }

      if (bary[3] < bary[0] && bary[3] < bary[2] && bary[3] < bary[1]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[2],
                                      nodes[1]),
            "0 2 1");
        continue;
      }

      /* less than or equal */
      if (bary[0] <= bary[1] && bary[0] <= bary[2] && bary[0] <= bary[3]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[1], nodes[2],
                                      nodes[3]),
            "1 2 3");
        continue;
      }

      if (bary[1] <= bary[0] && bary[1] <= bary[3] && bary[1] <= bary[2]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[3],
                                      nodes[2]),
            "0 3 2");
        continue;
      }

      if (bary[2] <= bary[0] && bary[2] <= bary[1] && bary[2] <= bary[3]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[1],
                                      nodes[3]),
            "0 1 3");
        continue;
      }

      if (bary[3] <= bary[0] && bary[3] <= bary[2] && bary[3] <= bary[1]) {
        RSS(ref_update_agent_tet_seed(ref_interp, id, nodes[0], nodes[2],
                                      nodes[1]),
            "0 2 1");
        continue;
      }
    }
    THROW("unable to find the next step");
  }

  /* steps reached limit */
  ref_agent_mode(ref_agents, id) = REF_AGENT_TERMINATED;

  return REF_SUCCESS;
}

REF_STATUS ref_interp_push_onto_queue(REF_INTERP ref_interp, REF_INT node) {
  REF_GRID ref_grid = ref_interp_to_grid(ref_interp);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_AGENTS ref_agents = ref_interp->ref_agents;
  REF_INT neighbor, nneighbor, neighbors[MAX_NODE_LIST];
  REF_INT id, other;

  RAS(ref_node_valid(ref_node, node), "invalid node");
  RAS(ref_node_owned(ref_node, node), "ghost node");

  RUS(REF_EMPTY, ref_interp->cell[node], "no cell for guess");

  RXS(ref_cell_node_list_around(ref_cell, node, MAX_NODE_LIST, &nneighbor,
                                neighbors),
      REF_INCREASE_LIMIT, "neighbors");
  for (neighbor = 0; neighbor < nneighbor; neighbor++) {
    other = neighbors[neighbor];
    if (ref_node_owned(ref_node, other)) {
      if (ref_interp->cell[other] == REF_EMPTY &&
          !(ref_interp->agent_hired[other])) {
        ref_interp->agent_hired[other] = REF_TRUE;
        RSS(ref_agents_push(ref_agents, other, ref_interp->part[node],
                            ref_interp->cell[node],
                            ref_node_xyz_ptr(ref_node, other), &id),
            "enque");
      }
    } else { /* add ghost seeding via REF_AGENT_SUGGESTION mode */
      RSS(ref_agents_push(ref_agents, other, ref_interp->part[node],
                          ref_interp->cell[node],
                          ref_node_xyz_ptr(ref_node, other), &id),
          "enque");
      ref_agent_mode(ref_agents, id) = REF_AGENT_SUGGESTION;
      ref_agent_node(ref_agents, id) = REF_EMPTY;
      ref_agent_home(ref_agents, id) = ref_node_part(ref_node, other);
      ref_agent_global(ref_agents, id) = ref_node_global(ref_node, other);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_process_agents(REF_INTERP ref_interp) {
  REF_NODE from_node = ref_grid_node(ref_interp_from_grid(ref_interp));
  REF_NODE to_node = ref_grid_node(ref_interp_to_grid(ref_interp));
  REF_CELL from_cell = ref_grid_tet(ref_interp_from_grid(ref_interp));
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_AGENTS ref_agents = ref_interp->ref_agents;
  REF_INT i, id, node;
  REF_INT n_agents;
  REF_INT sweep = 0;

  n_agents = ref_agents_n(ref_agents);
  RSS(ref_mpi_allsum(ref_mpi, &n_agents, 1, REF_INT_TYPE), "sum");

  while (n_agents > 0) {
    if (ref_interp->instrument && ref_mpi_once(ref_interp->ref_mpi))
      printf(" %2d sweep", sweep);
    if (ref_interp->instrument) ref_agents_population(ref_agents, "agent pop");
    sweep++;

    each_active_ref_agent(ref_agents, id) {
      if (REF_AGENT_WALKING == ref_agent_mode(ref_agents, id) &&
          ref_agent_part(ref_agents, id) == ref_mpi_rank(ref_mpi)) {
        RSS(ref_interp_walk_agent(ref_interp, id), "walking");
      }
    }

    RSS(ref_agents_migrate(ref_agents), "send it");

    each_active_ref_agent(ref_agents, id) {
      if (REF_AGENT_HOP_PART == ref_agent_mode(ref_agents, id) &&
          ref_agent_part(ref_agents, id) == ref_mpi_rank(ref_mpi)) {
        RSS(ref_node_local(from_node, ref_agent_global(ref_agents, id), &node),
            "localize");
        ref_agent_mode(ref_agents, id) = REF_AGENT_WALKING;
        /* pick best from orbit? */
        ref_agent_seed(ref_agents, id) = ref_cell_first_with(from_cell, node);
      }
    }

    each_active_ref_agent(ref_agents, id) {
      if (REF_AGENT_SUGGESTION == ref_agent_mode(ref_agents, id) &&
          ref_agent_home(ref_agents, id) == ref_mpi_rank(ref_mpi)) {
        RSS(ref_node_local(to_node, ref_agent_global(ref_agents, id), &node),
            "localize");
        if (REF_EMPTY != ref_interp->cell[node] ||
            ref_interp->agent_hired[node]) {
          RSS(ref_agents_remove(ref_interp->ref_agents, id), "already got one");
        } else {
          ref_agent_mode(ref_interp->ref_agents, id) = REF_AGENT_WALKING;
          ref_agent_node(ref_interp->ref_agents, id) = node;
          ref_agent_global(ref_interp->ref_agents, id) = REF_EMPTY;
          ref_interp->agent_hired[node] = REF_TRUE;
        }
      }
    }

    each_active_ref_agent(ref_agents, id) {
      if ((REF_AGENT_AT_BOUNDARY == ref_agent_mode(ref_agents, id) ||
           REF_AGENT_TERMINATED == ref_agent_mode(ref_agents, id)) &&
          ref_agent_home(ref_agents, id) == ref_mpi_rank(ref_mpi)) {
        node = ref_agent_node(ref_agents, id);
        RAS(ref_node_valid(to_node, node), "not vaild");
        RAS(ref_node_owned(to_node, node), "ghost, not owned");
        REIS(REF_EMPTY, ref_interp->cell[node], "already found?");
        RAS(ref_interp->agent_hired[node], "should have an agent");
        if (REF_AGENT_TERMINATED == ref_agent_mode(ref_agents, id)) {
          (ref_interp->walk_steps) += (ref_agent_step(ref_agents, id) + 1);
          (ref_interp->n_terminated)++;
        }
        ref_interp->agent_hired[node] = REF_FALSE; /* but nore more */
        RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");
      }
    }

    each_active_ref_agent(ref_agents, id) {
      if (REF_AGENT_ENCLOSING == ref_agent_mode(ref_agents, id) &&
          ref_agent_home(ref_agents, id) == ref_mpi_rank(ref_mpi)) {
        node = ref_agent_node(ref_agents, id);
        RAS(ref_node_valid(to_node, node), "not vaild");
        RAS(ref_node_owned(to_node, node), "ghost, not owned");
        REIS(REF_EMPTY, ref_interp->cell[node], "already found?");
        RAS(ref_interp->agent_hired[node], "should have an agent");

        ref_interp->cell[node] = ref_agent_seed(ref_agents, id);
        ref_interp->part[node] = ref_agent_part(ref_agents, id);
        for (i = 0; i < 4; i++)
          ref_interp->bary[i + 4 * node] = ref_agent_bary(ref_agents, i, id);
        (ref_interp->walk_steps) += (ref_agent_step(ref_agents, id) + 1);
        (ref_interp->n_walk)++;

        ref_interp->agent_hired[node] = REF_FALSE; /* but nore more */
        RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");
        RSS(ref_interp_push_onto_queue(ref_interp, node), "push");
      }
    }

    n_agents = ref_agents_n(ref_agents);
    RSS(ref_mpi_allsum(ref_mpi, &n_agents, 1, REF_INT_TYPE), "sum");
  }

  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->walk_steps), 1, REF_INT_TYPE),
      "sum");
  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_walk), 1, REF_INT_TYPE), "sum");
  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_terminated), 1, REF_INT_TYPE),
      "sum");

  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node))
      REIS(REF_FALSE, ref_interp->agent_hired[node], "should be done");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_geom_node_list(REF_GRID ref_grid, REF_LIST ref_list) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL edg = ref_grid_edg(ref_grid);
  REF_INT nfaceid, faceids[3];
  REF_INT nedgeid, edgeids[2];
  REF_INT node;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      RXS(ref_cell_id_list_around(edg, node, 2, &nedgeid, edgeids),
          REF_INCREASE_LIMIT, "count faceids");
      RXS(ref_cell_id_list_around(tri, node, 3, &nfaceid, faceids),
          REF_INCREASE_LIMIT, "count faceids");
      if (nfaceid >= 3 || nedgeid >= 2)
        RSS(ref_list_push(ref_list, node), "add geom node");
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_interp_geom_nodes(REF_INTERP ref_interp) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_LIST to_geom_list, from_geom_list;
  REF_INT to_geom_node, from_geom_node;
  REF_INT to_item, from_item;
  REF_DBL *xyz;
  REF_DBL *local_xyz, *global_xyz;
  REF_INT *local_node, *global_node;
  REF_INT total_node, *source, i, *best_node, *from_proc;
  REF_DBL dist, *best_dist;
  REF_INT nsend, nrecv;
  REF_INT *send_proc, *my_proc, *recv_proc;
  REF_INT *send_cell, *recv_cell;
  REF_INT *send_node, *recv_node;
  REF_DBL *send_bary, *recv_bary;

  RSS(ref_list_create(&to_geom_list), "create list");
  RSS(ref_list_create(&from_geom_list), "create list");
  RSS(ref_interp_geom_node_list(to_grid, to_geom_list), "to list");
  RSS(ref_interp_geom_node_list(from_grid, from_geom_list), "from list");

  ref_malloc(local_node, ref_list_n(to_geom_list), REF_INT);
  ref_malloc(local_xyz, 3 * ref_list_n(to_geom_list), REF_DBL);
  each_ref_list_item(to_geom_list, to_item) {
    to_geom_node = ref_list_value(to_geom_list, to_item);
    local_node[to_item] = to_geom_node;
    local_xyz[0 + 3 * to_item] = ref_node_xyz(to_node, 0, to_geom_node);
    local_xyz[1 + 3 * to_item] = ref_node_xyz(to_node, 1, to_geom_node);
    local_xyz[2 + 3 * to_item] = ref_node_xyz(to_node, 2, to_geom_node);
  }
  RSS(ref_mpi_allconcat(ref_mpi, 3, ref_list_n(to_geom_list), (void *)local_xyz,
                        &total_node, &source, (void **)&global_xyz,
                        REF_DBL_TYPE),
      "cat");
  ref_free(source);
  RSS(ref_mpi_allconcat(ref_mpi, 1, ref_list_n(to_geom_list),
                        (void *)local_node, &total_node, &source,
                        (void **)&global_node, REF_INT_TYPE),
      "cat");

  ref_malloc(best_dist, total_node, REF_DBL);
  ref_malloc(best_node, total_node, REF_INT);
  ref_malloc(from_proc, total_node, REF_INT);
  for (to_item = 0; to_item < total_node; to_item++) {
    xyz = &(global_xyz[3 * to_item]);
    best_dist[to_item] = 1.0e20;
    best_node[to_item] = REF_EMPTY;
    each_ref_list_item(from_geom_list, from_item) {
      from_geom_node = ref_list_value(from_geom_list, from_item);
      dist = pow(xyz[0] - ref_node_xyz(from_node, 0, from_geom_node), 2) +
             pow(xyz[1] - ref_node_xyz(from_node, 1, from_geom_node), 2) +
             pow(xyz[2] - ref_node_xyz(from_node, 2, from_geom_node), 2);
      dist = sqrt(dist);
      if (dist < best_dist[to_item] || 0 == from_item) {
        best_dist[to_item] = dist;
        best_node[to_item] = from_geom_node;
      }
    }
  }

  RSS(ref_mpi_allminwho(ref_mpi, best_dist, from_proc, total_node), "who");

  nsend = 0;
  for (to_item = 0; to_item < total_node; to_item++)
    if (ref_mpi_rank(ref_mpi) == from_proc[to_item]) nsend++;

  ref_malloc(send_bary, 4 * nsend, REF_DBL);
  ref_malloc(send_cell, nsend, REF_INT);
  ref_malloc(send_node, nsend, REF_INT);
  ref_malloc(send_proc, nsend, REF_INT);
  ref_malloc_init(my_proc, nsend, REF_INT, ref_mpi_rank(ref_mpi));

  nsend = 0;
  for (to_item = 0; to_item < total_node; to_item++)
    if (ref_mpi_rank(ref_mpi) == from_proc[to_item]) {
      RUS(REF_EMPTY, best_node[to_item], "no geom node");
      xyz = &(global_xyz[3 * to_item]);
      send_node[nsend] = global_node[to_item];
      send_proc[nsend] = source[to_item];
      RSS(ref_interp_exhaustive_tet_around_node(from_grid, best_node[to_item],
                                                xyz, &(send_cell[nsend]),
                                                &(send_bary[4 * nsend])),
          "tet around node");
      nsend++;
    }

  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_node, 1, nsend,
                        (void **)(&recv_node), &nrecv, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_cell, 1, nsend,
                        (void **)(&recv_cell), &nrecv, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)my_proc, 1, nsend,
                        (void **)(&recv_proc), &nrecv, REF_INT_TYPE),
      "blind send proc");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_bary, 4, nsend,
                        (void **)(&recv_bary), &nrecv, REF_DBL_TYPE),
      "blind send bary");

  for (from_item = 0; from_item < nrecv; from_item++) {
    if (recv_bary[0 + 4 * from_item] > ref_interp->inside &&
        recv_bary[1 + 4 * from_item] > ref_interp->inside &&
        recv_bary[2 + 4 * from_item] > ref_interp->inside &&
        recv_bary[3 + 4 * from_item] > ref_interp->inside) {
      ref_interp->n_geom++;
      to_geom_node = recv_node[from_item];
      REIS(REF_EMPTY, ref_interp->cell[to_geom_node], "geom already found?");
      if (ref_interp->agent_hired[to_geom_node]) { /* need to dequeue */
        RSS(ref_agents_delete(ref_interp->ref_agents, to_geom_node), "deq");
        ref_interp->agent_hired[to_geom_node] = REF_FALSE;
      }
      ref_interp->cell[to_geom_node] = recv_cell[from_item];
      ref_interp->part[to_geom_node] = recv_proc[from_item];
      for (i = 0; i < 4; i++)
        ref_interp->bary[i + 4 * to_geom_node] = recv_bary[i + 4 * from_item];
      RSS(ref_interp_push_onto_queue(ref_interp, to_geom_node), "push");
    } else {
      ref_interp->n_geom_fail++;
    }
  }

  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_geom), 1, REF_INT_TYPE), "as");
  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_geom_fail), 1, REF_INT_TYPE),
      "as");

  ref_free(recv_node);
  ref_free(recv_cell);
  ref_free(recv_proc);
  ref_free(recv_bary);

  ref_free(send_bary);
  ref_free(my_proc);
  ref_free(send_cell);
  ref_free(send_node);
  ref_free(send_proc);

  ref_free(from_proc);
  ref_free(best_node);
  ref_free(best_dist);

  ref_free(source);
  ref_free(global_node);
  ref_free(global_xyz);
  ref_free(local_xyz);
  ref_free(local_node);
  ref_list_free(from_geom_list);
  ref_list_free(to_geom_list);
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_tree(REF_INTERP ref_interp,
                                  REF_BOOL *increase_fuzz) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tet = ref_grid_tet(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_SEARCH ref_search = ref_interp_search(ref_interp);
  REF_DBL bary[4];
  REF_LIST ref_list;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, *best_node, *best_cell, *from_proc;
  REF_DBL *best_bary;
  REF_INT ntarget;
  REF_DBL *local_xyz, *global_xyz;
  REF_INT *local_node, *global_node;
  REF_INT *source, total_node;
  REF_INT nsend, nrecv;
  REF_INT *send_proc, *my_proc, *recv_proc;
  REF_INT *send_cell, *recv_cell;
  REF_INT *send_node, *recv_node;
  REF_DBL *send_bary, *recv_bary;
  REF_INT i, item;

  *increase_fuzz = REF_FALSE;

  RSS(ref_list_create(&ref_list), "create list");

  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    ntarget++;
  }
  ref_malloc(local_node, ntarget, REF_INT);
  ref_malloc(local_xyz, 3 * ntarget, REF_DBL);
  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    local_node[ntarget] = node;
    local_xyz[0 + 3 * ntarget] = ref_node_xyz(to_node, 0, node);
    local_xyz[1 + 3 * ntarget] = ref_node_xyz(to_node, 1, node);
    local_xyz[2 + 3 * ntarget] = ref_node_xyz(to_node, 2, node);
    ntarget++;
  }
  RSS(ref_mpi_allconcat(ref_mpi, 3, ntarget, (void *)local_xyz, &total_node,
                        &source, (void **)&global_xyz, REF_DBL_TYPE),
      "cat");
  ref_free(source);
  RSS(ref_mpi_allconcat(ref_mpi, 1, ntarget, (void *)local_node, &total_node,
                        &source, (void **)&global_node, REF_INT_TYPE),
      "cat");

  ref_malloc(best_bary, total_node, REF_DBL);
  ref_malloc(best_node, total_node, REF_INT);
  ref_malloc(best_cell, total_node, REF_INT);
  ref_malloc(from_proc, total_node, REF_INT);
  for (node = 0; node < total_node; node++) {
    best_node[node] = global_node[node];
    best_cell[node] = REF_EMPTY;
    best_bary[node] = 1.0e20; /* negative for min, until use max*/
    RSS(ref_search_touching(ref_search, ref_list, &(global_xyz[3 * node]),
                            ref_interp_search_fuzz(ref_interp)),
        "tch");
    if (ref_list_n(ref_list) > 0) {
      RSS(ref_interp_enclosing_tet_in_list(from_grid, ref_list,
                                           &(global_xyz[3 * node]),
                                           &(best_cell[node]), bary),
          "best in list");
      if (REF_EMPTY != best_cell[node]) {
        /* negative for min, until use max*/
        best_bary[node] = -MIN(MIN(bary[0], bary[1]), MIN(bary[2], bary[3]));
      }
    } else {
      best_cell[node] = REF_EMPTY;
    }
    (ref_interp->tree_cells) += ref_list_n(ref_list);
    RSS(ref_list_erase(ref_list), "reset list");
  }

  /* negative for min, until use max*/
  RSS(ref_mpi_allminwho(ref_mpi, best_bary, from_proc, total_node), "who");

  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) nsend++;

  ref_malloc(send_bary, 4 * nsend, REF_DBL);
  ref_malloc(send_cell, nsend, REF_INT);
  ref_malloc(send_node, nsend, REF_INT);
  ref_malloc(send_proc, nsend, REF_INT);
  ref_malloc_init(my_proc, nsend, REF_INT, ref_mpi_rank(ref_mpi));
  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) {
      send_proc[nsend] = source[node];
      send_node[nsend] = best_node[node];
      send_cell[nsend] = best_cell[node];
      if (REF_EMPTY != send_cell[nsend]) {
        RSB(ref_cell_nodes(from_tet, best_cell[node], nodes),
            "cell should be set and valid", {
              printf("global %d best cell %d best bary %e\n", best_node[node],
                     best_cell[node], best_bary[node]);
            });
        RSS(ref_node_bary4(from_node, nodes, &(global_xyz[3 * node]),
                           &(send_bary[4 * nsend])),
            "bary");
      } else {
        *increase_fuzz =
            REF_TRUE; /* candate not found, try again larger fuzz */
      }
      nsend++;
    }

  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_node, 1, nsend,
                        (void **)(&recv_node), &nrecv, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_cell, 1, nsend,
                        (void **)(&recv_cell), &nrecv, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)my_proc, 1, nsend,
                        (void **)(&recv_proc), &nrecv, REF_INT_TYPE),
      "blind send proc");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_bary, 4, nsend,
                        (void **)(&recv_bary), &nrecv, REF_DBL_TYPE),
      "blind send bary");

  for (item = 0; item < nrecv; item++) {
    (ref_interp->n_tree)++;
    node = recv_node[item];
    REIS(REF_EMPTY, ref_interp->cell[node], "tree already found?");
    if (ref_interp->agent_hired[node]) { /* need to dequeue */
      RSS(ref_agents_delete(ref_interp->ref_agents, node), "deq");
      ref_interp->agent_hired[node] = REF_FALSE;
    }
    ref_interp->cell[node] = recv_cell[item];
    ref_interp->part[node] = recv_proc[item];
    if (REF_EMPTY != recv_cell[item]) {
      for (i = 0; i < 4; i++)
        ref_interp->bary[i + 4 * node] = recv_bary[i + 4 * item];
    } else {
      (ref_interp->n_tree)--;
    }
  }

  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_tree), 1, REF_INT_TYPE), "as");
  RSS(ref_mpi_all_or(ref_mpi, increase_fuzz), "sync status");

  ref_free(recv_node);
  ref_free(recv_cell);
  ref_free(recv_proc);
  ref_free(recv_bary);

  ref_free(send_bary);
  ref_free(my_proc);
  ref_free(send_cell);
  ref_free(send_node);
  ref_free(send_proc);

  ref_free(from_proc);
  ref_free(best_cell);
  ref_free(best_node);
  ref_free(best_bary);

  ref_free(source);
  ref_free(global_node);
  ref_free(global_xyz);
  ref_free(local_xyz);
  ref_free(local_node);

  RSS(ref_list_free(ref_list), "free list");

  if (!(*increase_fuzz)) {
    each_ref_node_valid_node(to_node, node) {
      if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
        continue;
      RUS(REF_EMPTY, ref_interp->cell[node], "node missed by tree");
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_nearest_tri_in_tree(REF_INTERP ref_interp,
                                                 REF_SEARCH ref_search) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tet = ref_grid_tet(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_DBL bary[4];
  REF_LIST ref_list;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, *best_node, *best_cell, *from_proc;
  REF_DBL *best_bary;
  REF_INT ntarget;
  REF_DBL *local_xyz, *global_xyz;
  REF_INT *local_node, *global_node;
  REF_INT *source, total_node;
  REF_INT nsend, nrecv;
  REF_INT *send_proc, *my_proc, *recv_proc;
  REF_INT *send_cell, *recv_cell;
  REF_INT *send_node, *recv_node;
  REF_DBL *send_bary, *recv_bary;
  REF_INT i, item;

  RSS(ref_list_create(&ref_list), "create list");

  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    ntarget++;
  }
  ref_malloc(local_node, ntarget, REF_INT);
  ref_malloc(local_xyz, 3 * ntarget, REF_DBL);
  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    local_node[ntarget] = node;
    local_xyz[0 + 3 * ntarget] = ref_node_xyz(to_node, 0, node);
    local_xyz[1 + 3 * ntarget] = ref_node_xyz(to_node, 1, node);
    local_xyz[2 + 3 * ntarget] = ref_node_xyz(to_node, 2, node);
    ntarget++;
  }
  RSS(ref_mpi_allconcat(ref_mpi, 3, ntarget, (void *)local_xyz, &total_node,
                        &source, (void **)&global_xyz, REF_DBL_TYPE),
      "cat");
  ref_free(source);
  RSS(ref_mpi_allconcat(ref_mpi, 1, ntarget, (void *)local_node, &total_node,
                        &source, (void **)&global_node, REF_INT_TYPE),
      "cat");

  ref_malloc(best_bary, total_node, REF_DBL);
  ref_malloc(best_node, total_node, REF_INT);
  ref_malloc(best_cell, total_node, REF_INT);
  ref_malloc(from_proc, total_node, REF_INT);
  for (node = 0; node < total_node; node++) {
    best_node[node] = global_node[node];
    best_cell[node] = REF_EMPTY;
    best_bary[node] = 1.0e20; /* negative for min, until use max*/
    RSS(ref_search_nearest_candidates(ref_search, ref_list,
                                      &(global_xyz[3 * node])),
        "near canidates");
    if (ref_list_n(ref_list) > 0) {
      RSS(ref_interp_best_tri_in_list(from_grid, ref_list,
                                      &(global_xyz[3 * node]),
                                      &(best_cell[node]), bary),
          "best in list");
      if (REF_EMPTY != best_cell[node]) {
        /* negative for min, until use max*/
        best_bary[node] = -MIN(MIN(bary[0], bary[1]), MIN(bary[2], bary[3]));
      }
    } else {
      best_cell[node] = REF_EMPTY;
    }
    (ref_interp->tree_cells) += ref_list_n(ref_list);
    RSS(ref_list_erase(ref_list), "reset list");
  }

  /* negative for min, until use max*/
  RSS(ref_mpi_allminwho(ref_mpi, best_bary, from_proc, total_node), "who");

  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) nsend++;

  ref_malloc(send_bary, 4 * nsend, REF_DBL);
  ref_malloc(send_cell, nsend, REF_INT);
  ref_malloc(send_node, nsend, REF_INT);
  ref_malloc(send_proc, nsend, REF_INT);
  ref_malloc_init(my_proc, nsend, REF_INT, ref_mpi_rank(ref_mpi));
  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) {
      send_proc[nsend] = source[node];
      send_node[nsend] = best_node[node];
      send_cell[nsend] = best_cell[node];
      if (REF_EMPTY != send_cell[nsend]) {
        RSB(ref_cell_nodes(from_tet, best_cell[node], nodes),
            "cell should be set and valid", {
              printf("global %d best cell %d best bary %e\n", best_node[node],
                     best_cell[node], best_bary[node]);
            });
        RSS(ref_node_bary4(from_node, nodes, &(global_xyz[3 * node]),
                           &(send_bary[4 * nsend])),
            "bary");
      }
      nsend++;
    }

  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_node, 1, nsend,
                        (void **)(&recv_node), &nrecv, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_cell, 1, nsend,
                        (void **)(&recv_cell), &nrecv, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)my_proc, 1, nsend,
                        (void **)(&recv_proc), &nrecv, REF_INT_TYPE),
      "blind send proc");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_bary, 4, nsend,
                        (void **)(&recv_bary), &nrecv, REF_DBL_TYPE),
      "blind send bary");

  for (item = 0; item < nrecv; item++) {
    (ref_interp->n_tree)++;
    node = recv_node[item];
    REIS(REF_EMPTY, ref_interp->cell[node], "tree already found?");
    if (ref_interp->agent_hired[node]) { /* need to dequeue */
      RSS(ref_agents_delete(ref_interp->ref_agents, node), "deq");
      ref_interp->agent_hired[node] = REF_FALSE;
    }
    ref_interp->cell[node] = recv_cell[item];
    ref_interp->part[node] = recv_proc[item];
    if (REF_EMPTY != recv_cell[item]) {
      for (i = 0; i < 4; i++)
        ref_interp->bary[i + 4 * node] = recv_bary[i + 4 * item];
    } else {
      (ref_interp->n_tree)--;
    }
  }

  RSS(ref_mpi_allsum(ref_mpi, &(ref_interp->n_tree), 1, REF_INT_TYPE), "as");

  ref_free(recv_node);
  ref_free(recv_cell);
  ref_free(recv_proc);
  ref_free(recv_bary);

  ref_free(send_bary);
  ref_free(my_proc);
  ref_free(send_cell);
  ref_free(send_node);
  ref_free(send_proc);

  ref_free(from_proc);
  ref_free(best_cell);
  ref_free(best_node);
  ref_free(best_bary);

  ref_free(source);
  ref_free(global_node);
  ref_free(global_xyz);
  ref_free(local_xyz);
  ref_free(local_node);

  RSS(ref_list_free(ref_list), "free list");

  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      RUS(REF_EMPTY, ref_interp->cell[node], "node missed by tree");
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_seed_tree(REF_INTERP ref_interp) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tet = ref_grid_tet(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_SEARCH ref_search = ref_interp_search(ref_interp);
  REF_DBL bary[4];
  REF_LIST ref_list;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, *best_node, *best_cell, *from_proc;
  REF_DBL *best_bary;
  REF_INT ntarget, seed_target, stride;
  REF_DBL *local_xyz, *global_xyz;
  REF_INT *local_node, *global_node;
  REF_INT *source, total_node;
  REF_INT nsend, nrecv;
  REF_INT *send_proc, *my_proc, *recv_proc;
  REF_INT *send_cell, *recv_cell;
  REF_INT *send_node, *recv_node;
  REF_DBL *send_bary, *recv_bary;
  REF_INT i, item;

  RSS(ref_list_create(&ref_list), "create list");

  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    ntarget++;
  }
  ref_malloc(local_node, ntarget, REF_INT);
  ref_malloc(local_xyz, 3 * ntarget, REF_DBL);
  ntarget = 0;
  each_ref_node_valid_node(to_node, node) {
    if (!ref_node_owned(to_node, node) || REF_EMPTY != ref_interp->cell[node])
      continue;
    local_node[ntarget] = node;
    local_xyz[0 + 3 * ntarget] = ref_node_xyz(to_node, 0, node);
    local_xyz[1 + 3 * ntarget] = ref_node_xyz(to_node, 1, node);
    local_xyz[2 + 3 * ntarget] = ref_node_xyz(to_node, 2, node);
    ntarget++;
  }

  /* use approximately 1 percent of existing targets sampled uniformly,
     betwwen 10 or 100, less than or equal ntarget */
  seed_target = ntarget / 100;
  seed_target = MAX(seed_target, 10);
  seed_target = MIN(seed_target, 100);
  seed_target = MIN(seed_target, ntarget);
  if (seed_target > 0) {
    stride = 1 + ntarget / seed_target;
    seed_target = 0;
    for (node = 0; node < ntarget; node += stride) {
      local_node[seed_target] = local_node[node];
      local_xyz[0 + 3 * seed_target] = local_xyz[0 + 3 * node];
      local_xyz[1 + 3 * seed_target] = local_xyz[1 + 3 * node];
      local_xyz[2 + 3 * seed_target] = local_xyz[2 + 3 * node];
      seed_target++;
    }
    ntarget = seed_target;
  }

  RSS(ref_mpi_allconcat(ref_mpi, 3, ntarget, (void *)local_xyz, &total_node,
                        &source, (void **)&global_xyz, REF_DBL_TYPE),
      "cat");
  ref_free(source);
  RSS(ref_mpi_allconcat(ref_mpi, 1, ntarget, (void *)local_node, &total_node,
                        &source, (void **)&global_node, REF_INT_TYPE),
      "cat");

  ref_malloc(best_bary, total_node, REF_DBL);
  ref_malloc(best_node, total_node, REF_INT);
  ref_malloc(best_cell, total_node, REF_INT);
  ref_malloc(from_proc, total_node, REF_INT);
  for (node = 0; node < total_node; node++) {
    best_node[node] = global_node[node];
    best_cell[node] = REF_EMPTY;
    best_bary[node] = 1.0e20; /* negative for min, until use max*/
    RSS(ref_search_touching(ref_search, ref_list, &(global_xyz[3 * node]),
                            ref_interp_search_fuzz(ref_interp)),
        "tch");
    if (ref_list_n(ref_list) > 0) {
      RSS(ref_interp_enclosing_tet_in_list(from_grid, ref_list,
                                           &(global_xyz[3 * node]),
                                           &(best_cell[node]), bary),
          "best in list");
      if (REF_EMPTY != best_cell[node]) {
        /* negative for min, until use max*/
        best_bary[node] = -MIN(MIN(bary[0], bary[1]), MIN(bary[2], bary[3]));
      }
    } else {
      best_cell[node] = REF_EMPTY;
    }
    (ref_interp->tree_cells) += ref_list_n(ref_list);
    RSS(ref_list_erase(ref_list), "reset list");
  }

  /* negative for min, until use max*/
  RSS(ref_mpi_allminwho(ref_mpi, best_bary, from_proc, total_node), "who");

  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) nsend++;

  ref_malloc(send_bary, 4 * nsend, REF_DBL);
  ref_malloc(send_cell, nsend, REF_INT);
  ref_malloc(send_node, nsend, REF_INT);
  ref_malloc(send_proc, nsend, REF_INT);
  ref_malloc_init(my_proc, nsend, REF_INT, ref_mpi_rank(ref_mpi));
  nsend = 0;
  for (node = 0; node < total_node; node++)
    if (ref_mpi_rank(ref_mpi) == from_proc[node]) {
      send_proc[nsend] = source[node];
      send_node[nsend] = best_node[node];
      send_cell[nsend] = best_cell[node];
      if (REF_EMPTY != send_cell[nsend]) {
        RSB(ref_cell_nodes(from_tet, best_cell[node], nodes),
            "cell should be set and valid", {
              printf("global %d best cell %d best bary %e\n", best_node[node],
                     best_cell[node], best_bary[node]);
            });
        RSS(ref_node_bary4(from_node, nodes, &(global_xyz[3 * node]),
                           &(send_bary[4 * nsend])),
            "bary");
      }
      nsend++;
    }

  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_node, 1, nsend,
                        (void **)(&recv_node), &nrecv, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_cell, 1, nsend,
                        (void **)(&recv_cell), &nrecv, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)my_proc, 1, nsend,
                        (void **)(&recv_proc), &nrecv, REF_INT_TYPE),
      "blind send proc");
  RSS(ref_mpi_blindsend(ref_mpi, send_proc, (void *)send_bary, 4, nsend,
                        (void **)(&recv_bary), &nrecv, REF_DBL_TYPE),
      "blind send bary");

  for (item = 0; item < nrecv; item++) {
    (ref_interp->n_tree)++;
    node = recv_node[item];
    REIS(REF_EMPTY, ref_interp->cell[node], "tree already found?");
    if (ref_interp->agent_hired[node]) { /* need to dequeue */
      RSS(ref_agents_delete(ref_interp->ref_agents, node), "deq");
      ref_interp->agent_hired[node] = REF_FALSE;
    }
    ref_interp->cell[node] = recv_cell[item];
    ref_interp->part[node] = recv_proc[item];
    if (REF_EMPTY != recv_cell[item]) {
      for (i = 0; i < 4; i++)
        ref_interp->bary[i + 4 * node] = recv_bary[i + 4 * item];
      RSS(ref_interp_push_onto_queue(ref_interp, node), "queue neighbors");
    } else {
      (ref_interp->n_tree)--;
    }
  }

  ref_free(recv_node);
  ref_free(recv_cell);
  ref_free(recv_proc);
  ref_free(recv_bary);

  ref_free(send_bary);
  ref_free(my_proc);
  ref_free(send_cell);
  ref_free(send_node);
  ref_free(send_proc);

  ref_free(from_proc);
  ref_free(best_cell);
  ref_free(best_node);
  ref_free(best_bary);

  ref_free(source);
  ref_free(global_node);
  ref_free(global_xyz);
  ref_free(local_xyz);
  ref_free(local_node);

  RSS(ref_list_free(ref_list), "free list");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate(REF_INTERP ref_interp) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_BOOL increase_fuzz;
  REF_INT tries;

  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_start(ref_mpi), "locate clock");

  RSS(ref_interp_geom_nodes(ref_interp), "geom nodes");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "geom"), "locate clock");

  RSS(ref_interp_process_agents(ref_interp), "drain");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "drain"), "locate clock");

  increase_fuzz = REF_FALSE;
  for (tries = 0; tries < 12; tries++) {
    if (increase_fuzz) {
      ref_interp_search_fuzz(ref_interp) *= 10.0;
      if (ref_mpi_once(ref_mpi))
        printf("retry tree search with %e fuzz\n",
               ref_interp_search_fuzz(ref_interp));
    }
    RSS(ref_interp_tree(ref_interp, &increase_fuzz), "tree");
    if (ref_interp->instrument)
      RSS(ref_mpi_stopwatch_stop(ref_mpi, "tree"), "locate clock");
    if (!increase_fuzz) break;
  }
  REIS(REF_FALSE, increase_fuzz, "unable to grow fuzz to find tree candidate");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate_warm(REF_INTERP ref_interp) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE to_node = ref_grid_node(ref_interp_to_grid(ref_interp));
  REF_BOOL increase_fuzz;
  REF_INT tries;
  REF_INT node;

  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_start(ref_mpi), "locate clock");

  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node) && REF_EMPTY != ref_interp->cell[node]) {
      RSS(ref_interp_push_onto_queue(ref_interp, node), "queue neighbors");
    }
  }

  RSS(ref_interp_process_agents(ref_interp), "drain");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "drain"), "locate clock");

  increase_fuzz = REF_FALSE;
  for (tries = 0; tries < 12; tries++) {
    if (increase_fuzz) {
      ref_interp_search_fuzz(ref_interp) *= 10.0;
      if (ref_mpi_once(ref_mpi))
        printf("retry tree search with %e fuzz\n",
               ref_interp_search_fuzz(ref_interp));
    }
    RSS(ref_interp_tree(ref_interp, &increase_fuzz), "tree");
    if (ref_interp->instrument)
      RSS(ref_mpi_stopwatch_stop(ref_mpi, "tree"), "locate clock");
    if (!increase_fuzz) break;
  }
  REIS(REF_FALSE, increase_fuzz, "unable to grow fuzz to find tree candidate");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate_subset(REF_INTERP ref_interp) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_BOOL increase_fuzz;
  REF_INT tries;

  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_start(ref_mpi), "locate clock");

  RSS(ref_interp_seed_tree(ref_interp), "seed tree nodes");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "seed tree"), "locate clock");

  RSS(ref_interp_process_agents(ref_interp), "drain");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "drain"), "locate clock");

  increase_fuzz = REF_FALSE;
  for (tries = 0; tries < 12; tries++) {
    if (increase_fuzz) {
      ref_interp_search_fuzz(ref_interp) *= 10.0;
      if (ref_mpi_once(ref_mpi))
        printf("retry tree search with %e fuzz\n",
               ref_interp_search_fuzz(ref_interp));
    }
    RSS(ref_interp_tree(ref_interp, &increase_fuzz), "tree");
    if (ref_interp->instrument)
      RSS(ref_mpi_stopwatch_stop(ref_mpi, "tree"), "locate clock");
    if (!increase_fuzz) break;
  }
  REIS(REF_FALSE, increase_fuzz, "unable to grow fuzz to find tree candidate");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate_nearest(REF_INTERP ref_interp) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);

  REF_CELL from_tri = ref_grid_tri(from_grid);
  REF_NODE from_node = ref_grid_node(from_grid);

  REF_BOOL increase_fuzz;
  REF_SEARCH ref_search;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius;

  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_start(ref_mpi), "locate clock");

  RSS(ref_interp_seed_tree(ref_interp), "seed tree nodes");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "seed tree"), "locate clock");

  RSS(ref_interp_process_agents(ref_interp), "drain");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "drain"), "locate clock");

  increase_fuzz = REF_FALSE;
  RSS(ref_interp_tree(ref_interp, &increase_fuzz), "tree");
  if (ref_interp->instrument)
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "tree"), "locate clock");

  if (increase_fuzz) {
    RSS(ref_search_create(&ref_search, ref_cell_n(from_tri)), "create search");
    each_ref_cell_valid_cell_with_nodes(from_tri, cell, nodes) {
      RSS(ref_interp_bounding_sphere3(from_node, nodes, center, &radius), "b");
      RSS(ref_search_insert(ref_search, cell, center,
                            ref_interp_search_donor_scale(ref_interp) * radius),
          "ins");
    }
    RSS(ref_interp_nearest_tri_in_tree(ref_interp, ref_search), "near tri");
    RSS(ref_search_free(ref_search), "free search");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate_node(REF_INTERP ref_interp, REF_INT node) {
  REF_MPI ref_mpi;
  REF_NODE ref_node;
  REF_AGENTS ref_agents;
  REF_INT i, id;
  RNS(ref_interp, "ref_interp NULL");
  ref_mpi = ref_interp_mpi(ref_interp);

  RAS(node <= ref_interp_max(ref_interp), "more nodes added, should move only");

  /* no starting guess, skip */
  if (REF_EMPTY == ref_interp->cell[node] ||
      ref_mpi_rank(ref_mpi) != ref_interp->part[node])
    return REF_SUCCESS;

  ref_node = ref_grid_node(ref_interp_to_grid(ref_interp));
  ref_agents = ref_interp->ref_agents;
  REIS(0, ref_agents_n(ref_agents), "did not expect active agents");

  ref_interp->agent_hired[node] = REF_TRUE;
  RSS(ref_agents_push(ref_agents, node, ref_interp->part[node],
                      ref_interp->cell[node], ref_node_xyz_ptr(ref_node, node),
                      &id),
      "requeue");
  REIS(REF_AGENT_WALKING, ref_agent_mode(ref_agents, id), "should be walking");
  RSS(ref_interp_walk_agent(ref_interp, id), "walking");
  if (REF_AGENT_ENCLOSING == ref_agent_mode(ref_agents, id)) {
    ref_interp->cell[node] = ref_agent_seed(ref_agents, id);
    ref_interp->part[node] = ref_agent_part(ref_agents, id);
    for (i = 0; i < 4; i++)
      ref_interp->bary[i + 4 * node] = ref_agent_bary(ref_agents, i, id);
    (ref_interp->walk_steps) += (ref_agent_step(ref_agents, id) + 1);
    (ref_interp->n_walk)++;
    REIS(ref_mpi_rank(ref_mpi), ref_interp->part[node], "expected local");
    RAS(ref_cell_valid(ref_grid_tet(ref_interp_from_grid(ref_interp)),
                       ref_interp->cell[node]),
        "expected a valid cell");
  } else {
    /* new seed or go exhaustive for REF_AGENT_AT_BOUNDARY */
    /* what for parallel REF_AGENT_HOP_PART */
    ref_interp->cell[node] = REF_EMPTY;
  }
  ref_interp->agent_hired[node] = REF_FALSE; /* dismissed */
  RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");

  if (!ref_mpi_para(ref_node_mpi(ref_node)) &&
      REF_EMPTY == ref_interp->cell[node]) {
    REF_LIST ref_list;

    RSS(ref_list_create(&ref_list), "create list");
    RSS(ref_search_touching(ref_interp_search(ref_interp), ref_list,
                            ref_node_xyz_ptr(ref_node, node),
                            ref_interp_search_fuzz(ref_interp)),
        "tch");
    if (ref_list_n(ref_list) > 0) {
      RSS(ref_interp_enclosing_tet_in_list(
              ref_interp_from_grid(ref_interp), ref_list,
              ref_node_xyz_ptr(ref_node, node), &(ref_interp->cell[node]),
              &(ref_interp->bary[4 * node])),
          "best in list");
    }
    RSS(ref_list_free(ref_list), "free list");
  }

  if (REF_EMPTY == ref_interp->cell[node]) {
    return REF_NOT_FOUND;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate_between(REF_INTERP ref_interp, REF_INT node0,
                                     REF_INT node1, REF_INT new_node) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_MPI ref_mpi;
  REF_AGENTS ref_agents;
  REF_INT i, id;
  RNS(ref_interp, "ref_interp NULL");

  ref_mpi = ref_interp_mpi(ref_interp);
  ref_grid = ref_interp_to_grid(ref_interp);
  ref_node = ref_grid_node(ref_grid);
  if (new_node >= ref_interp_max(ref_interp)) {
    RSS(ref_interp_resize(ref_interp, ref_node_max(ref_node)), "resize");
  }
  ref_interp->cell[new_node] = REF_EMPTY; /* initialize new_node locate */

  ref_agents = ref_interp->ref_agents;
  REIS(0, ref_agents_n(ref_agents), "did not expect active agents");

  id = REF_EMPTY;

  if (id == REF_EMPTY && REF_EMPTY != ref_interp->cell[node0] &&
      ref_mpi_rank(ref_mpi) == ref_interp->part[node0]) {
    ref_interp->agent_hired[new_node] = REF_TRUE;
    RSS(ref_agents_push(ref_agents, new_node, ref_interp->part[node0],
                        ref_interp->cell[node0],
                        ref_node_xyz_ptr(ref_node, new_node), &id),
        "requeue");
    REIS(REF_AGENT_WALKING, ref_agent_mode(ref_agents, id), "not walking?");
    RSS(ref_interp_walk_agent(ref_interp, id), "walking");
    if (REF_AGENT_ENCLOSING != ref_agent_mode(ref_agents, id)) {
      ref_interp->agent_hired[new_node] = REF_FALSE; /* dismissed */
      RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");
      id = REF_EMPTY;
    }
  }

  if (id == REF_EMPTY && REF_EMPTY != ref_interp->cell[node1] &&
      ref_mpi_rank(ref_mpi) == ref_interp->part[node1]) {
    ref_interp->agent_hired[new_node] = REF_TRUE;
    RSS(ref_agents_push(ref_agents, new_node, ref_interp->part[node1],
                        ref_interp->cell[node1],
                        ref_node_xyz_ptr(ref_node, new_node), &id),
        "requeue");
    REIS(REF_AGENT_WALKING, ref_agent_mode(ref_agents, id), "not walking?");
    RSS(ref_interp_walk_agent(ref_interp, id), "walking");
    if (REF_AGENT_ENCLOSING != ref_agent_mode(ref_agents, id)) {
      ref_interp->agent_hired[new_node] = REF_FALSE; /* dismissed */
      RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");
      id = REF_EMPTY;
    }
  }

  if (id != REF_EMPTY) {
    if (REF_AGENT_ENCLOSING == ref_agent_mode(ref_agents, id)) {
      ref_interp->cell[new_node] = ref_agent_seed(ref_agents, id);
      ref_interp->part[new_node] = ref_agent_part(ref_agents, id);
      for (i = 0; i < 4; i++)
        ref_interp->bary[i + 4 * new_node] = ref_agent_bary(ref_agents, i, id);
      (ref_interp->walk_steps) += (ref_agent_step(ref_agents, id) + 1);
      (ref_interp->n_walk)++;
    } else {
      /* new seed or go exhaustive for REF_AGENT_AT_BOUNDARY */
      /* what for parallel REF_AGENT_HOP_PART */
      ref_interp->cell[new_node] = REF_EMPTY;
    }
    ref_interp->agent_hired[new_node] = REF_FALSE; /* dismissed */
    RSS(ref_agents_remove(ref_agents, id), "no longer neeeded");
  }

  if (!ref_mpi_para(ref_mpi) && REF_EMPTY == ref_interp->cell[new_node]) {
    REF_LIST ref_list;

    RSS(ref_list_create(&ref_list), "create list");
    RSS(ref_search_touching(ref_interp_search(ref_interp), ref_list,
                            ref_node_xyz_ptr(ref_node, new_node),
                            ref_interp_search_fuzz(ref_interp)),
        "tch");
    if (ref_list_n(ref_list) > 0) {
      RSS(ref_interp_enclosing_tet_in_list(
              ref_interp_from_grid(ref_interp), ref_list,
              ref_node_xyz_ptr(ref_node, new_node),
              &(ref_interp->cell[new_node]), &(ref_interp->bary[4 * new_node])),
          "best in list");
    }
    RSS(ref_list_free(ref_list), "free list");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_scalar(REF_INTERP ref_interp, REF_INT leading_dim,
                             REF_DBL *from_scalar, REF_DBL *to_scalar) {
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_MPI ref_mpi = ref_grid_mpi(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_INT node, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT receptor, n_recept, donation, n_donor;
  REF_DBL *recept_scalar, *donor_scalar, *recept_bary, *donor_bary;
  REF_INT *donor_node, *donor_ret, *donor_cell;
  REF_INT *recept_proc, *recept_ret, *recept_node, *recept_cell;

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      n_recept++;
    }
  }

  ref_malloc(recept_bary, 4 * n_recept, REF_DBL);
  ref_malloc(recept_cell, n_recept, REF_INT);
  ref_malloc(recept_node, n_recept, REF_INT);
  ref_malloc(recept_ret, n_recept, REF_INT);
  ref_malloc(recept_proc, n_recept, REF_INT);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      RUS(REF_EMPTY, ref_interp->cell[node], "node needs to be localized");
      RSS(ref_node_clip_bary4(&(ref_interp->bary[4 * node]),
                              &(recept_bary[4 * n_recept])),
          "clip");
      recept_proc[n_recept] = ref_interp->part[node];
      recept_cell[n_recept] = ref_interp->cell[node];
      recept_node[n_recept] = node;
      recept_ret[n_recept] = ref_mpi_rank(ref_mpi);
      n_recept++;
    }
  }

  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_cell, 1, n_recept,
                        (void **)(&donor_cell), &n_donor, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_ret, 1, n_recept,
                        (void **)(&donor_ret), &n_donor, REF_INT_TYPE),
      "blind send ret");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_node, 1, n_recept,
                        (void **)(&donor_node), &n_donor, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_bary, 4, n_recept,
                        (void **)(&donor_bary), &n_donor, REF_DBL_TYPE),
      "blind send bary");

  ref_free(recept_proc);
  ref_free(recept_ret);
  ref_free(recept_node);
  ref_free(recept_cell);
  ref_free(recept_bary);

  ref_malloc(donor_scalar, leading_dim * n_donor, REF_DBL);

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (ibary = 0; ibary < 4; ibary++) {
      for (im = 0; im < leading_dim; im++) {
        donor_scalar[im + leading_dim * donation] = 0.0;
        for (ibary = 0; ibary < 4; ibary++) {
          donor_scalar[im + leading_dim * donation] +=
              donor_bary[ibary + 4 * donation] *
              from_scalar[im + leading_dim * nodes[ibary]];
        }
      }
    }
  }
  ref_free(donor_cell);
  ref_free(donor_bary);

  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_scalar, leading_dim,
                        n_donor, (void **)(&recept_scalar), &n_recept,
                        REF_DBL_TYPE),
      "blind send bary");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_node, 1, n_donor,
                        (void **)(&recept_node), &n_recept, REF_INT_TYPE),
      "blind send node");
  ref_free(donor_scalar);
  ref_free(donor_node);
  ref_free(donor_ret);

  for (receptor = 0; receptor < n_recept; receptor++) {
    node = recept_node[receptor];
    for (im = 0; im < leading_dim; im++) {
      to_scalar[im + leading_dim * node] =
          recept_scalar[im + leading_dim * receptor];
    }
  }

  ref_free(recept_node);
  ref_free(recept_scalar);

  RSS(ref_node_ghost_dbl(to_node, to_scalar, leading_dim), "ghost");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_face_only(REF_INTERP ref_interp, REF_INT faceid,
                                REF_INT leading_dim, REF_DBL *from_scalar,
                                REF_DBL *to_scalar) {
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_MPI ref_mpi = ref_grid_mpi(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_INT node, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT receptor, n_recept, donation, n_donor;
  REF_DBL *recept_scalar, *donor_scalar, *recept_bary, *donor_bary;
  REF_INT *donor_node, *donor_ret, *donor_cell;
  REF_INT *recept_proc, *recept_ret, *recept_node, *recept_cell;
  REF_GLOB nnode, *l2c;
  REF_LONG ncell;
  REF_BOOL increase_fuzz;
  REF_INT tries;

  RSS(ref_grid_compact_cell_id_nodes(to_grid, ref_grid_tri(to_grid), faceid,
                                     &nnode, &ncell, &l2c),
      "l2c");

  each_ref_node_valid_node(to_node, node) {
    if (REF_EMPTY == l2c[node] || !ref_node_owned(to_node, node)) {
      ref_interp->cell[node] = -57; /* magic number to skip */
    }
  }

  increase_fuzz = REF_FALSE;
  for (tries = 0; tries < 12; tries++) {
    if (increase_fuzz) {
      ref_interp_search_fuzz(ref_interp) *= 10.0;
      if (ref_mpi_once(ref_mpi))
        printf("retry tree search with %e fuzz\n",
               ref_interp_search_fuzz(ref_interp));
    }
    RSS(ref_interp_tree(ref_interp, &increase_fuzz), "tree");
    if (ref_interp->instrument)
      RSS(ref_mpi_stopwatch_stop(ref_mpi, "tree"), "locate clock");
    if (!increase_fuzz) break;
  }
  REIS(REF_FALSE, increase_fuzz, "unable to grow fuzz to find tree candidate");

  each_ref_node_valid_node(to_node, node) {
    if (REF_EMPTY == l2c[node] || !ref_node_owned(to_node, node)) {
      ref_interp->cell[node] = REF_EMPTY; /* set back */
    }
  }

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (REF_EMPTY != l2c[node] && ref_node_owned(to_node, node)) {
      n_recept++;
    }
  }

  ref_malloc(recept_bary, 4 * n_recept, REF_DBL);
  ref_malloc(recept_cell, n_recept, REF_INT);
  ref_malloc(recept_node, n_recept, REF_INT);
  ref_malloc(recept_ret, n_recept, REF_INT);
  ref_malloc(recept_proc, n_recept, REF_INT);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (REF_EMPTY != l2c[node] && ref_node_owned(to_node, node)) {
      RUS(REF_EMPTY, ref_interp->cell[node], "node needs to be localized");
      RSS(ref_node_clip_bary4(&(ref_interp->bary[4 * node]),
                              &(recept_bary[4 * n_recept])),
          "clip");
      recept_proc[n_recept] = ref_interp->part[node];
      recept_cell[n_recept] = ref_interp->cell[node];
      recept_node[n_recept] = node;
      recept_ret[n_recept] = ref_mpi_rank(ref_mpi);
      n_recept++;
    }
  }

  ref_free(l2c);

  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_cell, 1, n_recept,
                        (void **)(&donor_cell), &n_donor, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_ret, 1, n_recept,
                        (void **)(&donor_ret), &n_donor, REF_INT_TYPE),
      "blind send ret");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_node, 1, n_recept,
                        (void **)(&donor_node), &n_donor, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_bary, 4, n_recept,
                        (void **)(&donor_bary), &n_donor, REF_DBL_TYPE),
      "blind send bary");

  ref_free(recept_proc);
  ref_free(recept_ret);
  ref_free(recept_node);
  ref_free(recept_cell);
  ref_free(recept_bary);

  ref_malloc(donor_scalar, leading_dim * n_donor, REF_DBL);

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (ibary = 0; ibary < 4; ibary++) {
      for (im = 0; im < leading_dim; im++) {
        donor_scalar[im + leading_dim * donation] = 0.0;
        for (ibary = 0; ibary < 4; ibary++) {
          donor_scalar[im + leading_dim * donation] +=
              donor_bary[ibary + 4 * donation] *
              from_scalar[im + leading_dim * nodes[ibary]];
        }
      }
    }
  }
  ref_free(donor_cell);
  ref_free(donor_bary);

  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_scalar, leading_dim,
                        n_donor, (void **)(&recept_scalar), &n_recept,
                        REF_DBL_TYPE),
      "blind send bary");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_node, 1, n_donor,
                        (void **)(&recept_node), &n_recept, REF_INT_TYPE),
      "blind send node");
  ref_free(donor_scalar);
  ref_free(donor_node);
  ref_free(donor_ret);

  for (receptor = 0; receptor < n_recept; receptor++) {
    node = recept_node[receptor];
    for (im = 0; im < leading_dim; im++) {
      to_scalar[im + leading_dim * node] =
          recept_scalar[im + leading_dim * receptor];
    }
  }

  ref_free(recept_node);
  ref_free(recept_scalar);

  RSS(ref_node_ghost_dbl(to_node, to_scalar, leading_dim), "ghost");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_min_bary(REF_INTERP ref_interp, REF_DBL *min_bary) {
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT node;
  REF_DBL this_bary;

  *min_bary = 1.0;

  each_ref_node_valid_node(to_node, node) if (ref_node_owned(to_node, node)) {
    RUS(REF_EMPTY, ref_interp->cell[node], "node needs to be localized");
    this_bary = MIN(
        MIN(ref_interp->bary[0 + 4 * node], ref_interp->bary[1 + 4 * node]),
        MIN(ref_interp->bary[2 + 4 * node], ref_interp->bary[3 + 4 * node]));
    *min_bary = MIN(*min_bary, this_bary);
  }
  this_bary = *min_bary;
  RSS(ref_mpi_min(ref_mpi, &this_bary, min_bary, REF_DBL_TYPE), "min");
  RSS(ref_mpi_bcast(ref_mpi, min_bary, 1, REF_DBL_TYPE), "bcast");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_max_error(REF_INTERP ref_interp, REF_DBL *max_error) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL error;
  REF_INT i, j;
  REF_INT receptor, n_recept, donation, n_donor;
  REF_DBL *recept_xyz, *donor_xyz, *recept_bary, *donor_bary;
  REF_INT *donor_node, *donor_ret, *donor_cell;
  REF_INT *recept_proc, *recept_ret, *recept_node, *recept_cell;

  if (ref_grid_twod(from_grid)) from_cell = ref_grid_tri(from_grid);

  *max_error = 0.0;

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node) && REF_EMPTY != ref_interp->cell[node]) {
      n_recept++;
    }
  }

  ref_malloc(recept_bary, 4 * n_recept, REF_DBL);
  ref_malloc(recept_cell, n_recept, REF_INT);
  ref_malloc(recept_node, n_recept, REF_INT);
  ref_malloc(recept_ret, n_recept, REF_INT);
  ref_malloc(recept_proc, n_recept, REF_INT);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node) && REF_EMPTY != ref_interp->cell[node]) {
      for (i = 0; i < 4; i++) {
        recept_bary[i + 4 * n_recept] = ref_interp->bary[i + 4 * node];
      }
      recept_proc[n_recept] = ref_interp->part[node];
      recept_cell[n_recept] = ref_interp->cell[node];
      recept_node[n_recept] = node;
      recept_ret[n_recept] = ref_mpi_rank(ref_mpi);
      n_recept++;
    }
  }

  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_cell, 1, n_recept,
                        (void **)(&donor_cell), &n_donor, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_ret, 1, n_recept,
                        (void **)(&donor_ret), &n_donor, REF_INT_TYPE),
      "blind send ret");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_node, 1, n_recept,
                        (void **)(&donor_node), &n_donor, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_bary, 4, n_recept,
                        (void **)(&donor_bary), &n_donor, REF_DBL_TYPE),
      "blind send bary");

  ref_free(recept_proc);
  ref_free(recept_ret);
  ref_free(recept_node);
  ref_free(recept_cell);
  ref_free(recept_bary);

  ref_malloc(donor_xyz, 3 * n_donor, REF_DBL);

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (i = 0; i < 3; i++) {
      donor_xyz[i + 3 * donation] = 0.0;
      for (j = 0; j < ref_cell_node_per(from_cell); j++) {
        donor_xyz[i + 3 * donation] +=
            donor_bary[j + 4 * donation] * ref_node_xyz(from_node, i, nodes[j]);
      }
    }
  }
  ref_free(donor_cell);
  ref_free(donor_bary);

  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_xyz, 3, n_donor,
                        (void **)(&recept_xyz), &n_recept, REF_DBL_TYPE),
      "blind send bary");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_node, 1, n_donor,
                        (void **)(&recept_node), &n_recept, REF_INT_TYPE),
      "blind send node");
  ref_free(donor_xyz);
  ref_free(donor_node);
  ref_free(donor_ret);

  for (receptor = 0; receptor < n_recept; receptor++) {
    node = recept_node[receptor];
    error =
        pow(recept_xyz[0 + 3 * receptor] - ref_node_xyz(to_node, 0, node), 2) +
        pow(recept_xyz[1 + 3 * receptor] - ref_node_xyz(to_node, 1, node), 2) +
        pow(recept_xyz[2 + 3 * receptor] - ref_node_xyz(to_node, 2, node), 2);
    *max_error = MAX(*max_error, sqrt(error));
  }
  ref_free(recept_node);
  ref_free(recept_xyz);

  error = *max_error;
  RSS(ref_mpi_max(ref_mpi, &error, max_error, REF_DBL_TYPE), "max");
  RSS(ref_mpi_bcast(ref_mpi, &max_error, 1, REF_DBL_TYPE), "max");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_stats(REF_INTERP ref_interp) {
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT extrapolate = 0;
  REF_INT node;
  REF_DBL this_bary, max_error, min_bary;

  if (ref_mpi_once(ref_mpi)) {
    if (ref_interp->n_tree > 0)
      printf("tree search: %d found, %.2f avg cells\n", ref_interp->n_tree,
             (REF_DBL)ref_interp->tree_cells / (REF_DBL)ref_interp->n_tree);
    if (ref_interp->n_walk > 0 || ref_interp->n_terminated > 0)
      printf("walks: %d successful, %.2f avg cells, %d terminated\n",
             ref_interp->n_walk,
             (REF_DBL)ref_interp->walk_steps / (REF_DBL)ref_interp->n_walk,
             ref_interp->n_terminated);
    printf("geom nodes: %d failed, %d successful\n", ref_interp->n_geom_fail,
           ref_interp->n_geom);
  }

  each_ref_node_valid_node(to_node, node) if (ref_node_owned(to_node, node)) {
    this_bary = MIN(
        MIN(ref_interp->bary[0 + 4 * node], ref_interp->bary[1 + 4 * node]),
        MIN(ref_interp->bary[2 + 4 * node], ref_interp->bary[3 + 4 * node]));
    if (this_bary < ref_interp->inside) extrapolate++;
  }
  node = extrapolate;
  RSS(ref_mpi_sum(ref_mpi, &node, &extrapolate, 1, REF_INT_TYPE), "sum");

  RSS(ref_interp_max_error(ref_interp, &max_error), "me");
  RSS(ref_interp_min_bary(ref_interp, &min_bary), "mb");

  if (ref_mpi_once(ref_mpi)) {
    printf("interp min bary %e max error %e extrap %d\n", min_bary, max_error,
           extrapolate);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_tec(REF_INTERP ref_interp, const char *filename) {
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_NODE ref_node = ref_grid_node(to_grid);
  FILE *file;
  REF_INT item;

  /* skip if noting to show */
  if (0 == ref_list_n(ref_interp->visualize)) return REF_SUCCESS;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine interp\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(file, "zone t=\"exhaust\", i=%d, datapacking=%s\n",
          ref_list_n(ref_interp->visualize), "point");

  each_ref_list_item(ref_interp->visualize, item) fprintf(
      file, "%.15e %.15e %.15e\n",
      ref_node_xyz(ref_node, 0, ref_list_value(ref_interp->visualize, item)),
      ref_node_xyz(ref_node, 1, ref_list_value(ref_interp->visualize, item)),
      ref_node_xyz(ref_node, 2, ref_list_value(ref_interp->visualize, item)));

  fclose(file);
  return REF_SUCCESS;
}
static int nq = 24;
static double xq[] = {
    -0.570794257481696, -0.287617227554912, -0.570794257481696,
    -0.570794257481696, -0.918652082930777, 0.755956248792332,
    -0.918652082930777, -0.918652082930777, -0.355324219715449,
    -0.934027340853653, -0.355324219715449, -0.355324219715449,
    -0.872677996249965, -0.872677996249965, -0.872677996249965,
    -0.872677996249965, -0.872677996249965, -0.872677996249965,
    -0.460655337083368, -0.460655337083368, -0.460655337083368,
    0.206011329583298,  0.206011329583298,  0.206011329583298};
static double yq[] = {
    -0.570794257481696, -0.570794257481696, -0.287617227554912,
    -0.570794257481696, -0.918652082930777, -0.918652082930777,
    0.755956248792332,  -0.918652082930777, -0.355324219715449,
    -0.355324219715449, -0.934027340853653, -0.355324219715449,
    -0.872677996249965, -0.460655337083368, -0.872677996249965,
    0.206011329583298,  -0.460655337083368, 0.206011329583298,
    -0.872677996249965, -0.872677996249965, 0.206011329583298,
    -0.872677996249965, -0.872677996249965, -0.460655337083368};
static double zq[] = {
    -0.570794257481696, -0.570794257481696, -0.570794257481696,
    -0.287617227554912, -0.918652082930777, -0.918652082930777,
    -0.918652082930777, 0.755956248792332,  -0.355324219715449,
    -0.355324219715449, -0.355324219715449, -0.934027340853653,
    -0.460655337083368, -0.872677996249965, 0.206011329583298,
    -0.872677996249965, 0.206011329583298,  -0.460655337083368,
    -0.872677996249965, 0.206011329583298,  -0.872677996249965,
    -0.460655337083368, -0.872677996249965, -0.872677996249965};
static double wq[] = {
    0.053230333677557, 0.053230333677557, 0.053230333677557, 0.053230333677557,
    0.013436281407094, 0.013436281407094, 0.013436281407094, 0.013436281407094,
    0.073809575391540, 0.073809575391540, 0.073809575391540, 0.073809575391540,
    0.064285714285714, 0.064285714285714, 0.064285714285714, 0.064285714285714,
    0.064285714285714, 0.064285714285714, 0.064285714285714, 0.064285714285714,
    0.064285714285714, 0.064285714285714, 0.064285714285714, 0.064285714285714};
REF_STATUS ref_interp_integrate(REF_GRID ref_grid, REF_DBL *canidate,
                                REF_DBL *truth, REF_INT norm_power,
                                REF_DBL *error) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT i, cell, cell_node, node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL diff, volume, total_volume, bary[4];
  REF_DBL canidate_at_gauss_point, truth_at_gauss_point;
  REF_INT part;
  *error = 0.0;
  total_volume = 0.0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "owner");
    if (part != ref_mpi_rank(ref_mpi)) continue;
    RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");
    total_volume += volume;
    for (i = 0; i < nq; i++) {
      bary[1] = 0.5 * (1.0 + xq[i]);
      bary[2] = 0.5 * (1.0 + yq[i]);
      bary[3] = 0.5 * (1.0 + zq[i]);
      bary[0] = 1.0 - bary[1] - bary[2] - bary[3];
      canidate_at_gauss_point = 0.0;
      truth_at_gauss_point = 0.0;
      each_ref_cell_cell_node(ref_cell, cell_node) {
        node = nodes[cell_node];
        canidate_at_gauss_point += bary[cell_node] * canidate[node];
        truth_at_gauss_point += bary[cell_node] * truth[node];
      }
      diff = ABS(canidate_at_gauss_point - truth_at_gauss_point);
      *error += (6.0 / 8.0) * wq[i] * volume * pow(diff, norm_power);
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, error, 1, REF_DBL_TYPE), "all sum");
  *error = pow(*error, 1.0 / ((REF_DBL)norm_power));
  RSS(ref_mpi_allsum(ref_mpi, &total_volume, 1, REF_DBL_TYPE), "all sum");
  *error /= total_volume;
  return REF_SUCCESS;
}

REF_STATUS ref_interp_convergence_rate(REF_DBL f3, REF_DBL h3, REF_DBL f2,
                                       REF_DBL h2, REF_DBL f1, REF_DBL h1,
                                       REF_DBL *rate) {
  REF_DBL e12, e23, r12, r23;
  REF_DBL beta, omega = 0.5;
  REF_DBL p, last_p;
  REF_INT i;
  /*   AIAA JOURNAL Vol. 36, No. 5, May 1998
       Verification of Codes and Calculations
       Patrick J. Roache */
  /*
  f3 = a coarse grid numerical solution obtained with grid spacing h3
  f2 = a medium grid numerical solution obtained with grid spacing h2
  f1 = a fine grid numerical solution obtained with grid spacing h1
  */
  *rate = -1.0;
  e12 = ABS(f1 - f2);
  e23 = ABS(f2 - f3);
  r12 = h2 / h1;
  r23 = h3 / h2;
  p = log(e23 / e12);
  for (i = 0; i < 20; i++) {
    last_p = p;
    if (!ref_math_divisible(e23, e12)) return REF_SUCCESS;
    beta = ((pow(r12, p) - 1.0) / (pow(r23, p) - 1.0)) * (e23 / e12);
    if (!ref_math_divisible(log(beta), log(r12))) return REF_SUCCESS;
    p = omega * p + (1.0 - omega) * log(beta) / log(r12);
    if (ABS(p - last_p) < 0.0001) break;
  }
  *rate = p;
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_plt_string(FILE *file, char *string, int maxlen) {
  int i, letter;
  for (i = 0; i < maxlen; i++) {
    REIS(1, fread(&letter, sizeof(int), 1, file), "plt string letter");
    string[i] = (char)letter;
    if (0 == letter) {
      return REF_SUCCESS;
    }
  }
  return REF_FAILURE;
}

static REF_STATUS ref_interp_plt_header(FILE *file, REF_INT *nvar,
                                        REF_LIST zone_nnode,
                                        REF_LIST zone_nelem) {
  char header[9];
  int endian, filetype;
  char title[1024], varname[1024], zonename[1024];
  int var, numvar;
  float zonemarker;
  int parent, strand, notused, zonetype, packing, location, neighbor;
  double solutiontime;
  int miscellaneous, i;
  int numpts, numelem;
  int dim, aux;

  RAS(header == fgets(header, 6, file), "header error");
  header[5] = '\0';
  REIS(0, strncmp(header, "#!TDV", 5), "header '#!TDV' missing")
  RAS(header == fgets(header, 4, file), "version error");
  header[4] = '\0';
  REIS(0, strncmp(header, "112", 5), "expected version '112'")

  REIS(1, fread(&endian, sizeof(int), 1, file), "magic");
  REIS(1, endian, "expected little endian plt");
  REIS(1, fread(&filetype, sizeof(int), 1, file), "filetype");
  REIS(0, filetype, "expected full filetype");

  RSS(ref_interp_plt_string(file, title, 1024), "read title");
  printf("plt title '%s'\n", title);

  REIS(1, fread(&numvar, sizeof(int), 1, file), "numvar");
  printf("plt number of variables %d\n", numvar);
  *nvar = numvar;

  for (var = 0; var < numvar; var++) {
    RSS(ref_interp_plt_string(file, varname, 1024), "read variable name");
    printf("plt varable name %d '%s'\n", var, varname);
  }

  REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  while (ABS(299.0 - zonemarker) < 1.0e-7) {
    RSS(ref_interp_plt_string(file, zonename, 1024), "read zonename");
    printf("plt zonename '%s'\n", zonename);

    REIS(1, fread(&parent, sizeof(int), 1, file), "parent");
    REIS(1, fread(&strand, sizeof(int), 1, file), "strand");
    REIS(1, fread(&solutiontime, sizeof(double), 1, file), "solutiontime");
    REIS(1, fread(&notused, sizeof(int), 1, file), "notused");
    REIS(-1, notused, "not unused shoud be -1 plt");
    REIS(1, fread(&zonetype, sizeof(int), 1, file), "zonetype");
    REIS(5, zonetype, "only FEBRICK plt zone implemented");
    REIS(1, fread(&packing, sizeof(int), 1, file), "packing");
    REIS(1, packing, "only point packing plt implemented");
    REIS(1, fread(&location, sizeof(int), 1, file), "location");
    REIS(0, location, "only node data location plt implemented");
    REIS(1, fread(&neighbor, sizeof(int), 1, file), "neighbor");
    REIS(0, neighbor, "no face nieghbor  plt implemented");

    for (i = 0; i < numvar; i++) {
      REIS(1, fread(&miscellaneous, sizeof(int), 1, file), "mystery");
      REIS(0, miscellaneous, "mystery data nonzero plt");
    }

    REIS(1, fread(&numpts, sizeof(int), 1, file), "numpts");
    REIS(1, fread(&numelem, sizeof(int), 1, file), "numelem");

    RSS(ref_list_push(zone_nnode, numpts), "save nnode");
    RSS(ref_list_push(zone_nelem, numelem), "save nelem");

    for (i = 0; i < 3; i++) {
      REIS(1, fread(&dim, sizeof(int), 1, file), "dim");
      REIS(0, dim, "dim nonzero plt");
    }
    REIS(1, fread(&aux, sizeof(int), 1, file), "aux");
    REIS(0, aux, "aux nonzero plt");

    REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  }

  RWDS(357.0, zonemarker, -1.0, "end of header marker expected");

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_plt_data(FILE *file, REF_INT nvar,
                                      REF_LIST zone_nnode, REF_LIST zone_nelem,
                                      REF_INT *length, REF_DBL **soln) {
  float zonemarker;
  int dataformat;
  REF_INT i, node, elem;
  int passive, sharing, conn, c2n;
  double minval, maxval;
  REF_INT nnode, nelem;
  float var;

  RSS(ref_list_shift(zone_nnode, &nnode), "zone node size");
  RSS(ref_list_shift(zone_nelem, &nelem), "zone elem size");
  *length = nnode;

  REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  RWDS(299.0, zonemarker, -1.0, "start of data header expected");

  for (i = 0; i < nvar; i++) {
    REIS(1, fread(&dataformat, sizeof(int), 1, file), "dim");
    REIS(1, dataformat, "dim nonzero plt");
  }
  REIS(1, fread(&passive, sizeof(int), 1, file), "dim");
  if (1 == passive) {
    for (i = 0; i < nvar; i++) {
      REIS(1, fread(&passive, sizeof(int), 1, file), "dim");
      REIS(0, passive, "passive varaible nonzero plt");
    }
  }
  REIS(1, fread(&sharing, sizeof(int), 1, file), "dim");
  if (1 == sharing) {
    for (i = 0; i < nvar; i++) {
      REIS(1, fread(&sharing, sizeof(int), 1, file), "dim");
      REIS(-1, sharing, "variable sharing not implemented plt");
    }
  }
  REIS(1, fread(&conn, sizeof(int), 1, file), "dim");
  REIS(-1, conn, "connectivity sharing not implemented plt");

  for (i = 0; i < nvar; i++) {
    REIS(1, fread(&minval, sizeof(double), 1, file), "dim");
    REIS(1, fread(&maxval, sizeof(double), 1, file), "dim");
  }

  ref_malloc(*soln, nvar * nnode, REF_DBL);

  for (i = 0; i < nvar; i++) {
    for (node = 0; node < nnode; node++) {
      REIS(1, fread(&var, sizeof(float), 1, file), "dim");
      (*soln)[i + nvar * node] = var;
    }
  }

  for (elem = 0; elem < nelem; elem++) {
    for (i = 0; i < 8; i++) {
      REIS(1, fread(&c2n, sizeof(int), 1, file), "dim");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_plt(REF_GRID ref_grid, const char *filename,
                          REF_INT *ldim, REF_DBL **scalar) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  FILE *file = NULL;
  REF_INT nvar;
  REF_LIST zone_nnode, zone_nelem, touching;
  REF_INT zone, nzone, length, node, i, point;
  REF_DBL *soln;
  REF_SEARCH ref_search;
  REF_DBL radius, position[3], dist, best_dist;
  REF_INT best, item;

  RSS(ref_search_create(&ref_search, ref_node_n(ref_node)), "create search");
  each_ref_node_valid_node(ref_node, node) {
    radius = 0.0;
    RSS(ref_search_insert(ref_search, node, ref_node_xyz_ptr(ref_node, node),
                          radius),
        "ins");
  }

  if (ref_mpi_once(ref_mpi)) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    RSS(ref_list_create(&zone_nnode), "nnode list");
    RSS(ref_list_create(&zone_nelem), "nelem list");

    RSS(ref_interp_plt_header(file, &nvar, zone_nnode, zone_nelem),
        "parse header");
    nzone = ref_list_n(zone_nnode);
  }
  RSS(ref_mpi_bcast(ref_mpi, &nvar, 1, REF_INT_TYPE), "b nvar");
  RSS(ref_mpi_bcast(ref_mpi, &nzone, 1, REF_INT_TYPE), "b nzone");

  *ldim = nvar - 3;
  ref_malloc_init(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL, -999.0);

  RSS(ref_list_create(&touching), "tounching list");
  for (zone = 0; zone < nzone; zone++) {
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_interp_plt_data(file, nvar, zone_nnode, zone_nelem, &length,
                              &soln),
          "read data");
      RSS(ref_mpi_bcast(ref_mpi, &length, 1, REF_INT_TYPE), "b length");
      RSS(ref_mpi_bcast(ref_mpi, soln, nvar * length, REF_DBL_TYPE), "b soln");
    } else {
      RSS(ref_mpi_bcast(ref_mpi, &length, 1, REF_INT_TYPE), "b length");
      ref_malloc(soln, nvar * length, REF_DBL);
      RSS(ref_mpi_bcast(ref_mpi, soln, nvar * length, REF_DBL_TYPE), "b soln");
    }
    for (point = 0; point < length; point++) {
      for (i = 0; i < 3; i++) {
        position[i] = soln[i + nvar * point];
      }
      /* single precision */
      radius =
          100.0 * 1.0e-8 *
          sqrt(pow(position[0], 2) + pow(position[1], 2) + pow(position[2], 2));
      RSS(ref_search_touching(ref_search, touching, position, radius),
          "search tree");
      best_dist = 1.0e+200;
      best = REF_EMPTY;
      each_ref_list_item(touching, item) {
        node = ref_list_value(touching, item);
        dist = sqrt(pow(ref_node_xyz(ref_node, 0, node) - position[0], 2) +
                    pow(ref_node_xyz(ref_node, 1, node) - position[1], 2) +
                    pow(ref_node_xyz(ref_node, 2, node) - position[2], 2));
        if (dist < best_dist) {
          best_dist = dist;
          best = node;
        }
      }
      if (REF_EMPTY != best) {
        for (i = 3; i < nvar; i++) {
          (*scalar)[(i - 3) + (*ldim) * best] = soln[i + nvar * point];
        }
      }
      RSS(ref_list_erase(touching), "erase");
    }
    free(soln);
  }
  RSS(ref_list_free(touching), "free touching");

  if (ref_mpi_once(ref_mpi)) {
    fclose(file);
    ref_list_free(zone_nnode);
    ref_list_free(zone_nelem);
  }

  RSS(ref_search_free(ref_search), "free search");

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_from_part_status(REF_INTERP ref_interp,
                                              REF_INT *from_part) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);

  REF_GLOB n_set, n_moving;
  REF_INT node;

  n_set = 0;
  n_moving = 0;
  each_ref_node_valid_node(from_node, node) {
    if (ref_node_owned(from_node, node) && REF_EMPTY != from_part[node]) {
      n_set++;
      if (from_part[node] != ref_node_part(from_node, node)) {
        n_moving++;
      }
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, &n_set, 1, REF_LONG_TYPE), "sum n set");
  RSS(ref_mpi_allsum(ref_mpi, &n_moving, 1, REF_LONG_TYPE), "sum n moving");
  if (ref_mpi_once(ref_mpi) && 0 < ref_node_n_global(from_node) && 0 < n_set) {
    printf(" %6.2f %% " REF_GLOB_FMT " set %6.2f %% " REF_GLOB_FMT
           " moving of " REF_GLOB_FMT " recept nodes\n",
           100.0 * (REF_DBL)n_set / (REF_DBL)ref_node_n_global(from_node),
           n_set, 100.0 * (REF_DBL)n_moving / (REF_DBL)n_set, n_moving,
           ref_node_n_global(from_node));
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_from_part_neighbor(REF_INTERP ref_interp,
                                                REF_INT *from_part,
                                                REF_INT node) {
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_CELL ref_cell = ref_grid_tet(from_grid);
  REF_INT item, cell, cell_node;

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    each_ref_cell_cell_node(ref_cell, cell_node) {
      if (REF_EMPTY != from_part[ref_cell_c2n(ref_cell, cell_node, cell)]) {
        from_part[node] = from_part[ref_cell_c2n(ref_cell, cell_node, cell)];
        return REF_SUCCESS;
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_interp_fill_empty_from_part(REF_INTERP ref_interp,
                                                  REF_INT *from_part) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_BOOL again;
  REF_INT nsweeps, node;
  nsweeps = 0;
  again = REF_TRUE;
  while (again) {
    nsweeps++;
    again = REF_FALSE;
    each_ref_node_valid_node(from_node, node) {
      if (ref_node_owned(from_node, node) && REF_EMPTY == from_part[node]) {
        RSS(ref_interp_from_part_neighbor(ref_interp, from_part, node), "fill");
        again = again || REF_EMPTY != from_part[node];
      }
    }

    RUS(200, nsweeps, "too many sweeps, stop inf loop");
    RSS(ref_mpi_all_or(ref_mpi, &again), "mpi all or");
    if (again) {
      RSS(ref_node_ghost_int(from_node, from_part, 1), "ghost from_part");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_from_part(REF_INTERP ref_interp, REF_INT *to_part) {
  REF_MPI ref_mpi = ref_interp_mpi(ref_interp);
  REF_GRID to_grid = ref_interp_to_grid(ref_interp);
  REF_GRID from_grid = ref_interp_from_grid(ref_interp);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_INT node, i, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *from_part;
  REF_INT recept, n_recept, donation, n_donor;
  REF_INT find, n_find, lookedup, n_lookedup;
  REF_INT *donor_ret, *donor_cell, *donor_donation, *donor_part,
      *donor_origpart;
  REF_INT *recept_part, *recept_ret, *recept_cell;
  REF_GLOB *recept_global, *donor_global, *donor_nodes;
  REF_DBL *recept_bary, *donor_bary;
  REF_INT *find_ret, *find_donation, *find_cell;
  REF_GLOB *find_nodes;
  REF_INT *lookedup_donation, *lookedup_cell;
  REF_DBL max_error;

  if (ref_grid_twod(from_grid)) from_cell = ref_grid_tri(from_grid);

  RSS(ref_interp_max_error(ref_interp, &max_error), "max error");
  if (ref_mpi_once(ref_grid_mpi(to_grid))) {
    printf("starting %e max error\n", max_error);
  }

  if (ref_node_max(to_node) > ref_interp_max(ref_interp)) {
    RSS(ref_interp_resize(ref_interp, ref_node_max(to_node)), "resize");
  }

  ref_malloc_init(from_part, ref_node_max(from_node), REF_INT, REF_EMPTY);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node) && REF_EMPTY != ref_interp->cell[node]) {
      n_recept++;
    }
  }

  ref_malloc(recept_bary, 4 * n_recept, REF_DBL);
  ref_malloc(recept_cell, n_recept, REF_INT);
  ref_malloc(recept_global, n_recept, REF_GLOB);
  ref_malloc(recept_ret, n_recept, REF_INT);
  ref_malloc(recept_part, n_recept, REF_INT);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node) && REF_EMPTY != ref_interp->cell[node]) {
      for (i = 0; i < 4; i++) {
        recept_bary[i + 4 * n_recept] = ref_interp->bary[i + 4 * node];
      }
      recept_cell[n_recept] = ref_interp->cell[node];
      recept_global[n_recept] = ref_node_global(to_node, node);
      recept_part[n_recept] = ref_interp->part[node];
      recept_ret[n_recept] = to_part[node];
      n_recept++;
    }
  }

  RSS(ref_mpi_blindsend(ref_mpi, recept_part, (void *)recept_cell, 1, n_recept,
                        (void **)(&donor_cell), &n_donor, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, recept_part, (void *)recept_ret, 1, n_recept,
                        (void **)(&donor_ret), &n_donor, REF_INT_TYPE),
      "blind send ret");
  RSS(ref_mpi_blindsend(ref_mpi, recept_part, (void *)recept_global, 1,
                        n_recept, (void **)(&donor_global), &n_donor,
                        REF_GLOB_TYPE),
      "blind send global");
  RSS(ref_mpi_blindsend(ref_mpi, recept_part, (void *)recept_bary, 4, n_recept,
                        (void **)(&donor_bary), &n_donor, REF_DBL_TYPE),
      "blind send bary");

  ref_malloc(donor_nodes, 4 * n_donor, REF_GLOB);
  ref_malloc(donor_donation, n_donor, REF_INT);
  ref_malloc(donor_part, n_donor, REF_INT);
  ref_malloc(donor_origpart, n_donor, REF_INT);

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (i = 0; i < ref_cell_node_per(from_cell); i++) {
      from_part[nodes[i]] = donor_ret[donation];
    }
  }
  RSS(ref_node_ghost_int(from_node, from_part, 1), "ghost from_part");
  RSS(ref_interp_from_part_status(ref_interp, from_part), "from part status");
  RSS(ref_interp_fill_empty_from_part(ref_interp, from_part), "fill part");
  RSS(ref_interp_from_part_status(ref_interp, from_part), "from part status");

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (i = 0; i < ref_cell_node_per(from_cell); i++) {
      donor_nodes[i + 4 * donation] = ref_node_global(from_node, nodes[i]);
    }
    donor_donation[donation] = donation;
    RSS(ref_cell_part_cell_node(from_cell, from_node, donor_cell[donation],
                                &cell_node),
        "part cell_node");
    donor_part[donation] = from_part[nodes[cell_node]];
    donor_origpart[donation] = ref_mpi_rank(ref_mpi);
  }

  /* set parts of from_node */
  for (node = 0; node < ref_node_max(from_node); node++)
    ref_node_part(from_node, node) = from_part[node];

  /* shuffle from_node */
  RSS(ref_migrate_shufflin(from_grid), "shufflin from grid");

  /* use new from part to translate nodes to cell (back and forth) */

  RSS(ref_mpi_blindsend(ref_mpi, donor_part, (void *)donor_nodes, 4, n_donor,
                        (void **)(&find_nodes), &n_find, REF_GLOB_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, donor_part, (void *)donor_donation, 1, n_donor,
                        (void **)(&find_donation), &n_find, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, donor_part, (void *)donor_origpart, 1, n_donor,
                        (void **)(&find_ret), &n_find, REF_INT_TYPE),
      "blind send cell");

  ref_malloc(find_cell, n_find, REF_INT);

  for (find = 0; find < n_find; find++) {
    for (i = 0; i < ref_cell_node_per(from_cell); i++) {
      RSS(ref_node_local(from_node, find_nodes[i + 4 * find], &(nodes[i])),
          "g2l");
    }
    RSS(ref_cell_with(from_cell, nodes, &(find_cell[find])),
        "find cell with nodes");
  }

  RSS(ref_mpi_blindsend(ref_mpi, find_ret, (void *)find_cell, 1, n_find,
                        (void **)(&lookedup_cell), &n_lookedup, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, find_ret, (void *)find_donation, 1, n_find,
                        (void **)(&lookedup_donation), &n_lookedup,
                        REF_INT_TYPE),
      "blind send cell");

  for (lookedup = 0; lookedup < n_lookedup; lookedup++) {
    donor_cell[lookedup_donation[lookedup]] = lookedup_cell[lookedup];
  }

  /* shuffle to */
  for (node = 0; node < ref_node_max(to_node); node++) {
    ref_node_part(to_node, node) = to_part[node];
  }

  RSS(ref_migrate_shufflin(to_grid), "shufflin to grid");

  /* return from data to to grid and refill ref_interp->data */
  RSS(ref_interp_reset(ref_interp), "ref_interp resize/reset");

  ref_free(recept_part);
  ref_free(recept_global);
  ref_free(recept_cell);
  ref_free(recept_bary);

  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_cell, 1, n_donor,
                        (void **)(&recept_cell), &n_recept, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_global, 1, n_donor,
                        (void **)(&recept_global), &n_recept, REF_GLOB_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_part, 1, n_donor,
                        (void **)(&recept_part), &n_recept, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_bary, 4, n_donor,
                        (void **)(&recept_bary), &n_recept, REF_DBL_TYPE),
      "blind send cell");

  for (recept = 0; recept < n_recept; recept++) {
    RSS(ref_node_local(to_node, recept_global[recept], &node), "g2l");
    ref_interp->cell[node] = recept_cell[recept];
    ref_interp->part[node] = recept_part[recept];
    for (i = 0; i < 4; i++) {
      ref_interp->bary[i + 4 * node] = recept_bary[i + 4 * recept];
    }
  }

  /* remake interp search tree */
  RSS(ref_search_free(ref_interp_search(ref_interp)), "free search tree");
  RSS(ref_interp_create_search(ref_interp), "(re)build search tree");

  ref_free(lookedup_donation);
  ref_free(lookedup_cell);

  ref_free(find_cell);
  ref_free(find_ret);
  ref_free(find_donation);
  ref_free(find_nodes);

  ref_free(donor_origpart);
  ref_free(donor_part);
  ref_free(donor_donation);
  ref_free(donor_nodes);
  ref_free(donor_bary);
  ref_free(donor_global);
  ref_free(donor_ret);
  ref_free(donor_cell);

  ref_free(recept_part);
  ref_free(recept_ret);
  ref_free(recept_global);
  ref_free(recept_cell);
  ref_free(recept_bary);

  ref_free(from_part);

  RSS(ref_interp_max_error(ref_interp, &max_error), "max error");
  if (ref_mpi_once(ref_grid_mpi(to_grid))) {
    printf("final %e max error\n", max_error);
  }

  return REF_SUCCESS;
}
