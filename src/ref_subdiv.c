
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

#include "ref_subdiv.h"

#include "ref_adj.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

#include "ref_part.h"

static REF_INT ref_subdiv_c2e(REF_SUBDIV ref_subdiv, REF_CELL ref_cell,
                              REF_INT cell_edge, REF_INT cell) {
  REF_INT node0, node1, edge;

  RAE(ref_cell_valid(ref_cell, cell), "invalid cell index");
  RAE(cell_edge < ref_cell_edge_per(ref_cell), "invalid edge index");

  node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
  node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);

  RSE(ref_edge_with(ref_subdiv_edge(ref_subdiv), node0, node1, &edge),
      "look up edge");

  return edge;
}

static REF_INT ref_subdiv_map(REF_SUBDIV ref_subdiv, REF_CELL ref_cell,
                              REF_INT cell) {
  REF_INT edge, map, bit;

  map = 0;
  bit = 1;
  for (edge = 0; edge < ref_cell_edge_per(ref_cell); edge++) {
    map +=
        bit * ref_subdiv_mark(ref_subdiv,
                              ref_subdiv_c2e(ref_subdiv, ref_cell, edge, cell));
    bit *= 2;
  }

  return map;
}

REF_STATUS ref_subdiv_inspect_cell(REF_SUBDIV ref_subdiv, REF_CELL ref_cell,
                                   REF_INT cell) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_INT cell_node, edge;

  for (edge = 0; edge < ref_cell_edge_per(ref_cell); edge++) {
    printf(" %d",
           ref_subdiv_mark(ref_subdiv,
                           ref_subdiv_c2e(ref_subdiv, ref_cell, edge, cell)));
  }
  printf(" cell %d rank %d", cell, ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)));
  each_ref_cell_cell_node(ref_cell, cell_node) {
    printf(" " REF_GLOB_FMT,
           ref_node_global(ref_node, ref_cell_c2n(ref_cell, cell_node, cell)));
  }
  printf("\n");

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_inspect_global(REF_SUBDIV ref_subdiv, REF_INT global0,
                                     REF_INT global1) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_CELL ref_cell;
  REF_INT node0, node1;
  REF_INT edge, part;
  REF_INT group, item, node, cell;

  RXS(ref_node_local(ref_node, global0, &node0), REF_NOT_FOUND, "g2l0");
  RXS(ref_node_local(ref_node, global1, &node1), REF_NOT_FOUND, "g2l1");
  if (REF_EMPTY == node0 || REF_EMPTY == node1) return REF_SUCCESS;
  RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), node0, node1, &edge), "with");
  RSS(ref_edge_part(ref_subdiv_edge(ref_subdiv), edge, &part), "part");

  printf(" mark %d edge %d rank %d part %d global %d %d\n",
         ref_subdiv_mark(ref_subdiv, edge), edge,
         ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)), part, global0, global1);

  each_ref_grid_ref_cell(ref_subdiv_grid(ref_subdiv), group, ref_cell) {
    each_ref_cell_having_node2(ref_cell, node0, node1, item, node, cell) {
      RSS(ref_subdiv_inspect_cell(ref_subdiv, ref_cell, cell), "insp");
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_map_to_edge(REF_INT map) {
  REF_INT edge, bit;

  bit = 2048;
  for (edge = 11; edge >= 0; edge--) {
    if (map >= bit) {
      map -= bit;
      printf("edge %d bit %d\n", edge, bit);
    }
    bit /= 2;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_create(REF_SUBDIV *ref_subdiv_ptr, REF_GRID ref_grid) {
  REF_SUBDIV ref_subdiv;

  ref_malloc(*ref_subdiv_ptr, 1, REF_SUBDIV_STRUCT);

  ref_subdiv = *ref_subdiv_ptr;

  ref_subdiv_grid(ref_subdiv) = ref_grid;

  RSS(ref_edge_create(&(ref_subdiv_edge(ref_subdiv)),
                      ref_subdiv_grid(ref_subdiv)),
      "create edge");

  ref_malloc_init(ref_subdiv->mark, ref_edge_n(ref_subdiv_edge(ref_subdiv)),
                  REF_INT, 0);
  ref_malloc_init(ref_subdiv->node, ref_edge_n(ref_subdiv_edge(ref_subdiv)),
                  REF_INT, REF_EMPTY);

  ref_subdiv->instrument = REF_FALSE;
  ref_subdiv->debug = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_free(REF_SUBDIV ref_subdiv) {
  if (NULL == (void *)ref_subdiv) return REF_NULL;

  ref_free(ref_subdiv->node);
  ref_free(ref_subdiv->mark);
  RSS(ref_edge_free(ref_subdiv_edge(ref_subdiv)), "free edge");

  ref_free(ref_subdiv);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_inspect(REF_SUBDIV ref_subdiv) {
  REF_INT group, cell, cell_edge, edge, map;
  REF_CELL ref_cell;

  each_ref_grid_ref_cell(ref_subdiv_grid(ref_subdiv), group, ref_cell)
      each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    printf(" group %d cell %d map %d\n", group, cell, map);
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      edge = ref_subdiv_c2e(ref_subdiv, ref_cell, cell_edge, cell);
      printf("  edge %d nodes %d %d mark %d\n", edge,
             ref_cell_e2n(ref_cell, 0, cell_edge, cell),
             ref_cell_e2n(ref_cell, 1, cell_edge, cell),
             ref_subdiv_mark(ref_subdiv, edge));
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_n(REF_SUBDIV ref_subdiv, REF_INT *n) {
  REF_INT edge, part;

  *n = 0;
  for (edge = 0; edge < ref_edge_n(ref_subdiv_edge(ref_subdiv)); edge++)
    if (0 != ref_subdiv_mark(ref_subdiv, edge)) {
      RSS(ref_edge_part(ref_subdiv_edge(ref_subdiv), edge, &part), "edge part")
      if (part == ref_mpi_rank(ref_subdiv_mpi(ref_subdiv))) (*n)++;
    }
  RSS(ref_mpi_allsum(ref_subdiv_mpi(ref_subdiv), n, 1, REF_INT_TYPE), "allsum");

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_to_split(REF_SUBDIV ref_subdiv, REF_INT node0,
                                    REF_INT node1) {
  REF_INT edge;

  RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), node0, node1, &edge),
      "missing edge");

  ref_subdiv_mark(ref_subdiv, edge) = 1;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_all(REF_SUBDIV ref_subdiv) {
  REF_INT edge;

  for (edge = 0; edge < ref_edge_n(ref_subdiv_edge(ref_subdiv)); edge++)
    ref_subdiv_mark(ref_subdiv, edge) = 1;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_prism_by_metric(REF_SUBDIV ref_subdiv) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_CELL ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  REF_INT cell, cell_edge;
  REF_INT node0, node1;
  REF_DBL ratio, ratio_limit;

  REF_INT tri_cell_edge;
  REF_INT pri_tri_cell_edge[] = {0, 1, 3, 6, 7, 8};

  ratio_limit = sqrt(2.0);

  each_ref_cell_valid_cell(ref_cell, cell) {
    for (tri_cell_edge = 0; tri_cell_edge < 6; tri_cell_edge++) {
      cell_edge = pri_tri_cell_edge[tri_cell_edge];
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
      if (ratio > ratio_limit)
        RSS(ref_subdiv_mark_to_split(ref_subdiv, node0, node1), "sp");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_prism_by_ratio(REF_SUBDIV ref_subdiv,
                                          REF_DBL *node_ratio) {
  REF_CELL ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  REF_INT cell, cell_edge;
  REF_INT node0, node1;
  REF_DBL ratio_limit;

  REF_INT tri_cell_edge;
  REF_INT pri_tri_cell_edge[] = {0, 1, 3, 6, 7, 8};

  ratio_limit = 0.5 * sqrt(2.0);

  each_ref_cell_valid_cell(ref_cell, cell) {
    for (tri_cell_edge = 0; tri_cell_edge < 6; tri_cell_edge++) {
      cell_edge = pri_tri_cell_edge[tri_cell_edge];
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      if (node_ratio[node0] < ratio_limit && node_ratio[node1] < ratio_limit)
        RSS(ref_subdiv_mark_to_split(ref_subdiv, node0, node1), "sp");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_prism_sides(REF_SUBDIV ref_subdiv) {
  REF_CELL ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  REF_INT cell, cell_edge;
  REF_INT node0, node1;

  REF_INT pri_side_edge;
  REF_INT pri_side_cell_edge[] = {2, 4, 5};

  each_ref_cell_valid_cell(ref_cell, cell) {
    for (pri_side_edge = 0; pri_side_edge < 3; pri_side_edge++) {
      cell_edge = pri_side_cell_edge[pri_side_edge];
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      RSS(ref_subdiv_mark_to_split(ref_subdiv, node0, node1), "sd");
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_node_between(REF_SUBDIV ref_subdiv, REF_INT node0,
                                          REF_INT node1, REF_INT *new_node) {
  REF_INT edge;

  (*new_node) = REF_EMPTY;

  RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), node0, node1, &edge),
      "missing edge");

  (*new_node) = ref_subdiv_node(ref_subdiv, edge);

  return REF_SUCCESS;
}

#define edge_or(ce0, ce1)                                  \
  {                                                        \
    REF_INT ge0, ge1;                                      \
    ge0 = ref_subdiv_c2e(ref_subdiv, ref_cell, ce0, cell); \
    ge1 = ref_subdiv_c2e(ref_subdiv, ref_cell, ce1, cell); \
    if (ref_subdiv_mark(ref_subdiv, ge0) !=                \
        ref_subdiv_mark(ref_subdiv, ge1)) {                \
      again = REF_TRUE;                                    \
      ref_subdiv_mark(ref_subdiv, ge0) = 1;                \
      ref_subdiv_mark(ref_subdiv, ge1) = 1;                \
    }                                                      \
  }

#define promote_2_3(ce0, ce1, ce2)                                             \
  {                                                                            \
    REF_INT ge0, ge1, ge2, sum;                                                \
    ge0 = ref_subdiv_c2e(ref_subdiv, ref_cell, ce0, cell);                     \
    ge1 = ref_subdiv_c2e(ref_subdiv, ref_cell, ce1, cell);                     \
    ge2 = ref_subdiv_c2e(ref_subdiv, ref_cell, ce2, cell);                     \
    sum = ref_subdiv_mark(ref_subdiv, ge0) +                                   \
          ref_subdiv_mark(ref_subdiv, ge1) + ref_subdiv_mark(ref_subdiv, ge2); \
    if (2 == sum) {                                                            \
      again = REF_TRUE;                                                        \
      ref_subdiv_mark(ref_subdiv, ge0) = 1;                                    \
      ref_subdiv_mark(ref_subdiv, ge1) = 1;                                    \
      ref_subdiv_mark(ref_subdiv, ge2) = 1;                                    \
    }                                                                          \
  }

#define promote_2_all()                                                       \
  {                                                                           \
    REF_INT ge0, ge1, ge2, ge3, ge4, ge5, sum;                                \
    ge0 = ref_subdiv_c2e(ref_subdiv, ref_cell, 0, cell);                      \
    ge1 = ref_subdiv_c2e(ref_subdiv, ref_cell, 1, cell);                      \
    ge2 = ref_subdiv_c2e(ref_subdiv, ref_cell, 2, cell);                      \
    ge3 = ref_subdiv_c2e(ref_subdiv, ref_cell, 3, cell);                      \
    ge4 = ref_subdiv_c2e(ref_subdiv, ref_cell, 4, cell);                      \
    ge5 = ref_subdiv_c2e(ref_subdiv, ref_cell, 5, cell);                      \
    sum =                                                                     \
        ref_subdiv_mark(ref_subdiv, ge0) + ref_subdiv_mark(ref_subdiv, ge1) + \
        ref_subdiv_mark(ref_subdiv, ge2) + ref_subdiv_mark(ref_subdiv, ge3) + \
        ref_subdiv_mark(ref_subdiv, ge4) + ref_subdiv_mark(ref_subdiv, ge5);  \
    if (2 == sum)                                                             \
      if ((ge0 > 0 && ge5 > 0) || (ge1 > 0 && ge4 > 0) ||                     \
          (ge2 > 0 && ge1 > 3)) {                                             \
        again = REF_TRUE;                                                     \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 0, cell)) = 1;   \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 1, cell)) = 1;   \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 2, cell)) = 1;   \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 3, cell)) = 1;   \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 4, cell)) = 1;   \
        ref_subdiv_mark(ref_subdiv,                                           \
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 5, cell)) = 1;   \
      }                                                                       \
  }

REF_STATUS ref_subdiv_mark_relax(REF_SUBDIV ref_subdiv) {
  REF_INT group, cell, nsweeps, nmark;
  REF_CELL ref_cell;
  REF_BOOL again;

  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  nsweeps = 0;
  again = REF_TRUE;
  while (again) {
    nsweeps++;
    again = REF_FALSE;

    each_ref_grid_ref_cell(ref_subdiv_grid(ref_subdiv), group, ref_cell)
        each_ref_cell_valid_cell(ref_cell, cell) {
      switch (ref_cell_node_per(ref_cell)) {
        case 4:
          promote_2_3(3, 4, 5);
          promote_2_3(1, 2, 5);
          promote_2_3(0, 2, 4);
          promote_2_3(0, 1, 3);
          promote_2_all();
          break;
        case 5:
          edge_or(0, 7); /* opposite quad edges */
          edge_or(2, 4); /* opposite quad edges */
          break;
        case 6:
          edge_or(0, 6);
          edge_or(3, 8);
          edge_or(1, 7);
          promote_2_3(0, 1, 3);
          promote_2_3(6, 7, 8);
          break;
        default:
          RSS(REF_IMPLEMENT, "implement cell type");
          break;
      }
    }

    RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                           ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
        "ghost mark");

    RSS(ref_mpi_all_or(ref_subdiv_mpi(ref_subdiv), &again), "mpi all or");
  }

  if (ref_subdiv->instrument) {
    RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
    if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
      printf(" %d edges marked after %d relaxations\n", nmark, nsweeps);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_one_of_two(REF_SUBDIV ref_subdiv, REF_INT e0,
                                        REF_INT e1) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_EDGE ref_edge = ref_subdiv_edge(ref_subdiv);
  REF_GLOB e0min, e0max;
  REF_GLOB e1min, e1max;

  /* unmark the edge with the smallest edge globals */
  e0min = MIN(ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, e0)),
              ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, e0)));
  e0max = MAX(ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, e0)),
              ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, e0)));
  e1min = MIN(ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, e1)),
              ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, e1)));
  e1max = MAX(ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, e1)),
              ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, e1)));
  if (e0min == e1min) {
    if (e0max < e1max) {
      ref_subdiv_mark(ref_subdiv, e1) = 0;
    } else {
      ref_subdiv_mark(ref_subdiv, e0) = 0;
    }
  } else {
    if (e0min < e1min) {
      ref_subdiv_mark(ref_subdiv, e1) = 0;
    } else {
      ref_subdiv_mark(ref_subdiv, e0) = 0;
    }
  }

  if (ref_subdiv->debug)
    printf("unmark one of two proc %d edges %d %d\n",
           ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)), e0, e1);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_tet_face(REF_SUBDIV ref_subdiv, REF_CELL ref_cell,
                                      REF_INT cell, REF_BOOL *again, REF_INT s0,
                                      REF_INT s1, REF_INT s2) {
  REF_INT e0, e1, e2;
  e0 = ref_subdiv_c2e(ref_subdiv, ref_cell, s0, cell);
  e1 = ref_subdiv_c2e(ref_subdiv, ref_cell, s1, cell);
  e2 = ref_subdiv_c2e(ref_subdiv, ref_cell, s2, cell);
  if (ref_subdiv_mark(ref_subdiv, e0) && ref_subdiv_mark(ref_subdiv, e1) &&
      !ref_subdiv_mark(ref_subdiv, e2)) {
    ref_subdiv_unmark_one_of_two(ref_subdiv, e0, e1);
    *again = REF_TRUE;
  }
  if (ref_subdiv_mark(ref_subdiv, e0) && !ref_subdiv_mark(ref_subdiv, e1) &&
      ref_subdiv_mark(ref_subdiv, e2)) {
    ref_subdiv_unmark_one_of_two(ref_subdiv, e0, e2);
    *again = REF_TRUE;
  }
  if (!ref_subdiv_mark(ref_subdiv, e0) && ref_subdiv_mark(ref_subdiv, e1) &&
      ref_subdiv_mark(ref_subdiv, e2)) {
    ref_subdiv_unmark_one_of_two(ref_subdiv, e1, e2);
    *again = REF_TRUE;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_tet_opp_edge(REF_SUBDIV ref_subdiv,
                                          REF_CELL ref_cell, REF_INT cell,
                                          REF_BOOL *again, REF_INT s0,
                                          REF_INT s1) {
  REF_INT e0, e1;
  e0 = ref_subdiv_c2e(ref_subdiv, ref_cell, s0, cell);
  e1 = ref_subdiv_c2e(ref_subdiv, ref_cell, s1, cell);
  if (ref_subdiv_mark(ref_subdiv, e0) && ref_subdiv_mark(ref_subdiv, e1)) {
    ref_subdiv_unmark_one_of_two(ref_subdiv, e0, e1);
    if (ref_subdiv->debug)
      printf("unmark tet opp proc %d cell %d sides %d %d edges %d %d\n",
             ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)), cell, s0, s1, e0, e1);
    *again = REF_TRUE;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_tet(REF_SUBDIV ref_subdiv, REF_INT cell,
                                 REF_BOOL *again) {
  REF_CELL ref_cell = ref_grid_tet(ref_subdiv_grid(ref_subdiv));
  REF_INT sum;

  RSS(ref_subdiv_unmark_tet_face(ref_subdiv, ref_cell, cell, again, 3, 4, 5),
      "face 0");
  RSS(ref_subdiv_unmark_tet_face(ref_subdiv, ref_cell, cell, again, 1, 2, 5),
      "face 1");
  RSS(ref_subdiv_unmark_tet_face(ref_subdiv, ref_cell, cell, again, 0, 4, 2),
      "face 2");
  RSS(ref_subdiv_unmark_tet_face(ref_subdiv, ref_cell, cell, again, 0, 1, 3),
      "face 3");

  sum = ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 0, cell)) +
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 1, cell)) +
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 2, cell)) +
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 3, cell)) +
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 4, cell)) +
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, 5, cell));

  if (2 == sum) {
    RSS(ref_subdiv_unmark_tet_opp_edge(ref_subdiv, ref_cell, cell, again, 0, 5),
        "edges 0-5");
    RSS(ref_subdiv_unmark_tet_opp_edge(ref_subdiv, ref_cell, cell, again, 1, 4),
        "edges 1-4");
    RSS(ref_subdiv_unmark_tet_opp_edge(ref_subdiv, ref_cell, cell, again, 2, 3),
        "edges 2-3");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_relax(REF_SUBDIV ref_subdiv) {
  REF_INT group, cell, nsweeps, nmark;
  REF_CELL ref_cell;
  REF_BOOL again;

  /* make sure consistent before starting */
  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  nsweeps = 0;
  again = REF_TRUE;
  while (again) {
    nsweeps++;
    again = REF_FALSE;

    each_ref_grid_ref_cell(ref_subdiv_grid(ref_subdiv), group, ref_cell) {
      each_ref_cell_valid_cell(ref_cell, cell) {
        switch (ref_cell_node_per(ref_cell)) {
          case 4:

            RSS(ref_subdiv_unmark_tet(ref_subdiv, cell, &again), "unmark tet");

            break;
          default:
            /* RSS(REF_IMPLEMENT,"implement cell type"); */
            break;
        }
      }
    }

    /* most conservative, unmark if any ghosts unmarked */
    RSS(ref_edge_ghost_min_int(ref_subdiv_edge(ref_subdiv),
                               ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
        "ghost mark");

    if (nsweeps > 5) {
      RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
      if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
        printf(" %d edges marked after %d unmark relaxations\n", nmark,
               nsweeps);
    }

    RUS(200, nsweeps, "too many sweeps, stop inf loop");

    RSS(ref_mpi_all_or(ref_subdiv_mpi(ref_subdiv), &again), "mpi all or");
  }

  if (ref_subdiv->instrument) {
    RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
    if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
      printf(" %d edges marked after %d unmark relaxations\n", nmark, nsweeps);
  }

  return REF_SUCCESS;
}

#define node_swap(nodes, a, b)   \
  {                              \
    REF_INT nst;                 \
    nst = (nodes)[(a)];          \
    (nodes)[(a)] = (nodes)[(b)]; \
    (nodes)[(b)] = nst;          \
  }

REF_STATUS ref_subdiv_unmark_neg_tet_geom_support(REF_SUBDIV ref_subdiv,
                                                  REF_BOOL *again) {
  REF_GRID ref_grid = ref_subdiv_grid(ref_subdiv);
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_CELL ref_cell = ref_grid_tet(ref_subdiv_grid(ref_subdiv));
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT map;
  REF_DBL xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  REF_DBL *xyzs[4] = {xyz0, xyz1, xyz2, xyz3};
  REF_INT node, i;
  REF_INT edge, split_edge, n0, n1;
  REF_DBL volume;
  REF_BOOL unmark_cell;

  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    unmark_cell = REF_FALSE;
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    switch (map) {
      case 0: /* don't split */
        break;
      case 1:
      case 2:
      case 4:
      case 8:
      case 16:
      case 32:
        split_edge = REF_EMPTY;
        for (edge = 0; edge < ref_cell_edge_per(ref_cell); edge++)
          if (ref_subdiv_mark(ref_subdiv,
                              ref_subdiv_c2e(ref_subdiv, ref_cell, edge, cell)))
            split_edge = edge;
        RAS(REF_EMPTY != split_edge, "edge not found");
        n0 = ref_cell_e2n_gen(ref_cell, 0, split_edge);
        n1 = ref_cell_e2n_gen(ref_cell, 1, split_edge);

        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[n0], nodes[n1], xyzs[n0]),
            "b0");
        RSS(ref_node_xyz_vol(xyzs, &volume), "edge split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;

        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[n0], nodes[n1], xyzs[n1]),
            "b1");
        RSS(ref_node_xyz_vol(xyzs, &volume), "edge split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;

        break;
      case 11:
      case 56:
      case 38:
      case 21:
        /* orient cell for other cases */
        if (56 == map) {
          node_swap(nodes, 0, 3);
          node_swap(nodes, 1, 2);
        }
        if (38 == map) {
          node_swap(nodes, 1, 3);
          node_swap(nodes, 0, 2);
        }
        if (21 == map) {
          node_swap(nodes, 2, 3);
          node_swap(nodes, 0, 1);
        }

        /* near node 0 */
        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[0], nodes[1], xyzs[1]), "s0");
        RSS(ref_geom_xyz_between(ref_grid, nodes[0], nodes[2], xyzs[2]), "s0");
        RSS(ref_node_xyz_vol(xyzs, &volume), "face split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;

        /* near node 1 */
        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[1], nodes[0], xyzs[0]), "s1");
        RSS(ref_geom_xyz_between(ref_grid, nodes[1], nodes[2], xyzs[2]), "s1");
        RSS(ref_node_xyz_vol(xyzs, &volume), "face split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;

        /* near node 2 */
        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[2], nodes[0], xyzs[0]), "s2");
        RSS(ref_geom_xyz_between(ref_grid, nodes[2], nodes[1], xyzs[1]), "s2");
        RSS(ref_node_xyz_vol(xyzs, &volume), "face split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;

        /* center */
        for (node = 0; node < 4; node++) {
          for (i = 0; i < 3; i++) {
            xyzs[node][i] = ref_node_xyz(ref_node, i, nodes[node]);
          }
        }
        RSS(ref_geom_xyz_between(ref_grid, nodes[0], nodes[1], xyzs[0]), "c0");
        RSS(ref_geom_xyz_between(ref_grid, nodes[1], nodes[2], xyzs[1]), "c1");
        RSS(ref_geom_xyz_between(ref_grid, nodes[2], nodes[0], xyzs[2]), "c2");
        RSS(ref_node_xyz_vol(xyzs, &volume), "face split vol");
        if (ref_node_min_volume(ref_node) > volume) unmark_cell = REF_TRUE;
    }
    if (unmark_cell) {
      *again = REF_TRUE;
      for (edge = 0; edge < ref_cell_edge_per(ref_cell); edge++) {
        ref_subdiv_mark(ref_subdiv,
                        ref_subdiv_c2e(ref_subdiv, ref_cell, edge, cell)) = 0;
      }
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_neg_tet_relax(REF_SUBDIV ref_subdiv) {
  REF_INT nsweeps, nmark;
  REF_BOOL again;

  /* make sure consistent before starting */
  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  nsweeps = 0;
  again = REF_TRUE;
  while (again) {
    nsweeps++;
    again = REF_FALSE;

    RSS(ref_subdiv_unmark_neg_tet_geom_support(ref_subdiv, &again), "neg geom");

    /* most conservative, unmark if any ghosts unmarked */
    RSS(ref_edge_ghost_min_int(ref_subdiv_edge(ref_subdiv),
                               ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
        "ghost mark");

    if (nsweeps > 5 || ref_subdiv->instrument) {
      RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
      if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
        printf(" %d edges marked after %d unmark neg geom tet relaxations\n",
               nmark, nsweeps);
    }

    RUS(200, nsweeps, "too many unmark neg geom tet sweeps, stop inf loop");

    RSS(ref_mpi_all_or(ref_subdiv_mpi(ref_subdiv), &again), "mpi all or");
  }

  /* not be required but here for safety? */
  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  if (ref_subdiv->instrument) {
    RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
    if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
      printf(" %d edges marked after %d neg geom tet unmark relaxations\n",
             nmark, nsweeps);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_unmark_geom_support(REF_SUBDIV ref_subdiv) {
  REF_EDGE ref_edge = ref_subdiv_edge(ref_subdiv);
  REF_INT edge;
  REF_BOOL needs_support;

  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    if (ref_subdiv_mark(ref_subdiv, edge)) {
      RSS(ref_geom_support_between(
              ref_subdiv_grid(ref_subdiv), ref_edge_e2n(ref_edge, 0, edge),
              ref_edge_e2n(ref_edge, 1, edge), &needs_support),
          "support check");
      if (needs_support) ref_subdiv_mark(ref_subdiv, edge) = 0;
    }
  }

  /* not be required but here for safety? */
  RSS(ref_edge_ghost_int(ref_subdiv_edge(ref_subdiv),
                         ref_subdiv_mpi(ref_subdiv), ref_subdiv->mark),
      "ghost mark");

  if (ref_subdiv->instrument) {
    REF_INT nmark;
    RSS(ref_subdiv_mark_n(ref_subdiv, &nmark), "count");
    if (ref_mpi_once(ref_subdiv_mpi(ref_subdiv)))
      printf(" %d edges marked without geom support\n", nmark);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_new_node(REF_SUBDIV ref_subdiv) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_EDGE ref_edge = ref_subdiv_edge(ref_subdiv);
  REF_MPI ref_mpi = ref_subdiv_mpi(ref_subdiv);
  REF_INT edge, node, i;
  REF_GLOB global;
  REF_INT part;

  REF_GLOB *edge_global;
  REF_INT *edge_part;
  REF_DBL *edge_real;
  REF_DBL *edge_aux;

  RSS(ref_node_synchronize_globals(ref_node), "sync glob");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    if (ref_subdiv_mark(ref_subdiv, edge)) {
      RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
      if (ref_mpi_rank(ref_mpi) == part) {
        RSS(ref_node_next_global(ref_node, &global), "next global");
        RSS(ref_node_add(ref_node, global, &node), "add node");
        ref_subdiv_node(ref_subdiv, edge) = node;
        RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                      ref_edge_e2n(ref_edge, 1, edge), 0.5,
                                      node),
            "new node");
      }
    }
  }

  RSS(ref_node_shift_new_globals(ref_node), "shift glob");

  ref_malloc_init(edge_global, ref_edge_n(ref_edge), REF_GLOB, REF_EMPTY);
  ref_malloc_init(edge_part, ref_edge_n(ref_edge), REF_INT, REF_EMPTY);
  ref_malloc_init(edge_real, REF_NODE_REAL_PER * ref_edge_n(ref_edge), REF_DBL,
                  -999.0);
  edge_aux = NULL;
  if (ref_node_naux(ref_node) > 0)
    ref_malloc_init(edge_aux, ref_node_naux(ref_node) * ref_edge_n(ref_edge),
                    REF_DBL, 0.0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node = ref_subdiv_node(ref_subdiv, edge);
    if (REF_EMPTY != node) {
      edge_global[edge] = ref_node_global(ref_node, node);
      edge_part[edge] = ref_node_part(ref_node, node);
      for (i = 0; i < REF_NODE_REAL_PER; i++)
        edge_real[i + REF_NODE_REAL_PER * edge] =
            ref_node_real(ref_node, i, node);
      for (i = 0; i < ref_node_naux(ref_node); i++)
        edge_aux[i + ref_node_naux(ref_node) * edge] =
            ref_node_aux(ref_node, i, node);
    }
  }

  RSS(ref_edge_ghost_glob(ref_edge, ref_mpi, edge_global), "global ghost");
  RSS(ref_edge_ghost_int(ref_edge, ref_mpi, edge_part), "part ghost");
  RSS(ref_edge_ghost_dbl(ref_edge, ref_mpi, edge_real, REF_NODE_REAL_PER),
      "xyz ghost");
  if (ref_node_naux(ref_node) > 0)
    RSS(ref_edge_ghost_dbl(ref_edge, ref_mpi, edge_aux,
                           ref_node_naux(ref_node)),
        "aux ghost");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node = ref_subdiv_node(ref_subdiv, edge);
    global = edge_global[edge];
    if (REF_EMPTY == node && REF_EMPTY != global) {
      RSS(ref_node_add(ref_node, global, &node), "add node");
      ref_subdiv_node(ref_subdiv, edge) = node;
      ref_node_part(ref_node, node) = edge_part[edge];
      for (i = 0; i < REF_NODE_REAL_PER; i++)
        ref_node_real(ref_node, i, node) =
            edge_real[i + REF_NODE_REAL_PER * edge];
      for (i = 0; i < ref_node_naux(ref_node); i++)
        ref_node_aux(ref_node, i, node) =
            edge_aux[i + ref_node_naux(ref_node) * edge];
    }
  }

  ref_free(edge_aux);
  ref_free(edge_real);
  ref_free(edge_part);
  ref_free(edge_global);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_add_local_cell(REF_SUBDIV ref_subdiv,
                                            REF_CELL ref_cell, REF_INT *nodes) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_BOOL has_local;
  REF_INT node, new_cell;

  has_local = REF_FALSE;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    has_local = has_local || (ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)) ==
                              ref_node_part(ref_node, nodes[node]));

  if (has_local) RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add");

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_qua(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01, edge12, edge23, edge30;

  ref_cell = ref_grid_qua(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(ref_cell), REF_INT, 0);

  RSS(ref_cell_create(&ref_cell_split, ref_cell_node_per(ref_cell),
                      ref_cell_last_node_is_an_id(ref_cell)),
      "temp cell");

  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[0], nodes[1], &edge01),
        "e01");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[1], nodes[2], &edge12),
        "e12");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[2], nodes[3], &edge23),
        "e23");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[3], nodes[0], &edge30),
        "e30");
    if (ref_subdiv_mark(ref_subdiv, edge01) &&
        ref_subdiv_mark(ref_subdiv, edge23) &&
        ref_subdiv_mark(ref_subdiv, edge12) &&
        ref_subdiv_mark(ref_subdiv, edge30))
      RSS(REF_IMPLEMENT, "all quad edges");
    if (ref_subdiv_mark(ref_subdiv, edge01) !=
            ref_subdiv_mark(ref_subdiv, edge23) ||
        ref_subdiv_mark(ref_subdiv, edge12) !=
            ref_subdiv_mark(ref_subdiv, edge30))
      RSS(REF_IMPLEMENT, "quad edges not paired");

    if (ref_subdiv_mark(ref_subdiv, edge01) &&
        ref_subdiv_mark(ref_subdiv, edge23)) {
      marked_for_removal[cell] = 1;

      RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[3],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

      RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[3],
                                  &(new_nodes[3])),
          "mis");
      RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
    }

    if (ref_subdiv_mark(ref_subdiv, edge12) &&
        ref_subdiv_mark(ref_subdiv, edge30)) {
      marked_for_removal[cell] = 1;

      RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[0],
                                  &(new_nodes[3])),
          "mis");
      RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

      RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
    }
  }

  for (cell = 0; cell < ref_cell_max(ref_cell); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(ref_cell, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(ref_cell_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes), "add local");

  RSS(ref_cell_free(ref_cell_split), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_tri(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL tri, qua;
  REF_CELL tri_split, qua_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01, edge12, edge20;

  tri = ref_grid_tri(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(tri), REF_INT, 0);

  RSS(ref_cell_create(&tri_split, ref_cell_node_per(tri),
                      ref_cell_last_node_is_an_id(tri)),
      "temp tri");

  qua = ref_grid_qua(ref_subdiv_grid(ref_subdiv));
  RSS(ref_cell_create(&qua_split, ref_cell_node_per(qua),
                      ref_cell_last_node_is_an_id(qua)),
      "temp qua");

  each_ref_cell_valid_cell(tri, cell) {
    RSS(ref_cell_nodes(tri, cell, nodes), "nodes");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[0], nodes[1], &edge01),
        "e01");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[1], nodes[2], &edge12),
        "e12");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[2], nodes[0], &edge20),
        "e20");

    if (ref_subdiv_mark(ref_subdiv, edge01) &&
        ref_subdiv_mark(ref_subdiv, edge12) &&
        ref_subdiv_mark(ref_subdiv, edge20)) {
      marked_for_removal[cell] = 1;
      /* near node 0 */
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      /* near node 1 */
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      /* near node 2 */
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      /* center */
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
      continue;
    }

    /*
      2   3-2
      |\  | |
      0-1 0-1
     */

    if (ref_subdiv_mark(ref_subdiv, edge01) &&
        ref_subdiv_mark(ref_subdiv, edge12) &&
        !ref_subdiv_mark(ref_subdiv, edge20)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      new_nodes[4] = new_nodes[3]; /* faceid */
      new_nodes[3] = new_nodes[2]; /* last node */
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[0],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_cell_add(qua_split, new_nodes, &new_cell), "add");
      continue;
    }

    if (ref_subdiv_mark(ref_subdiv, edge12) &&
        ref_subdiv_mark(ref_subdiv, edge20) &&
        !ref_subdiv_mark(ref_subdiv, edge01)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      new_nodes[4] = new_nodes[3]; /* faceid */
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[1],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[3])),
          "mis");
      RSS(ref_cell_add(qua_split, new_nodes, &new_cell), "add");
      continue;
    }

    if (ref_subdiv_mark(ref_subdiv, edge20) &&
        ref_subdiv_mark(ref_subdiv, edge01) &&
        !ref_subdiv_mark(ref_subdiv, edge12)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");

      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      new_nodes[4] = new_nodes[3]; /* faceid */
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                  &(new_nodes[3])),
          "mis");
      RSS(ref_cell_add(qua_split, new_nodes, &new_cell), "add");
      continue;
    }

    if (ref_subdiv_mark(ref_subdiv, edge01)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
    }

    if (ref_subdiv_mark(ref_subdiv, edge12)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
    }

    if (ref_subdiv_mark(ref_subdiv, edge20)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[2])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
      RSS(ref_cell_nodes(tri, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(tri_split, new_nodes, &new_cell), "add");
    }
  }

  for (cell = 0; cell < ref_cell_max(tri); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(tri, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(tri_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, tri, nodes), "add local");

  each_ref_cell_valid_cell_with_nodes(qua_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, qua, nodes), "add local");

  RSS(ref_cell_free(tri_split), "temp tri free");
  RSS(ref_cell_free(qua_split), "temp qua free");

  free(marked_for_removal);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_edg(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL edg;
  REF_CELL edg_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01;

  edg = ref_grid_edg(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(edg), REF_INT, 0);

  RSS(ref_cell_create(&edg_split, ref_cell_node_per(edg),
                      ref_cell_last_node_is_an_id(edg)),
      "temp edg");

  each_ref_cell_valid_cell(edg, cell) {
    RSS(ref_cell_nodes(edg, cell, nodes), "nodes");
    RSS(ref_edge_with(ref_subdiv_edge(ref_subdiv), nodes[0], nodes[1], &edge01),
        "e01");

    if (ref_subdiv_mark(ref_subdiv, edge01)) {
      marked_for_removal[cell] = 1;
      RSS(ref_cell_nodes(edg, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[0])),
          "mis");
      RSS(ref_cell_add(edg_split, new_nodes, &new_cell), "add");
      RSS(ref_cell_nodes(edg, cell, new_nodes), "nodes");
      RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                  &(new_nodes[1])),
          "mis");
      RSS(ref_cell_add(edg_split, new_nodes, &new_cell), "add");
    }
  }

  for (cell = 0; cell < ref_cell_max(edg); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(edg, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(edg_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, edg, nodes), "add local");

  RSS(ref_cell_free(edg_split), "temp edg free");

  free(marked_for_removal);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_pri(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;

  ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(ref_cell), REF_INT, 0);

  RSS(ref_cell_create(&ref_cell_split, ref_cell_node_per(ref_cell),
                      ref_cell_last_node_is_an_id(ref_cell)),
      "temp cell");
  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    switch (map) {
      case 0: /* don't split */
        break;
      case 52: /* prism split top and bottom, edges 2,4,5 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[3],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[4],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[5],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[3],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[5],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
        break;
      case 65: /* prism split edges 0, 6 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
        break;
      case 130: /* prism split edges 1, 7 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[5],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[5],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
        break;
      case 264: /* prism split edges 3, 8 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[5],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[5],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");
        break;
      case 459: /* prism split */
        marked_for_removal[cell] = 1;

        /* near edge 0-3 */
        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[5],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* near edge 1-4 */
        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[0],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[3],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[5],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* near edge 2-5 */
        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[5], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[5], nodes[3],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* center */
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                    &(new_nodes[2])),
            "mis");

        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[5],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[5], nodes[3],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("pri %d, map %d\n", cell, map);
        RSS(REF_IMPLEMENT, "map not implemented yet")
    }
  }

  for (cell = 0; cell < ref_cell_max(ref_cell); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(ref_cell, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(ref_cell_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes), "add local");

  RSS(ref_cell_free(ref_cell_split), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

#define add_cell_with(fnnw0, fnnw1, fnnw2, fnnw3) \
  new_nodes[0] = (fnnw0);                         \
  new_nodes[1] = (fnnw1);                         \
  new_nodes[2] = (fnnw2);                         \
  new_nodes[3] = (fnnw3);                         \
  RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

static REF_STATUS ref_subdiv_split_tet(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;

  REF_INT edge, split_edge, global_edge;

  REF_INT node;

  ref_cell = ref_grid_tet(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(ref_cell), REF_INT, 0);

  RSS(ref_cell_create(&ref_cell_split, ref_cell_node_per(ref_cell),
                      ref_cell_last_node_is_an_id(ref_cell)),
      "temp cell");
  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    switch (map) {
      case 0: /* don't split */
        break;
      case 1:
      case 2:
      case 4:
      case 8:
      case 16:
      case 32:
        split_edge = REF_EMPTY;
        for (edge = 0; edge < ref_cell_edge_per(ref_cell); edge++)
          if (ref_subdiv_mark(ref_subdiv,
                              ref_subdiv_c2e(ref_subdiv, ref_cell, edge, cell)))
            split_edge = edge;
        RAS(REF_EMPTY != split_edge, "edge not found");
        global_edge = ref_subdiv_c2e(ref_subdiv, ref_cell, split_edge, cell);

        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        new_nodes[ref_cell_e2n_gen(ref_cell, 0, split_edge)] =
            ref_subdiv_node(ref_subdiv, global_edge);
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(ref_cell, cell, new_nodes), "nodes");
        new_nodes[ref_cell_e2n_gen(ref_cell, 1, split_edge)] =
            ref_subdiv_node(ref_subdiv, global_edge);
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        break;
      case 11:
      case 56:
      case 38:
      case 21:
        /* orient cell for other cases */

        marked_for_removal[cell] = 1;

        if (56 == map) {
          node_swap(nodes, 0, 3);
          node_swap(nodes, 1, 2);
        }
        if (38 == map) {
          node_swap(nodes, 1, 3);
          node_swap(nodes, 0, 2);
        }
        if (21 == map) {
          node_swap(nodes, 2, 3);
          node_swap(nodes, 0, 1);
        }

        /* near node 0 */
        for (node = 0; node < 4; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* near node 1 */
        for (node = 0; node < 4; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[0],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* near node 2 */
        for (node = 0; node < 4; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        /* center */
        for (node = 0; node < 4; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(ref_cell_split, new_nodes, &new_cell), "add");

        break;
      /*
                          inode3------5------inode2
                             / \              . /
                            /   \          .   /
                           /     \      .     /
                          /       \  .       /
                         /        .\        /
                        2      1    4      3
                       /    .        \    /
                      /  .            \  /
                     /.                \/
                  inode0------0------inode1
       */
      case 63:
        marked_for_removal[cell] = 1;
        {
          REF_INT n0, n1, n2, n3;
          REF_INT e0, e1, e2, e3, e4, e5;
          n0 = nodes[0];
          n1 = nodes[1];
          n2 = nodes[2];
          n3 = nodes[3];
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1], &e0),
              "e");
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2], &e1),
              "e");
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[3], &e2),
              "e");
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2], &e3),
              "e");
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[3], &e4),
              "e");
          RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[3], &e5),
              "e");
          add_cell_with(e0, e2, e1, n0);
          add_cell_with(e0, e3, e4, n1);
          add_cell_with(e1, e5, e3, n2);
          add_cell_with(e2, e4, e5, n3);
          add_cell_with(e0, e5, e1, e2);
          add_cell_with(e0, e5, e2, e4);
          add_cell_with(e0, e5, e4, e3);
          add_cell_with(e0, e5, e3, e1);
        }
        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("tet %d, map %d\n", cell, map);
        RSS(REF_IMPLEMENT, "map not implemented yet")
    }
  }

  for (cell = 0; cell < ref_cell_max(ref_cell); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(ref_cell, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(ref_cell_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes), "add local");

  RSS(ref_cell_free(ref_cell_split), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}
static REF_STATUS ref_subdiv_split_pyr(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL pyr, pri, tet;
  REF_CELL pyr_split, pri_split, tet_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;
  REF_INT node;

  pyr = ref_grid_pyr(ref_subdiv_grid(ref_subdiv));

  ref_malloc_init(marked_for_removal, ref_cell_max(pyr), REF_INT, 0);

  RSS(ref_cell_create(&pyr_split, ref_cell_node_per(pyr),
                      ref_cell_last_node_is_an_id(pyr)),
      "temp pyr");

  pri = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  RSS(ref_cell_create(&pri_split, ref_cell_node_per(pri),
                      ref_cell_last_node_is_an_id(pri)),
      "temp pri");

  tet = ref_grid_tet(ref_subdiv_grid(ref_subdiv));
  RSS(ref_cell_create(&tet_split, ref_cell_node_per(tet),
                      ref_cell_last_node_is_an_id(tet)),
      "temp pri");

  each_ref_cell_valid_cell(pyr, cell) {
    map = ref_subdiv_map(ref_subdiv, pyr, cell);
    RSS(ref_cell_nodes(pyr, cell, nodes), "nodes");
    switch (map) {
      case 0: /* don't split */
        break;
      case 129: /* split into two pyr on edges 0, 7*/
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");
        break;
      case 20: /* split into two pyr on edges 2, 4*/
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[3],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[4],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[3],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");
        break;
      case 72: /* split into pyr and pri edge 3-6 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(pri_split, new_nodes, &new_cell), "add");
        break;
      case 34: /* split into pyr and pri edge 1-5 */
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[2],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(pri_split, new_nodes, &new_cell), "add");
        break;
      case 235: /* split into 1 pyr, 2 pri*/
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[0],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[2], nodes[3],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(pri_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[2],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(pri_split, new_nodes, &new_cell), "add");

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[5])),
            "mis");
        RSS(ref_cell_add(pri_split, new_nodes, &new_cell), "add");

        break;
      case 64: /* split into 1 pyr, 3 tet, on edge 6, nodes 2, 4*/
        marked_for_removal[cell] = 1;

        RSS(ref_cell_nodes(pyr, cell, new_nodes), "nodes");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        new_nodes[0] = nodes[0];
        new_nodes[1] = nodes[1];
        new_nodes[2] = nodes[2];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        new_nodes[0] = nodes[0];
        new_nodes[1] = nodes[2];
        new_nodes[2] = nodes[3];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[4], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        break;
      case 139:
      case 225:
      case 92:

        if (225 == map) {
          for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
          nodes[0] = new_nodes[4];
          nodes[1] = new_nodes[3];
          nodes[2] = new_nodes[2];
          nodes[3] = new_nodes[1];
          nodes[4] = new_nodes[0];
        }

        if (92 == map) {
          for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
          nodes[0] = new_nodes[1];
          nodes[1] = new_nodes[4];
          nodes[2] = new_nodes[2];
          nodes[3] = new_nodes[0];
          nodes[4] = new_nodes[3];
        }

        marked_for_removal[cell] = 1;

        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[4])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_cell_add(pyr_split, new_nodes, &new_cell), "add");

        /* top sides */
        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        new_nodes[0] = new_nodes[4];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[0])),
            "mis");
        new_nodes[1] = new_nodes[3];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        /* center */
        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        for (node = 0; node < 5; node++) new_nodes[node] = nodes[node];
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[2],
                                    &(new_nodes[0])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[0], nodes[1],
                                    &(new_nodes[1])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[1], nodes[2],
                                    &(new_nodes[2])),
            "mis");
        RSS(ref_subdiv_node_between(ref_subdiv, nodes[3], nodes[4],
                                    &(new_nodes[3])),
            "mis");
        RSS(ref_cell_add(tet_split, new_nodes, &new_cell), "add");

        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("pyr %d, map %d\n", cell, map);
        RSS(REF_IMPLEMENT, "map not implemented yet")
    }
  }

  for (cell = 0; cell < ref_cell_max(pyr); cell++)
    if (1 == marked_for_removal[cell])
      RSS(ref_cell_remove(pyr, cell), "remove");

  each_ref_cell_valid_cell_with_nodes(pyr_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, pyr, nodes), "add local");

  each_ref_cell_valid_cell_with_nodes(pri_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, pri, nodes), "add local");

  each_ref_cell_valid_cell_with_nodes(tet_split, cell, nodes)
      RSS(ref_subdiv_add_local_cell(ref_subdiv, tet, nodes), "add local");

  RSS(ref_cell_free(pyr_split), "temp pyr free");
  RSS(ref_cell_free(pri_split), "temp pri free");
  RSS(ref_cell_free(tet_split), "temp pri free");

  free(marked_for_removal);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split(REF_SUBDIV ref_subdiv) {
  REF_GRID ref_grid = ref_subdiv_grid(ref_subdiv);
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_INT node;

  RSS(ref_node_synchronize_globals(ref_node), "sync glob for mark relax");

  RSS(ref_subdiv_unmark_neg_tet_relax(ref_subdiv), "geom neg marks");
  RSS(ref_subdiv_unmark_geom_support(ref_subdiv), "geom marks");
  RSS(ref_subdiv_mark_relax(ref_subdiv), "relax marks");
  RSS(ref_subdiv_unmark_neg_tet_relax(ref_subdiv), "geom neg marks");
  RSS(ref_subdiv_unmark_geom_support(ref_subdiv), "geom marks");
  RSS(ref_subdiv_unmark_relax(ref_subdiv), "relax marks");

  RSS(ref_subdiv_test_impossible_marks(ref_subdiv), "possible");

  RSS(ref_subdiv_new_node(ref_subdiv), "new nodes");

  RSS(ref_subdiv_split_tet(ref_subdiv), "split tet");
  RSS(ref_subdiv_split_pri(ref_subdiv), "split pri");
  /* pyr comes last, it can make other elements too */
  RSS(ref_subdiv_split_pyr(ref_subdiv), "split pyr");

  RSS(ref_subdiv_split_qua(ref_subdiv), "split qua");
  /* tri comes last, it can make qua elements too */
  RSS(ref_subdiv_split_tri(ref_subdiv), "split tri");

  RSS(ref_subdiv_split_edg(ref_subdiv), "split edg");

  /* remove unused nodes on partition boundaries */
  each_ref_node_valid_node(
      ref_node,
      node) if (ref_adj_empty(ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
                ref_adj_empty(ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
                ref_adj_empty(ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
                ref_adj_empty(ref_cell_adj(ref_grid_hex(ref_grid)), node)) {
    if (ref_mpi_rank(ref_subdiv_mpi(ref_subdiv)) ==
        ref_node_part(ref_node, node))
      RSS(REF_FAILURE, "unused local node");
    if (!ref_adj_empty(ref_cell_adj(ref_grid_tri(ref_grid)), node) ||
        !ref_adj_empty(ref_cell_adj(ref_grid_qua(ref_grid)), node))
      RSS(REF_FAILURE, "boundary face node not in vol cells");
    RSS(ref_node_remove_without_global(ref_node, node), "rm");
    RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), node), "rm");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_tec_zone(REF_SUBDIV ref_subdiv, REF_CELL ref_cell,
                                      REF_INT cell, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_EDGE ref_edge = ref_subdiv_edge(ref_subdiv);
  REF_INT cell_edge, edge, node;

  fprintf(file,
          "zone t=\"celledge\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          2 * ref_cell_edge_per(ref_cell), ref_cell_edge_per(ref_cell), "point",
          "felineseg");

  for (cell_edge = 0; cell_edge < ref_cell_edge_per(ref_cell); cell_edge++) {
    edge = ref_subdiv_c2e(ref_subdiv, ref_cell, cell_edge, cell);
    node = ref_edge_e2n(ref_edge, 0, edge);
    fprintf(file, " %.16e %.16e %.16e %d\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            ref_subdiv_mark(ref_subdiv, edge));
    node = ref_edge_e2n(ref_edge, 1, edge);
    fprintf(file, " %.16e %.16e %.16e %d\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            ref_subdiv_mark(ref_subdiv, edge));
  }

  for (cell_edge = 0; cell_edge < ref_cell_edge_per(ref_cell); cell_edge++)
    fprintf(file, " %d %d\n", 1 + 2 * cell_edge, 2 + 2 * cell_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_verify(REF_SUBDIV ref_subdiv) {
  REF_GRID ref_grid = ref_subdiv_grid(ref_subdiv);
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT map;

  FILE *file;

  file = fopen("ref_subdiv_verify.tec", "w");
  if (NULL == (void *)file)
    printf("unable to open %s\n", "ref_subdiv_verify.tec");
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"s\"\n");

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    switch (map) {
      case 0:
      case 1:
      case 2:
      case 4:
      case 8:
      case 16:
      case 32:
      case 11:
      case 56:
      case 38:
      case 21:
      case 63:
        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("tet %d, map %d\n", cell, map);
        RSS(ref_subdiv_tec_zone(ref_subdiv, ref_cell, cell, file), "zone");
    }
  }

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    switch (map) {
      case 0:
      case 129:
      case 20:
      case 72:
      case 34:
      case 235:
      case 139:
      case 225:
      case 92:
        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("pyr %d, map %d\n", cell, map);
        RSS(ref_subdiv_tec_zone(ref_subdiv, ref_cell, cell, file), "zone");
    }
  }

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    switch (map) {
      case 0:
      case 65:
      case 130:
      case 264:
      case 459:
        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("pri %d, map %d\n", cell, map);
        RSS(ref_subdiv_tec_zone(ref_subdiv, ref_cell, cell, file), "zone");
    }
  }

  fclose(file);

  return REF_SUCCESS;
}

#define fill_pri_xyz(ref_node, nodes, xyz)           \
  (xyz)[0][0] = ref_node_xyz(ref_node, 0, nodes[0]); \
  (xyz)[0][1] = ref_node_xyz(ref_node, 1, nodes[0]); \
  (xyz)[0][2] = ref_node_xyz(ref_node, 2, nodes[0]); \
  (xyz)[1][0] = ref_node_xyz(ref_node, 0, nodes[1]); \
  (xyz)[1][1] = ref_node_xyz(ref_node, 1, nodes[1]); \
  (xyz)[1][2] = ref_node_xyz(ref_node, 2, nodes[1]); \
  (xyz)[2][0] = ref_node_xyz(ref_node, 0, nodes[2]); \
  (xyz)[2][1] = ref_node_xyz(ref_node, 1, nodes[2]); \
  (xyz)[2][2] = ref_node_xyz(ref_node, 2, nodes[2]); \
  (xyz)[3][0] = ref_node_xyz(ref_node, 0, nodes[3]); \
  (xyz)[3][1] = ref_node_xyz(ref_node, 1, nodes[3]); \
  (xyz)[3][2] = ref_node_xyz(ref_node, 2, nodes[3]); \
  (xyz)[4][0] = ref_node_xyz(ref_node, 0, nodes[4]); \
  (xyz)[4][1] = ref_node_xyz(ref_node, 1, nodes[4]); \
  (xyz)[4][2] = ref_node_xyz(ref_node, 2, nodes[4]); \
  (xyz)[5][0] = ref_node_xyz(ref_node, 0, nodes[5]); \
  (xyz)[5][1] = ref_node_xyz(ref_node, 1, nodes[5]); \
  (xyz)[5][2] = ref_node_xyz(ref_node, 2, nodes[5]);

#define replace_xyz0_avg(xyz, n0, n1)                 \
  (xyz)[n0][0] = 0.5 * ((xyz)[n0][0] + (xyz)[n1][0]); \
  (xyz)[n0][1] = 0.5 * ((xyz)[n0][1] + (xyz)[n1][1]); \
  (xyz)[n0][2] = 0.5 * ((xyz)[n0][2] + (xyz)[n1][2]);

#define fill_xyz_avg(xyz, n, n0, n1)                          \
  (xyz)[n][0] = 0.5 * (ref_node_xyz(ref_node, 0, nodes[n0]) + \
                       ref_node_xyz(ref_node, 0, nodes[n1])); \
  (xyz)[n][1] = 0.5 * (ref_node_xyz(ref_node, 1, nodes[n0]) + \
                       ref_node_xyz(ref_node, 1, nodes[n1])); \
  (xyz)[n][2] = 0.5 * (ref_node_xyz(ref_node, 2, nodes[n0]) + \
                       ref_node_xyz(ref_node, 2, nodes[n1]));

static REF_STATUS ref_subdiv_test_pri(REF_DBL xyz[6][3], REF_BOOL *possible) {
  REF_INT n1, n2, n3;
  REF_DBL xnorm, ynorm, znorm;
  REF_DBL dx, dy, dz;
  REF_DBL crdot;

  REF_INT i;

  *possible = REF_TRUE;

  n1 = 0;
  n2 = 1;
  n3 = 2;

  xnorm = 0.5 * ((xyz[n2][1] - xyz[n1][1]) * (xyz[n3][2] - xyz[n1][2]) -
                 (xyz[n2][2] - xyz[n1][2]) * (xyz[n3][1] - xyz[n1][1]));
  ynorm = -0.5 * ((xyz[n2][0] - xyz[n1][0]) * (xyz[n3][2] - xyz[n1][2]) -
                  (xyz[n2][2] - xyz[n1][2]) * (xyz[n3][0] - xyz[n1][0]));
  znorm = 0.5 * ((xyz[n2][0] - xyz[n1][0]) * (xyz[n3][1] - xyz[n1][1]) -
                 (xyz[n2][1] - xyz[n1][1]) * (xyz[n3][0] - xyz[n1][0]));

  dx = (xyz[3][0] + xyz[4][0] + xyz[5][0]) / 3.0 -
       (xyz[0][0] + xyz[1][0] + xyz[2][0]) / 3.0;
  dy = (xyz[3][1] + xyz[4][1] + xyz[5][1]) / 3.0 -
       (xyz[0][1] + xyz[1][1] + xyz[2][1]) / 3.0;
  dz = (xyz[3][2] + xyz[4][2] + xyz[5][2]) / 3.0 -
       (xyz[0][2] + xyz[1][2] + xyz[2][2]) / 3.0;

  crdot = dx * xnorm + dy * ynorm + dz * znorm;

  if (crdot < 0.0) {
    printf("bad pri base\n");
    *possible = REF_FALSE;
  }

  n1 = 3;
  n2 = 5;
  n3 = 4;

  xnorm = 0.5 * ((xyz[n2][1] - xyz[n1][1]) * (xyz[n3][2] - xyz[n1][2]) -
                 (xyz[n2][2] - xyz[n1][2]) * (xyz[n3][1] - xyz[n1][1]));
  ynorm = -0.5 * ((xyz[n2][0] - xyz[n1][0]) * (xyz[n3][2] - xyz[n1][2]) -
                  (xyz[n2][2] - xyz[n1][2]) * (xyz[n3][0] - xyz[n1][0]));
  znorm = 0.5 * ((xyz[n2][0] - xyz[n1][0]) * (xyz[n3][1] - xyz[n1][1]) -
                 (xyz[n2][1] - xyz[n1][1]) * (xyz[n3][0] - xyz[n1][0]));

  dx = -(xyz[3][0] + xyz[4][0] + xyz[5][0]) / 3.0 +
       (xyz[0][0] + xyz[1][0] + xyz[2][0]) / 3.0;
  dy = -(xyz[3][1] + xyz[4][1] + xyz[5][1]) / 3.0 +
       (xyz[0][1] + xyz[1][1] + xyz[2][1]) / 3.0;
  dz = -(xyz[3][2] + xyz[4][2] + xyz[5][2]) / 3.0 +
       (xyz[0][2] + xyz[1][2] + xyz[2][2]) / 3.0;

  crdot = dx * xnorm + dy * ynorm + dz * znorm;

  if (crdot < 0.0) {
    printf("bad pri top\n");
    *possible = REF_FALSE;
  }

  if (!*possible) {
    printf(
        "zone t=\"bad\", nodes=6, elements=1, datapacking=point, "
        "zonetype=febrick");
    for (i = 0; i < 6; i++)
      printf(" %.16e %.16e %.16e\n", xyz[i][0], xyz[i][1], xyz[i][2]);
    printf(" 1 2 3 3 4 5 6 6\n");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_test_impossible_marks(REF_SUBDIV ref_subdiv) {
  REF_INT cell;
  REF_CELL ref_cell;
  REF_NODE ref_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[6][3];
  REF_INT map;
  REF_BOOL possible;
  ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));

  ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));

  each_ref_cell_valid_cell(ref_cell, cell) {
    map = ref_subdiv_map(ref_subdiv, ref_cell, cell);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    switch (map) {
      case 0: /* don't split */
        break;
      case 52: /* prism split top and bottom, edges 2,4,5 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 0, 3);
        replace_xyz0_avg(xyz, 1, 4);
        replace_xyz0_avg(xyz, 2, 5);

        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 3, 0);
        replace_xyz0_avg(xyz, 4, 1);
        replace_xyz0_avg(xyz, 5, 2);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        break;
      case 65: /* prism split edges 0, 6 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 0, 1);
        replace_xyz0_avg(xyz, 3, 4);

        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 1, 0);
        replace_xyz0_avg(xyz, 4, 3);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        break;
      case 130: /* prism split edges 1, 7 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 0, 2);
        replace_xyz0_avg(xyz, 3, 5);

        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 2, 0);
        replace_xyz0_avg(xyz, 5, 3);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        break;
      case 264: /* prism split edges 3, 8 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 1, 2);
        replace_xyz0_avg(xyz, 4, 5);

        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 2, 1);
        replace_xyz0_avg(xyz, 5, 4);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        break;
      case 459: /* prism split */
        /* near edge 0-3 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 1, 0);
        replace_xyz0_avg(xyz, 2, 0);
        replace_xyz0_avg(xyz, 4, 3);
        replace_xyz0_avg(xyz, 5, 3);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        /* near edge 1-4 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 0, 1);
        replace_xyz0_avg(xyz, 2, 1);
        replace_xyz0_avg(xyz, 3, 4);
        replace_xyz0_avg(xyz, 5, 4);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        /* near edge 2-5 */
        fill_pri_xyz(ref_node, nodes, xyz);
        replace_xyz0_avg(xyz, 0, 2);
        replace_xyz0_avg(xyz, 1, 2);
        replace_xyz0_avg(xyz, 3, 5);
        replace_xyz0_avg(xyz, 4, 5);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        /* center */
        fill_xyz_avg(xyz, 0, 0, 1);
        fill_xyz_avg(xyz, 1, 1, 2);
        fill_xyz_avg(xyz, 2, 2, 0);

        fill_xyz_avg(xyz, 3, 3, 4);
        fill_xyz_avg(xyz, 4, 4, 5);
        fill_xyz_avg(xyz, 5, 5, 0);
        RSS(ref_subdiv_test_pri(xyz, &possible), "test pri");

        break;
      default:
        RSS(ref_subdiv_map_to_edge(map), "map2edge");
        printf("pri %d, map %d\n", cell, map);
        RSS(REF_IMPLEMENT, "map not implemented yet")
    }
  }

  return REF_SUCCESS;
}
