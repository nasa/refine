
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

#include <stdio.h>
#include <stdlib.h>

#include "ref_swap.h"

/* parallel requirement, all local */
REF_STATUS ref_swap_remove_two_face_cell(REF_GRID ref_grid, REF_INT cell) {
  REF_CELL ref_cell;
  REF_INT cell_face;
  REF_INT node;
  REF_INT cell_nodes[4];
  REF_INT face_nodes[4];
  REF_INT found;
  REF_INT face0, face1;
  REF_INT cell_face0, cell_face1;
  REF_INT faceid0, faceid1;
  REF_INT temp, face;

  ref_cell = ref_grid_tet(ref_grid);

  face0 = REF_EMPTY;
  face1 = REF_EMPTY;
  faceid0 = REF_EMPTY;
  faceid1 = REF_EMPTY;
  cell_face0 = REF_EMPTY;
  cell_face1 = REF_EMPTY;
  for (cell_face = 0; cell_face < 4; cell_face++) {
    for (node = 0; node < 4; node++)
      cell_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
    if (REF_SUCCESS ==
        ref_cell_with(ref_grid_tri(ref_grid), cell_nodes, &found)) {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid), found, face_nodes), "tri");
      if (REF_EMPTY == face0) {
        face0 = found;
        faceid0 = face_nodes[3];
        cell_face0 = cell_face;
      } else {
        if (REF_EMPTY != face1) {
          RSS(REF_INVALID, "three or more faces detected");
        }
        face1 = found;
        faceid1 = face_nodes[3];
        cell_face1 = cell_face;
      }
    }
  }

  if (REF_EMPTY == face0) return REF_INVALID;
  if (REF_EMPTY == face1) return REF_INVALID;
  if (faceid0 != faceid1) return REF_INVALID;

  RSS(ref_cell_remove(ref_grid_tri(ref_grid), face0), "remove tri0");
  RSS(ref_cell_remove(ref_grid_tri(ref_grid), face1), "remove tri1");

  for (cell_face = 0; cell_face < 4; cell_face++) {
    for (node = 0; node < 4; node++)
      cell_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
    if (cell_face != cell_face0 && cell_face != cell_face1) {
      cell_nodes[3] = faceid0;
      temp = cell_nodes[0];
      cell_nodes[0] = cell_nodes[1];
      cell_nodes[1] = temp;
      RSS(ref_cell_add(ref_grid_tri(ref_grid), cell_nodes, &face), "add tri");
    }
  }

  RSS(ref_cell_remove(ref_cell, cell), "remove tet");

  return REF_SUCCESS;
}

REF_STATUS ref_swap_remove_three_face_cell(REF_GRID ref_grid, REF_INT cell) {
  REF_CELL ref_cell;
  REF_INT cell_face;
  REF_INT node;
  REF_INT cell_nodes[4];
  REF_INT cell_face_nodes[4];
  REF_INT face_nodes[4];
  REF_INT found;
  REF_INT face0, face1, face2;
  REF_INT cell_face0, cell_face1, cell_face2;
  REF_INT faceid0, faceid1, faceid2;
  REF_INT temp, face;
  REF_INT remove_this_node;

  ref_cell = ref_grid_tet(ref_grid);

  face0 = REF_EMPTY;
  face1 = REF_EMPTY;
  face2 = REF_EMPTY;
  faceid0 = REF_EMPTY;
  faceid1 = REF_EMPTY;
  faceid2 = REF_EMPTY;
  cell_face0 = REF_EMPTY;
  cell_face1 = REF_EMPTY;
  cell_face2 = REF_EMPTY;
  for (cell_face = 0; cell_face < 4; cell_face++) {
    for (node = 0; node < 4; node++)
      cell_face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
    if (REF_SUCCESS ==
        ref_cell_with(ref_grid_tri(ref_grid), cell_face_nodes, &found)) {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid), found, face_nodes), "tri");
      if (REF_EMPTY == face0) {
        face0 = found;
        faceid0 = face_nodes[3];
        cell_face0 = cell_face;
      } else if (REF_EMPTY == face1) {
        face1 = found;
        faceid1 = face_nodes[3];
        cell_face1 = cell_face;
      } else {
        if (REF_EMPTY != face2) {
          RSS(REF_INVALID, "four faces detected");
        }
        face2 = found;
        faceid2 = face_nodes[3];
        cell_face2 = cell_face;
      }
    }
  }

  if (REF_EMPTY == face0) return REF_INVALID;
  if (REF_EMPTY == face1) return REF_INVALID;
  if (REF_EMPTY == face2) return REF_INVALID;
  if (faceid0 != faceid1 || faceid0 != faceid2) return REF_INVALID;

  RSS(ref_cell_remove(ref_grid_tri(ref_grid), face0), "remove tri0");
  RSS(ref_cell_remove(ref_grid_tri(ref_grid), face1), "remove tri1");
  RSS(ref_cell_remove(ref_grid_tri(ref_grid), face2), "remove tri1");

  cell_face = 0 + 1 + 2 + 3 - cell_face0 - cell_face1 - cell_face2;

  for (node = 0; node < 4; node++)
    face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);

  face_nodes[3] = faceid0;
  temp = face_nodes[0];
  face_nodes[0] = face_nodes[1];
  face_nodes[1] = temp;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), face_nodes, &face), "add tri");

  RSS(ref_cell_nodes(ref_cell, cell, cell_nodes), "tet");

  remove_this_node = cell_nodes[0] + cell_nodes[1] + cell_nodes[2] +
                     cell_nodes[3] - face_nodes[0] - face_nodes[1] -
                     face_nodes[2];

  RSS(ref_cell_remove(ref_cell, cell), "remove tet");

  RSS(ref_node_remove(ref_grid_node(ref_grid), remove_this_node),
      "remove node");

  return REF_SUCCESS;
}

REF_STATUS ref_swap_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL tet = ref_grid_tet(ref_grid);
  REF_INT tri_index, tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT i, face_nodes[4];
  REF_INT tet0, tet1;
  REF_INT faceid0, faceid1;
  REF_BOOL more_than_two;
  REF_INT cell_face, node, cell_nodes[4];
  REF_INT found, found_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT rank;

  each_ref_cell_valid_cell_with_nodes(tri, tri_index, tri_nodes) {
    for (i = 0; i < 3; i++) face_nodes[i] = tri_nodes[i];
    face_nodes[3] = face_nodes[0];
    RSS(ref_cell_with_face(tet, face_nodes, &tet0, &tet1),
        "unable to find tets with face");
    if (REF_EMPTY == tet0) THROW("boundry tet missing");
    if (REF_EMPTY != tet1) THROW("boundry tri has two tets, not manifold");

    /* must be local to swap */
    rank = ref_mpi_rank(ref_grid_mpi(ref_grid));
    if (rank != ref_node_part(ref_node, ref_cell_c2n(tet, 0, tet0)) ||
        rank != ref_node_part(ref_node, ref_cell_c2n(tet, 1, tet0)) ||
        rank != ref_node_part(ref_node, ref_cell_c2n(tet, 2, tet0)) ||
        rank != ref_node_part(ref_node, ref_cell_c2n(tet, 3, tet0)))
      continue;

    faceid0 = REF_EMPTY;
    faceid1 = REF_EMPTY;
    more_than_two = REF_FALSE;
    for (cell_face = 0; cell_face < 4; cell_face++) {
      for (node = 0; node < 4; node++)
        cell_nodes[node] = ref_cell_f2n(tet, node, cell_face, tet0);
      if (REF_SUCCESS == ref_cell_with(tri, cell_nodes, &found)) {
        RSS(ref_cell_nodes(tri, found, found_nodes), "tri");
        if (REF_EMPTY == faceid0) {
          faceid0 = found_nodes[3];
        } else {
          if (REF_EMPTY != faceid1) {
            more_than_two = REF_TRUE;
          }
          faceid1 = found_nodes[3];
        }
      }
    }
    if (REF_EMPTY == faceid0 || REF_EMPTY == faceid1 || faceid0 != faceid1 ||
        more_than_two)
      continue;
    RSS(ref_swap_remove_two_face_cell(ref_grid, tet0), "remove it");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_swap_node23(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                           REF_INT *node2, REF_INT *node3) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT ncell, cell_to_swap[2];
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  *node2 = REF_EMPTY;
  *node3 = REF_EMPTY;
  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cell_to_swap),
      "more then two");
  REIS(2, ncell, "there should be two triangles for manifold");

  RSS(ref_cell_nodes(ref_cell, cell_to_swap[0], nodes), "nodes tri0");
  if (node0 == nodes[0] && node1 == nodes[1]) *node2 = nodes[2];
  if (node0 == nodes[1] && node1 == nodes[2]) *node2 = nodes[0];
  if (node0 == nodes[2] && node1 == nodes[0]) *node2 = nodes[1];
  if (node1 == nodes[0] && node0 == nodes[1]) *node3 = nodes[2];
  if (node1 == nodes[1] && node0 == nodes[2]) *node3 = nodes[0];
  if (node1 == nodes[2] && node0 == nodes[0]) *node3 = nodes[1];
  RSS(ref_cell_nodes(ref_cell, cell_to_swap[1], nodes), "nodes tri0");
  if (node0 == nodes[0] && node1 == nodes[1]) *node2 = nodes[2];
  if (node0 == nodes[1] && node1 == nodes[2]) *node2 = nodes[0];
  if (node0 == nodes[2] && node1 == nodes[0]) *node2 = nodes[1];
  if (node1 == nodes[0] && node0 == nodes[1]) *node3 = nodes[2];
  if (node1 == nodes[1] && node0 == nodes[2]) *node3 = nodes[0];
  if (node1 == nodes[2] && node0 == nodes[0]) *node3 = nodes[1];

  RUS(REF_EMPTY, *node2, "node2 not found");
  RUS(REF_EMPTY, *node3, "node3 not found");

  return REF_SUCCESS;
}

REF_STATUS ref_swap_same_faceid(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_BOOL *allowed) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT ncell;
  REF_INT cell_to_swap[2];
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT id0, id1;

  *allowed = REF_FALSE;

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cell_to_swap),
      "more then two");

  if (0 == ncell) { /* away from boundary */
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }
  REIS(2, ncell, "there should be zero or two triangles for manifold");

  RSS(ref_cell_nodes(ref_cell, cell_to_swap[0], nodes), "nodes tri0");
  id0 = nodes[ref_cell_node_per(ref_cell)];
  RSS(ref_cell_nodes(ref_cell, cell_to_swap[1], nodes), "nodes tri1");
  id1 = nodes[ref_cell_node_per(ref_cell)];

  if (id0 == id1) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_swap_local_cell(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                               REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, search_node, test_node;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    printf("cell %d\n", cell);
    for (search_node = 0; search_node < ref_cell_node_per(ref_cell);
         search_node++) {
      if (node1 == ref_cell_c2n(ref_cell, search_node, cell)) {
        for (test_node = 0; test_node < ref_cell_node_per(ref_cell);
             test_node++) {
          printf("node %d\n", ref_cell_c2n(ref_cell, test_node, cell));
          if (!ref_node_owned(ref_node,
                              ref_cell_c2n(ref_cell, test_node, cell))) {
            return REF_SUCCESS;
          }
        }
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_swap_normdev(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                            REF_BOOL *allowed) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT ncell, cell_to_swap[2];
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node2, node3;
  REF_BOOL node0_support, node1_support;
  REF_DBL normdev0, normdev1, normdev2, normdev3;

  *allowed = REF_FALSE;

  RSS(ref_geom_supported(ref_geom, node0, &node0_support), "support0");
  RSS(ref_geom_supported(ref_geom, node1, &node1_support), "support1");
  if (!ref_geom_model_loaded(ref_geom) || !node0_support || !node1_support) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_node23(ref_grid, node0, node1, &node2, &node3), "other nodes");

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cell_to_swap),
      "more then two");
  REIS(2, ncell, "there should be two triangles for manifold");

  RSS(ref_cell_nodes(ref_cell, cell_to_swap[0], nodes), "nodes tri0");
  RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev0), "nd0");
  RSS(ref_cell_nodes(ref_cell, cell_to_swap[1], nodes), "nodes tri1");
  RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev1), "nd1");
  nodes[0] = node0;
  nodes[1] = node3;
  nodes[2] = node2;
  RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev2), "nd2");
  nodes[0] = node1;
  nodes[1] = node2;
  nodes[2] = node3;
  RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev3), "nd3");

  if ((MIN(normdev2, normdev3) > MIN(normdev0, normdev1)) ||
      (normdev2 > ref_grid_adapt(ref_grid, post_min_normdev) &&
       normdev3 > ref_grid_adapt(ref_grid, post_min_normdev))) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_swap_ratio(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node2, node3;
  REF_DBL ratio;

  *allowed = REF_FALSE;

  RSS(ref_swap_node23(ref_grid, node0, node1, &node2, &node3), "other nodes");
  RSS(ref_node_ratio(ref_node, node2, node3, &ratio), "ratio");

  if (ref_grid_adapt(ref_grid, post_min_ratio) < ratio &&
      ratio < ref_grid_adapt(ref_grid, post_max_ratio)) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}
