
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

#include "ref_iso.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_node.h"

static REF_STATUS ref_iso_ghost(REF_GRID iso_grid, REF_EDGE ref_edge,
                                REF_INT *new_node) {
  REF_MPI ref_mpi = ref_grid_mpi(iso_grid);
  REF_GLOB *edge_global;
  REF_INT *edge_part;
  REF_DBL *edge_real;
  REF_INT edge, node, part, i;
  REF_GLOB global;

  ref_malloc_init(edge_global, ref_edge_n(ref_edge), REF_GLOB, REF_EMPTY);
  ref_malloc_init(edge_part, ref_edge_n(ref_edge), REF_INT, REF_EMPTY);
  ref_malloc_init(edge_real, REF_NODE_REAL_PER * ref_edge_n(ref_edge), REF_DBL,
                  -999.0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node = new_node[edge];
    if (REF_EMPTY != node) {
      RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
      if (ref_mpi_rank(ref_mpi) == part) {
        edge_global[edge] = ref_node_global(ref_grid_node(iso_grid), node);
        edge_part[edge] = ref_node_part(ref_grid_node(iso_grid), node);
        for (i = 0; i < REF_NODE_REAL_PER; i++)
          edge_real[i + REF_NODE_REAL_PER * edge] =
              ref_node_real(ref_grid_node(iso_grid), i, node);
      }
    }
  }

  RSS(ref_edge_ghost_glob(ref_edge, ref_mpi, edge_global), "global ghost");
  RSS(ref_edge_ghost_int(ref_edge, ref_mpi, edge_part), "part ghost");
  RSS(ref_edge_ghost_dbl(ref_edge, ref_mpi, edge_real, REF_NODE_REAL_PER),
      "xyz ghost");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    global = edge_global[edge];
    if (REF_EMPTY != global) {
      /* finds local if already inserted */
      RSS(ref_node_add(ref_grid_node(iso_grid), global, &node), "add node");
      new_node[edge] = node;
      ref_node_part(ref_grid_node(iso_grid), node) = edge_part[edge];
      for (i = 0; i < REF_NODE_REAL_PER; i++)
        ref_node_real(ref_grid_node(iso_grid), i, node) =
            edge_real[i + REF_NODE_REAL_PER * edge];
    }
  }

  ref_free(edge_real);
  ref_free(edge_part);
  ref_free(edge_global);

  return REF_SUCCESS;
}

REF_STATUS ref_iso_insert(REF_GRID *iso_grid_ptr, REF_GRID ref_grid,
                          REF_DBL *field) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_GRID iso_grid;
  REF_EDGE ref_edge;
  REF_INT *new_node;
  REF_INT edge, part, new_cell;
  REF_INT node0, node1, node2;
  REF_INT edge0, edge1, edge2;
  REF_GLOB global;
  REF_DBL t1, t0, d;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT id = 1;

  RSS(ref_node_synchronize_globals(ref_grid_node(ref_grid)), "sync glob");

  RSS(ref_grid_create(iso_grid_ptr, ref_mpi), "create");
  iso_grid = *iso_grid_ptr;
  RSS(ref_node_initialize_n_global(ref_grid_node(iso_grid), 0), "zero glob");

  RSS(ref_edge_create(&ref_edge, ref_grid), "create edge");
  ref_malloc_init(new_node, ref_edge_n(ref_edge), REF_INT, REF_EMPTY);

  each_ref_edge(ref_edge, edge) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (ref_mpi_rank(ref_mpi) == part) {
      if (MAX(field[node0], field[node1]) > 0.0 &&
          MIN(field[node0], field[node1]) < 0.0) {
        RSS(ref_node_next_global(ref_grid_node(iso_grid), &global),
            "next global");
        RSS(ref_node_add(ref_grid_node(iso_grid), global, &(new_node[edge])),
            "add node");
        d = field[node1] - field[node0];
        t1 = 0.5;
        if (ref_math_divisible(field[node0], d)) {
          t1 = -field[node0] / d;
        }
        t0 = 1.0 - t1;
        /*printf("%f %f d %f t %f\n",field[node0],field[node1],d,t);*/
        ref_node_xyz(ref_grid_node(iso_grid), 0, new_node[edge]) =
            t1 * ref_node_xyz(ref_grid_node(ref_grid), 0, node1) +
            t0 * ref_node_xyz(ref_grid_node(ref_grid), 0, node0);
        ref_node_xyz(ref_grid_node(iso_grid), 1, new_node[edge]) =
            t1 * ref_node_xyz(ref_grid_node(ref_grid), 1, node1) +
            t0 * ref_node_xyz(ref_grid_node(ref_grid), 1, node0);
        ref_node_xyz(ref_grid_node(iso_grid), 2, new_node[edge]) =
            t1 * ref_node_xyz(ref_grid_node(ref_grid), 2, node1) +
            t0 * ref_node_xyz(ref_grid_node(ref_grid), 2, node0);
      }
    }
  }
  RSS(ref_node_shift_new_globals(ref_grid_node(iso_grid)), "shift iso glob");

  RSS(ref_iso_ghost(iso_grid, ref_edge, new_node), "new ghost");

  if (ref_grid_twod(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_edge_with(ref_edge, nodes[1], nodes[2], &edge0), "e0");
      RSS(ref_edge_with(ref_edge, nodes[2], nodes[0], &edge1), "e1");
      RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge2), "e2");
      node0 = new_node[edge0];
      node1 = new_node[edge1];
      node2 = new_node[edge2];
      if (REF_EMPTY == node0 && REF_EMPTY != node1 && REF_EMPTY != node2) {
        new_nodes[0] = node1;
        new_nodes[1] = node2;
        new_nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(iso_grid), new_nodes, &new_cell), "add");
      }
      if (REF_EMPTY != node0 && REF_EMPTY == node1 && REF_EMPTY != node2) {
        new_nodes[0] = node2;
        new_nodes[1] = node0;
        new_nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(iso_grid), new_nodes, &new_cell), "add");
      }
      if (REF_EMPTY != node0 && REF_EMPTY != node1 && REF_EMPTY == node2) {
        new_nodes[0] = node0;
        new_nodes[1] = node1;
        new_nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(iso_grid), new_nodes, &new_cell), "add");
      }
    }
  } else {
    REF_INT cell_edge;
    REF_INT edge_nodes[6];
    REF_INT nedge;
    ref_cell = ref_grid_tet(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      nedge = 0;
      each_ref_cell_cell_edge(ref_cell, cell_edge) {
        RSS(ref_edge_with(ref_edge, ref_cell_e2n(ref_cell, 0, cell_edge, cell),
                          ref_cell_e2n(ref_cell, 1, cell_edge, cell), &edge),
            "cell edge");
        edge_nodes[cell_edge] = new_node[edge];
        if (REF_EMPTY != edge_nodes[cell_edge]) {
          new_nodes[nedge] = edge_nodes[cell_edge];
          nedge++;
        }
      }
      switch (nedge) {
        case 0:
          /* cell is not intersected by iso surface */
          break;
        case 3:
          RSS(ref_cell_add(ref_grid_tri(iso_grid), new_nodes, &new_cell),
              "add");
          break;
        default:
          printf("nedge %d\n", nedge);
          RSS(REF_IMPLEMENT, "implement iso surface tet intersection");
          break;
      }
    }
  }

  ref_free(new_node);
  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_iso_distance(REF_GRID ref_grid, REF_DBL *field,
                            REF_DBL *distance) {
  REF_GRID iso_grid;
  REF_SEARCH ref_search;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], cand[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius, dist;
  REF_DBL scale = 1.0 + 1.0e-8;
  REF_LIST ref_list;
  REF_INT node, item, candidate;

  RSS(ref_iso_insert(&iso_grid, ref_grid, field), "iso");

  ref_cell = ref_grid_edg(iso_grid);
  RSS(ref_search_create(&ref_search, ref_cell_n(ref_cell)), "create search");
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_bounding_sphere(ref_grid_node(iso_grid), nodes, 2, center,
                                 &radius),
        "b");
    RSS(ref_search_insert(ref_search, cell, center, scale * radius), "ins");
  }
  RSS(ref_list_create(&ref_list), "create list");
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    distance[node] = REF_DBL_MAX;
    RSS(ref_search_nearest_candidates(
            ref_search, ref_list,
            ref_node_xyz_ptr(ref_grid_node(ref_grid), node)),
        "candidates");
    each_ref_list_item(ref_list, item) {
      candidate = ref_list_value(ref_list, item);
      RSS(ref_cell_nodes(ref_cell, candidate, cand), "cell");
      RSS(ref_search_distance2(
              ref_node_xyz_ptr(ref_grid_node(iso_grid), cand[0]),
              ref_node_xyz_ptr(ref_grid_node(iso_grid), cand[1]),
              ref_node_xyz_ptr(ref_grid_node(ref_grid), node), &dist),
          "dist2");
      distance[node] = MIN(distance[node], dist);
    }

    RSS(ref_list_erase(ref_list), "reset list");
  }
  RSS(ref_list_free(ref_list), "free");

  RSS(ref_grid_free(iso_grid), "iso free");

  return REF_SUCCESS;
}
