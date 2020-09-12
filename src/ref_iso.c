
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
  REF_INT edge, part;
  REF_INT node0, node1;
  REF_GLOB global;

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
        /* set xyz */
      }
    }
  }
  RSS(ref_node_shift_new_globals(ref_grid_node(iso_grid)), "shift iso glob");

  RSS(ref_iso_ghost(iso_grid, ref_edge, new_node), "new ghost");

  ref_free(new_node);
  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}
