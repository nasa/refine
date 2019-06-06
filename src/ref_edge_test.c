
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
#include <string.h>

#include "ref_edge.h"

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_node.h"
#include "ref_sort.h"

#include "ref_fixture.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

#include "ref_import.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  if (2 == argc) {
    REF_GRID ref_grid;
    REF_EDGE ref_edge;
    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
    RSS(ref_edge_create(&ref_edge, ref_grid), "create");
    ref_mpi_stopwatch_stop(ref_mpi, "create");
    RSS(ref_edge_tec_fill(ref_edge, "ref_edge_test_efill_raw.tec"),
        "plot edge fill");
    RSS(ref_edge_free(ref_edge), "free");
    RSS(ref_cell_tec_fill(ref_grid_tet(ref_grid),
                          "ref_edge_test_cfill_raw.tec"),
        "plot cell fill");
    RSS(ref_adj_tec_fill(ref_cell_adj(ref_grid_tet(ref_grid)),
                         "ref_edge_test_afill_raw.tec"),
        "plot adj fill");
    ref_mpi_stopwatch_stop(ref_mpi, "plot fill");
    ref_grid_pack(ref_grid);
    ref_mpi_stopwatch_stop(ref_mpi, "pack");
    RSS(ref_edge_create(&ref_edge, ref_grid), "create");
    ref_mpi_stopwatch_stop(ref_mpi, "create");
    RSS(ref_edge_tec_fill(ref_edge, "ref_edge_test_efill_rcm.tec"),
        "plot fill");
    RSS(ref_cell_tec_fill(ref_grid_tet(ref_grid),
                          "ref_edge_test_cfill_rcm.tec"),
        "plot cell fill");
    RSS(ref_adj_tec_fill(ref_cell_adj(ref_grid_tet(ref_grid)),
                         "ref_edge_test_afill_rcm.tec"),
        "plot adj fill");
    ref_mpi_stopwatch_stop(ref_mpi, "plot fill");
    RSS(ref_edge_free(ref_edge), "free");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    return 0;
  }

  { /* make edges shared by two elements */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &cell), "add pri");

    nodes[0] = 3;
    nodes[1] = 4;
    nodes[2] = 5;
    nodes[3] = 6;
    RSS(ref_cell_add(ref_grid_tet(ref_grid), nodes, &cell), "add tet");

    RSS(ref_edge_create(&ref_edge, ref_grid), "create");

    REIS(12, ref_edge_n(ref_edge), "check total edges");

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* find edge with nodes */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT edge;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_edge_create(&ref_edge, ref_grid), "create");

    RSS(ref_edge_with(ref_edge, 0, 1, &edge), "find");
    REIS(0, edge, "right one");
    REIS(REF_NOT_FOUND, ref_edge_with(ref_edge, 0, 4, &edge), "not found");
    REIS(REF_EMPTY, edge, "not found");

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* ref_edge_ghost */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT edge, part;
    REF_INT data[9] = {REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY,
                       REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY};

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_edge_create(&ref_edge, ref_grid), "create");

    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
      if (ref_mpi_rank(ref_mpi) == part) data[edge] = edge;
    }

    RSS(ref_edge_ghost_int(ref_edge, ref_mpi, data), "ghost");

    for (edge = 0; edge < ref_edge_n(ref_edge); edge++)
      REIS(edge, data[edge], "ghost");

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* ref_edge_ghost_min_int */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT edge, part, magic, nglobal;
    REF_INT data[9] = {REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY,
                       REF_EMPTY, REF_EMPTY, REF_EMPTY, REF_EMPTY};

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_edge_create(&ref_edge, ref_grid), "create");

    nglobal = (REF_INT)ref_node_n_global(ref_edge_node(ref_edge));
    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
      magic = MIN((REF_INT)ref_node_global(ref_edge_node(ref_edge),
                                           ref_edge_e2n(ref_edge, 0, edge)),
                  (REF_INT)ref_node_global(ref_edge_node(ref_edge),
                                           ref_edge_e2n(ref_edge, 1, edge)));
      if (ref_mpi_rank(ref_mpi) == part) {
        data[edge] = nglobal;
      } else {
        data[edge] = magic;
      }
    }

    RSS(ref_edge_ghost_min_int(ref_edge, ref_mpi, data), "ghost");

    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      magic = MIN((REF_INT)ref_node_global(ref_edge_node(ref_edge),
                                           ref_edge_e2n(ref_edge, 0, edge)),
                  (REF_INT)ref_node_global(ref_edge_node(ref_edge),
                                           ref_edge_e2n(ref_edge, 1, edge)));
      RAS(magic == data[edge] || nglobal == data[edge], "not expected");
    }

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* uniq */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_node_add(ref_node, 0, &node), "first add");
    RSS(ref_node_add(ref_node, 1, &node), "first add");

    RSS(ref_edge_create(&ref_edge, ref_grid), "create");
    RSS(ref_edge_uniq(ref_edge, 0, 1), "uniq");
    REIS(1, ref_edge_n(ref_edge), "nedge");

    RSS(ref_edge_uniq(ref_edge, 0, 1), "uniq");
    REIS(1, ref_edge_n(ref_edge), "nedge");

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* rcm */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT *o2n, *n2o;
    REF_INT node;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_node_add(ref_node, 0, &node), "first add");
    RSS(ref_node_add(ref_node, 1, &node), "first add");
    RSS(ref_node_add(ref_node, 2, &node), "first add");
    RSS(ref_node_add(ref_node, 3, &node), "first add");
    RSS(ref_node_add(ref_node, 4, &node), "first add");

    /*  3 0 4 1 2
     *  1 - 2
     *   \ /
     *    0 - 4
     *    |
     *    3
     */

    RSS(ref_edge_create(&ref_edge, ref_grid), "create");
    RSS(ref_edge_uniq(ref_edge, 1, 2), "uniq");
    RSS(ref_edge_uniq(ref_edge, 1, 0), "uniq");
    RSS(ref_edge_uniq(ref_edge, 2, 0), "uniq");
    RSS(ref_edge_uniq(ref_edge, 3, 0), "uniq");
    RSS(ref_edge_uniq(ref_edge, 4, 0), "uniq");

    RSS(ref_edge_rcm(ref_edge, &o2n, &n2o), "create");

    REIS(2, n2o[0], "n2o[0]");
    REIS(1, n2o[1], "n2o[1]");
    REIS(4, n2o[2], "n2o[2]");
    REIS(0, n2o[3], "n2o[3]");
    REIS(3, n2o[4], "n2o[4]");

    REIS(3, o2n[0], "o2n[0]");
    REIS(1, o2n[1], "o2n[1]");
    REIS(0, o2n[2], "o2n[2]");
    REIS(4, o2n[3], "o2n[3]");
    REIS(2, o2n[4], "o2n[4]");

    ref_free(n2o);
    ref_free(o2n);

    RSS(ref_edge_free(ref_edge), "edge");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
