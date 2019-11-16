
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

#include "ref_swap.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  { /* leave single faces alone */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    REIS(REF_INVALID, ref_swap_remove_two_face_cell(ref_grid, 0), "cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* leave different 2 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0, 3, 1, 50}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    REIS(REF_INVALID, ref_swap_remove_two_face_cell(ref_grid, 0), "cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* remove same 2 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0, 3, 1, 10}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    RSS(ref_swap_remove_two_face_cell(ref_grid, 0), "removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* leave different 3 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    nodes[0] = 0;
    nodes[1] = 3;
    nodes[2] = 1;
    nodes[3] = 50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    REIS(REF_INVALID, ref_swap_remove_three_face_cell(ref_grid, 0), "cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* remove same 3 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    nodes[0] = 0;
    nodes[1] = 3;
    nodes[2] = 1;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "other tri");

    RSS(ref_swap_remove_three_face_cell(ref_grid, 0), "removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)), "nodes");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* swap no face allowed */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "set up");

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_same_faceid(ref_grid, node0, node1, &allowed), "no face");

    REIS(REF_TRUE, allowed, "yes");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* swap same faceid allowed */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "set up");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri0");
    nodes[0] = 1;
    nodes[1] = 0;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri1");

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_same_faceid(ref_grid, node0, node1, &allowed), "no face");

    REIS(REF_TRUE, allowed, "yes");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* swap different faceid not allowed */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "set up");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri0");
    nodes[0] = 1;
    nodes[1] = 0;
    nodes[2] = 3;
    nodes[3] = 20;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri1");

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_same_faceid(ref_grid, node0, node1, &allowed), "no face");

    REIS(REF_FALSE, allowed, "no");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* local: swap allowed? */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "set up");
    RSS(ref_node_add(ref_grid_node(ref_grid), 0, &node0), "add node 0");
    RSS(ref_node_add(ref_grid_node(ref_grid), 1, &node0), "add node 1");
    RSS(ref_node_add(ref_grid_node(ref_grid), 2, &node0), "add node 2");
    RSS(ref_node_add(ref_grid_node(ref_grid), 3, &node0), "add node 3");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri0");

    nodes[0] = 1;
    nodes[1] = 0;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri1");

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_local_cell(ref_grid, node0, node1, &allowed), "no face");
    REIS(REF_TRUE, allowed, "yes");

    ref_node_part(ref_grid_node(ref_grid), 2) =
        ref_node_part(ref_grid_node(ref_grid), 2) + 1;

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_local_cell(ref_grid, node0, node1, &allowed), "no face");
    REIS(REF_FALSE, allowed, "no");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* swap nodes 23 */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT node0, node1;
    REF_INT node2, node3;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "set up");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri0");
    nodes[0] = 1;
    nodes[1] = 0;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "tri1");

    node0 = 0;
    node1 = 1;
    RSS(ref_swap_node23(ref_grid, node0, node1, &node2, &node3), "no others");
    REIS(2, node2, "yes");
    REIS(3, node3, "yes");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
