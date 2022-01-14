
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

#include "ref_axi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.c"
#include "ref_grid.h"
#include "ref_import.c"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

static REF_STATUS ref_quad_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                REF_DBL z0, REF_DBL z1) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[5] = {0, 1, 2, 3, 10};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node, 0, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_node_add(ref_node, 1, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 2, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 3, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  return REF_SUCCESS;
}

static REF_STATUS ref_prism_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 REF_DBL z0, REF_DBL z1, REF_DBL z2) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[6] = {0, 1, 2, 3, 4, 5};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node, 0, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_node_add(ref_node, 1, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 2, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.5;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z2;

  RSS(ref_node_add(ref_node, 3, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_node_add(ref_node, 4, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 5, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.5;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z2;

  RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &cell), "add prism");

  return REF_SUCCESS;
}

static REF_STATUS ref_hex_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                               REF_DBL z0, REF_DBL z1, REF_DBL z2, REF_DBL z3) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node, 0, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_node_add(ref_node, 1, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 2, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z1;

  RSS(ref_node_add(ref_node, 3, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z0;

  RSS(ref_node_add(ref_node, 4, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z3;

  RSS(ref_node_add(ref_node, 5, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = z2;

  RSS(ref_node_add(ref_node, 6, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z2;

  RSS(ref_node_add(ref_node, 7, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = z3;

  RSS(ref_cell_add(ref_grid_hex(ref_grid), nodes, &cell), "add prism");

  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");
  if (argc == 3) {
    REF_GRID ref_grid;

    printf("importing %s\n", argv[1]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "from msh");
    printf("complete.\n");

    RSS(ref_grid_inspect(ref_grid), "inspection");

    RSS(ref_axi_wedge(ref_grid), "axi wedge");

    RSS(ref_grid_inspect(ref_grid), "inspection");

    printf("exporting %s\n", argv[2]);
    RSS(ref_export_by_extension(ref_grid, argv[2]), "to file");

    RSS(ref_grid_free(ref_grid), "free");
    printf("done.\n");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* collapse node pairs : node 0 */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node0, node1;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_node_add(ref_node, 0, &node0), "add node");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 1.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    RSS(ref_node_add(ref_node, 1, &node1), "add node");
    ref_node_xyz(ref_node, 0, node1) = 0.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(1, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(REF_FALSE, ref_node_valid(ref_node, node0), "val");
    REIS(REF_TRUE, ref_node_valid(ref_node, node1), "val");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse node pairs : node 1 */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node0, node1;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_node_add(ref_node, 0, &node0), "add node");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    RSS(ref_node_add(ref_node, 1, &node1), "add node");
    ref_node_xyz(ref_node, 0, node1) = 0.0;
    ref_node_xyz(ref_node, 1, node1) = 1.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(1, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(REF_TRUE, ref_node_valid(ref_node, node0), "val");
    REIS(REF_FALSE, ref_node_valid(ref_node, node1), "val");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse quad 0 */
    REF_GRID ref_grid;

    RSS(ref_quad_grid(&ref_grid, ref_mpi, 0.0, 0.0), "quad fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(2, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse quad 1a */
    REF_GRID ref_grid;

    RSS(ref_quad_grid(&ref_grid, ref_mpi, 0.0, 0.1), "quad fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse quad 1b */
    REF_GRID ref_grid;

    RSS(ref_quad_grid(&ref_grid, ref_mpi, 0.1, 0.0), "quad fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse quad 2 */
    REF_GRID ref_grid;

    RSS(ref_quad_grid(&ref_grid, ref_mpi, 0.1, 0.1), "quad fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(4, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(1, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to prism */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.1, 0.1, 0.1), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(6, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(1, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to pyramid 0 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.0, 0.1, 0.1), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(5, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to pyramid 1 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.1, 0.0, 0.1), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(5, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to pyramid 2 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.1, 0.1, 0.0), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(5, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to tet 0 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.1, 0.0, 0.0), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(4, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to tet 1 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.0, 0.1, 0.0), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(4, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* collapse prism to tet 2 */
    REF_GRID ref_grid;

    RSS(ref_prism_grid(&ref_grid, ref_mpi, 0.0, 0.0, 0.1), "prism fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(4, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* don't collapse hex */
    REF_GRID ref_grid;

    RSS(ref_hex_grid(&ref_grid, ref_mpi, 1.0, 1.0, 2.0, 2.0), "hex fixture");
    RSS(ref_axi_wedge(ref_grid), "wedge");

    REIS(8, ref_node_n(ref_grid_node(ref_grid)), "node n");

    REIS(1, ref_cell_n(ref_grid_hex(ref_grid)), "hex n");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri n");
    REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr n");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet n");

    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
