
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

#include "ref_dist.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;

    REF_INT global, node, n;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    global = 0;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 1;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 1.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 2;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 1.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "tri");

    RSS(ref_dist_collisions(ref_grid, REF_FALSE, &n), "collisions");
    REIS(0, n, "no collisions expected");
    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;

    REF_INT global, node, n;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    global = 0;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 1;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 1.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 2;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 1.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "tri");
    global = 3;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.2;
    ref_node_xyz(ref_node, 1, node) = 0.2;
    ref_node_xyz(ref_node, 2, node) = -0.5;
    global = 4;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.2;
    ref_node_xyz(ref_node, 1, node) = 0.2;
    ref_node_xyz(ref_node, 2, node) = 0.5;
    global = 5;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.2;
    ref_node_xyz(ref_node, 1, node) = -0.8;
    ref_node_xyz(ref_node, 2, node) = -0.5;
    nodes[0] = 3;
    nodes[1] = 4;
    nodes[2] = 5;
    nodes[3] = 20;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "tri");

    RSS(ref_dist_collisions(ref_grid, REF_FALSE, &n), "collisions");
    REIS(2, n, "pair of collisions expected");
    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;

    REF_INT global, node, n;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    global = 0;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 1;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 1.0;
    ref_node_xyz(ref_node, 1, node) = 0.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    global = 2;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 0.0;
    ref_node_xyz(ref_node, 1, node) = 1.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "tri");
    global = 3;
    RSS(ref_node_add(ref_node, global, &node), "node");
    ref_node_xyz(ref_node, 0, node) = 1.0;
    ref_node_xyz(ref_node, 1, node) = 1.0;
    ref_node_xyz(ref_node, 2, node) = 0.0;
    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "tri");

    RSS(ref_dist_collisions(ref_grid, REF_FALSE, &n), "collisions");
    REIS(0, n, "no collisions expected, adjacent");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
