
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

#include "ref_fixture.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_math.h"
#include "ref_mpi.h"
#include "ref_part.h"

#define add_that_node(node, x, y, z)                                         \
  RSS(ref_node_add(ref_node, global[(node)], &(local[(node)])), "add node"); \
  ref_node_xyz(ref_node, 0, local[(node)]) = (x);                            \
  ref_node_xyz(ref_node, 1, local[(node)]) = (y);                            \
  ref_node_xyz(ref_node, 2, local[(node)]) = (z);                            \
  ref_node_part(ref_node, local[(node)]) = ref_part_implicit(                \
      nnodesg, ref_mpi_n(ref_mpi), ref_node_global(ref_node, local[(node)]));

REF_STATUS ref_fixture_tri_surf_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT nnodesg = 3;
  REF_INT cell;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_surf(ref_grid) = REF_TRUE;

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 1;
  global[1] = 2;
  global[3] = 20;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_edg(ref_grid), local, &cell), "add edg");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_tet_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT nnodesg = 4;
  REF_INT cell;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;

  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);
    add_that_node(3, 0.0, 0.0, 1.0);

    RSS(ref_cell_add(ref_grid_tet(ref_grid), local, &cell), "add tet");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 1;
  global[1] = 2;
  global[3] = 20;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_edg(ref_grid), local, &cell), "add edg");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_tet2_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT nnodesg = 5;
  REF_INT cell;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;

  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);
    add_that_node(3, 0.0, 0.0, 1.0);

    RSS(ref_cell_add(ref_grid_tet(ref_grid), local, &cell), "add tet");
  }

  global[0] = 4;
  global[1] = 2;
  global[2] = 1;
  global[3] = 3;

  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 1.0, 1.0, 1.0);
    add_that_node(1, 0.0, 1.0, 0.0);
    add_that_node(2, 1.0, 0.0, 0.0);
    add_that_node(3, 0.0, 0.0, 1.0);

    RSS(ref_cell_add(ref_grid_tet(ref_grid), local, &cell), "add tet");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pyr_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[5] = {0, 1, 2, 3, 4};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node, 0, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 1, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 2, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 3, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_add(ref_node, 4, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_initialize_n_global(ref_node, 5), "init glob");

  RSS(ref_cell_add(ref_grid_pyr(ref_grid), nodes, &cell), "add pyr");

  nodes[0] = 0;
  nodes[1] = 3;
  nodes[2] = 4;
  nodes[3] = 1;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add qua");

  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add tri");

  nodes[0] = 1;
  nodes[1] = 4;
  nodes[2] = 2;
  nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add tri");

  nodes[0] = 4;
  nodes[1] = 3;
  nodes[2] = 2;
  nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add tri");

  nodes[0] = 0;
  nodes[1] = 2;
  nodes[2] = 3;
  nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add tri");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 6;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;
  global[4] = 4;
  global[5] = 5;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 0.0, 0.0, 1.0);
    add_that_node(2, 1.0, 0.0, 0.0);
    add_that_node(3, 0.0, 1.0, 0.0);
    add_that_node(4, 0.0, 1.0, 1.0);
    add_that_node(5, 1.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 0;
  global[1] = 3;
  global[2] = 4;
  global[3] = 1;
  global[4] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  global[0] = 3;
  global[1] = 5;
  global[2] = 4;
  global[3] = 100;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  return REF_SUCCESS;
}

/*
     14
   / |  \
67 - 03 - 25
 */

REF_STATUS ref_fixture_pri2_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 8;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;
  global[4] = 4;
  global[5] = 5;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 0.0, 0.0, 1.0);
    add_that_node(2, 1.0, 0.0, 0.0);
    add_that_node(3, 0.0, 1.0, 0.0);
    add_that_node(4, 0.0, 1.0, 1.0);
    add_that_node(5, 1.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  global[0] = 0;
  global[1] = 6;
  global[2] = 1;
  global[3] = 3;
  global[4] = 7;
  global[5] = 4;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, -1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 0.0, 1.0);
    add_that_node(3, 0.0, 1.0, 0.0);
    add_that_node(4, -1.0, 1.0, 0.0);
    add_that_node(5, 0.0, 1.0, 1.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 2;
  global[1] = 5;
  global[2] = 3;
  global[3] = 0;
  global[4] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }
  global[0] = 0;
  global[1] = 3;
  global[2] = 7;
  global[3] = 6;
  global[4] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  global[0] = 3;
  global[1] = 5;
  global[2] = 4;
  global[3] = 100;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }
  global[0] = 3;
  global[1] = 4;
  global[2] = 7;
  global[3] = 100;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }
  global[0] = 0;
  global[1] = 6;
  global[2] = 1;
  global[3] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_tet_cap_grid(REF_GRID *ref_grid_ptr,
                                        REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 6;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;
  global[4] = 4;
  global[5] = 5;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);
    add_that_node(3, 0.0, 0.0, 1.0);
    add_that_node(4, 1.0, 0.0, 1.0);
    add_that_node(5, 0.0, 1.0, 1.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  global[0] = 3;
  global[1] = 4;
  global[2] = 5;
  global[3] = 6;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 0.0, 1.0);
    add_that_node(1, 1.0, 0.0, 1.0);
    add_that_node(2, 0.0, 1.0, 1.0);
    add_that_node(3, 0.3, 0.3, 2.0);

    RSS(ref_cell_add(ref_grid_tet(ref_grid), local, &cell), "add prism");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 0;
  global[1] = 3;
  global[2] = 4;
  global[3] = 1;
  global[4] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  global[0] = 3;
  global[1] = 5;
  global[2] = 4;
  global[3] = 100;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_stack_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_GLOB global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_GLOB nnodesg = 12;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;
  global[4] = 4;
  global[5] = 5;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);
    add_that_node(3, 0.0, 0.0, 1.0);
    add_that_node(4, 1.0, 0.0, 1.0);
    add_that_node(5, 0.0, 1.0, 1.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  global[0] = 3;
  global[1] = 4;
  global[2] = 5;
  global[3] = 6;
  global[4] = 7;
  global[5] = 8;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 1.0);
    add_that_node(1, 1.0, 0.0, 1.0);
    add_that_node(2, 0.0, 1.0, 1.0);
    add_that_node(3, 0.0, 0.0, 2.0);
    add_that_node(4, 1.0, 0.0, 2.0);
    add_that_node(5, 0.0, 1.0, 2.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  global[0] = 6;
  global[1] = 7;
  global[2] = 8;
  global[3] = 9;
  global[4] = 10;
  global[5] = 11;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[4]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[5])) {
    add_that_node(0, 0.0, 0.0, 2.0);
    add_that_node(1, 1.0, 0.0, 2.0);
    add_that_node(2, 0.0, 1.0, 2.0);
    add_that_node(3, 0.0, 0.0, 3.0);
    add_that_node(4, 1.0, 0.0, 3.0);
    add_that_node(5, 0.0, 1.0, 3.0);

    RSS(ref_cell_add(ref_grid_pri(ref_grid), local, &cell), "add prism");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "glob");

  global[0] = 1;
  global[1] = 0;
  global[2] = 3;
  global[3] = 4;
  global[4] = 20;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = (REF_INT)global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  global[0] = 4;
  global[1] = 3;
  global[2] = 6;
  global[3] = 7;
  global[4] = 20;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = (REF_INT)global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  global[0] = 7;
  global[1] = 6;
  global[2] = 9;
  global[3] = 10;
  global[4] = 20;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    RSS(ref_node_local(ref_node, global[2], &(local[2])), "loc");
    RSS(ref_node_local(ref_node, global[3], &(local[3])), "loc");
    local[4] = (REF_INT)global[4];
    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add quad");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_hex_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
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
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 1, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 2, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 3, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 0.0;

  RSS(ref_node_add(ref_node, 4, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_add(ref_node, 5, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_add(ref_node, 6, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_add(ref_node, 7, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 1.0;

  RSS(ref_node_initialize_n_global(ref_node, 8), "init glob");

  RSS(ref_cell_add(ref_grid_hex(ref_grid), nodes, &cell), "add prism");

  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 4;
  nodes[1] = 7;
  nodes[2] = 6;
  nodes[3] = 5;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 1;
  nodes[1] = 5;
  nodes[2] = 6;
  nodes[3] = 2;
  nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 0;
  nodes[1] = 3;
  nodes[2] = 7;
  nodes[3] = 4;
  nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 0;
  nodes[1] = 4;
  nodes[2] = 5;
  nodes[3] = 1;
  nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 2;
  nodes[1] = 6;
  nodes[2] = 7;
  nodes[3] = 3;
  nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_tet_brick_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], tet[4], cell;
  REF_INT quad[5], tri[4];

  REF_INT l = 4, m = 4, n = 4;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  REF_DBL dx, dy, dz;

  dx = (x1 - x0) / ((REF_DBL)(l - 1));
  dy = (y1 - y0) / ((REF_DBL)(m - 1));
  dz = (z1 - z0) / ((REF_DBL)(n - 1));

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i, j, k, l, m, n) ((i) + (j) * (l) + (k) * (l) * (m))

  for (k = 0; k < n; k++)
    for (j = 0; j < m; j++)
      for (i = 0; i < l; i++) {
        global = ijk2node(i, j, k, l, m, n);
        RSS(ref_node_add(ref_node, global, &node), "node");
        ref_node_xyz(ref_node, 0, node) = x0 + dx * (REF_DBL)i;
        ref_node_xyz(ref_node, 1, node) = y0 + dy * (REF_DBL)j;
        ref_node_xyz(ref_node, 2, node) = z0 + dz * (REF_DBL)k;
      }

#define ijk2hex(i, j, k, l, m, n, hex)                     \
  (hex)[0] = ijk2node((i)-1, (j)-1, (k)-1, (l), (m), (n)); \
  (hex)[1] = ijk2node((i), (j)-1, (k)-1, (l), (m), (n));   \
  (hex)[2] = ijk2node((i), (j), (k)-1, (l), (m), (n));     \
  (hex)[3] = ijk2node((i)-1, (j), (k)-1, (l), (m), (n));   \
  (hex)[4] = ijk2node((i)-1, (j)-1, (k), (l), (m), (n));   \
  (hex)[5] = ijk2node((i), (j)-1, (k), (l), (m), (n));     \
  (hex)[6] = ijk2node((i), (j), (k), (l), (m), (n));       \
  (hex)[7] = ijk2node((i)-1, (j), (k), (l), (m), (n));

  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++)
      for (i = 1; i < l; i++) {
        ijk2hex(i, j, k, l, m, n, hex);
        tet[0] = hex[0];
        tet[1] = hex[1];
        tet[2] = hex[3];
        tet[3] = hex[4];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
        tet[0] = hex[1];
        tet[1] = hex[3];
        tet[2] = hex[4];
        tet[3] = hex[5];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
        tet[0] = hex[3];
        tet[1] = hex[4];
        tet[2] = hex[5];
        tet[3] = hex[7];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
        tet[0] = hex[1];
        tet[1] = hex[2];
        tet[2] = hex[3];
        tet[3] = hex[5];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
        tet[0] = hex[2];
        tet[1] = hex[3];
        tet[2] = hex[5];
        tet[3] = hex[7];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
        tet[0] = hex[2];
        tet[1] = hex[7];
        tet[2] = hex[5];
        tet[3] = hex[6];
        RSS(ref_cell_add(ref_grid_tet(ref_grid), tet, &cell), "tet");
      }

  quad[4] = 1;
  tri[3] = quad[4];
  i = 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[0];
      quad[1] = hex[3];
      quad[2] = hex[7];
      quad[3] = hex[4];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[1];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  quad[4] = 2;
  tri[3] = quad[4];
  i = l - 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[2];
      quad[1] = hex[1];
      quad[2] = hex[5];
      quad[3] = hex[6];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[2];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[0];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  quad[4] = 3;
  tri[3] = quad[4];
  j = 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[1];
      quad[1] = hex[0];
      quad[2] = hex[4];
      quad[3] = hex[5];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[2];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[0];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  quad[4] = 4;
  tri[3] = quad[4];
  j = m - 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[3];
      quad[1] = hex[2];
      quad[2] = hex[6];
      quad[3] = hex[7];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[1];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  quad[4] = 5;
  tri[3] = quad[4];
  k = 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[0];
      quad[1] = hex[1];
      quad[2] = hex[2];
      quad[3] = hex[3];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[1];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  quad[4] = 6;
  tri[3] = quad[4];
  k = n - 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[5];
      quad[1] = hex[4];
      quad[2] = hex[7];
      quad[3] = hex[6];
      tri[0] = quad[0];
      tri[1] = quad[1];
      tri[2] = quad[2];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = quad[0];
      tri[1] = quad[2];
      tri[2] = quad[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  RSS(ref_node_initialize_n_global(ref_node, l * m * n), "init glob");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_twod_brick_grid(REF_GRID *ref_grid_ptr,
                                       REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], pri[6], cell;
  REF_INT qua[5], tri[4];

  REF_INT l = 4, m = 2, n = 4;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  REF_DBL dx, dy, dz;

  dx = (x1 - x0) / ((REF_DBL)(l - 1));
  dy = (y1 - y0) / ((REF_DBL)(m - 1));
  dz = (z1 - z0) / ((REF_DBL)(n - 1));

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  /*
                               inode7-----11-----inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             9    .        /  10    6
                            /     7           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4-8----------inode5     |
                         |      inode3.....|...5..inode2
                         |       .         |       /
                         |      .          |      /
                         |     .           |     /
                         2    1    /       4    3
  z                       |   .             |   /
  ^   y                   |  .              |  /
  |  /                    | .               | /
  | /                     |.                |/
  |/                    inode0------0-----inode1
  +----> x
  */

#define ijk2node(i, j, k, l, m, n) ((i) + (j) * (l) + (k) * (l) * (m))

  for (k = 0; k < n; k++)
    for (j = 0; j < m; j++)
      for (i = 0; i < l; i++) {
        global = ijk2node(i, j, k, l, m, n);
        RSS(ref_node_add(ref_node, global, &node), "node");
        ref_node_xyz(ref_node, 0, node) = x0 + dx * (REF_DBL)i;
        ref_node_xyz(ref_node, 1, node) = y0 + dy * (REF_DBL)j;
        ref_node_xyz(ref_node, 2, node) = z0 + dz * (REF_DBL)k;
      }

#define ijk2hex(i, j, k, l, m, n, hex)                     \
  (hex)[0] = ijk2node((i)-1, (j)-1, (k)-1, (l), (m), (n)); \
  (hex)[1] = ijk2node((i), (j)-1, (k)-1, (l), (m), (n));   \
  (hex)[2] = ijk2node((i), (j), (k)-1, (l), (m), (n));     \
  (hex)[3] = ijk2node((i)-1, (j), (k)-1, (l), (m), (n));   \
  (hex)[4] = ijk2node((i)-1, (j)-1, (k), (l), (m), (n));   \
  (hex)[5] = ijk2node((i), (j)-1, (k), (l), (m), (n));     \
  (hex)[6] = ijk2node((i), (j), (k), (l), (m), (n));       \
  (hex)[7] = ijk2node((i)-1, (j), (k), (l), (m), (n));

  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++)
      for (i = 1; i < l; i++) {
        ijk2hex(i, j, k, l, m, n, hex);
        pri[0] = hex[0];
        pri[1] = hex[4];
        pri[2] = hex[5];
        pri[3] = hex[3];
        pri[4] = hex[7];
        pri[5] = hex[6];
        RSS(ref_cell_add(ref_grid_pri(ref_grid), pri, &cell), "pri");
        pri[0] = hex[0];
        pri[1] = hex[5];
        pri[2] = hex[1];
        pri[3] = hex[3];
        pri[4] = hex[6];
        pri[5] = hex[2];
        RSS(ref_cell_add(ref_grid_pri(ref_grid), pri, &cell), "pri");
      }

  qua[4] = 3;
  i = 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[0];
      qua[1] = hex[3];
      qua[2] = hex[7];
      qua[3] = hex[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), qua, &cell), "qua");
    }

  qua[4] = 4;
  i = l - 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[2];
      qua[1] = hex[1];
      qua[2] = hex[5];
      qua[3] = hex[6];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), qua, &cell), "qua");
    }

  qua[4] = 1;
  tri[3] = qua[4];
  j = 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[1];
      qua[1] = hex[0];
      qua[2] = hex[4];
      qua[3] = hex[5];
      tri[0] = qua[0];
      tri[1] = qua[1];
      tri[2] = qua[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = qua[1];
      tri[1] = qua[2];
      tri[2] = qua[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  qua[4] = 2;
  tri[3] = qua[4];
  j = m - 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[3];
      qua[1] = hex[2];
      qua[2] = hex[6];
      qua[3] = hex[7];
      tri[0] = qua[0];
      tri[1] = qua[2];
      tri[2] = qua[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
      tri[0] = qua[0];
      tri[1] = qua[1];
      tri[2] = qua[2];
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "qua");
    }

  qua[4] = 5;
  k = 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[0];
      qua[1] = hex[1];
      qua[2] = hex[2];
      qua[3] = hex[3];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), qua, &cell), "qua");
    }

  qua[4] = 6;
  k = n - 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      qua[0] = hex[5];
      qua[1] = hex[4];
      qua[2] = hex[7];
      qua[3] = hex[6];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), qua, &cell), "qua");
    }

  return REF_SUCCESS;
}
