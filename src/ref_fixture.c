
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

REF_STATUS ref_fixture_tri_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 3;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

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
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 0.0, 1.0, 0.0);
    add_that_node(2, 1.0, 0.0, 0.0);
    local[3] = global[3];
    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 0;
  global[1] = 1;
  global[2] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    local[2] = global[2];
    RSS(ref_cell_add(ref_grid_edg(ref_grid), local, &cell), "add edg");
  }

  return REF_SUCCESS;
}

/*
     1
   / |  \
 3 - 0 - 2
 */

REF_STATUS ref_fixture_tri2_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 4;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

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
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 0.0, 1.0, 0.0);
    add_that_node(2, 1.0, 0.0, 0.0);

    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 0;
  global[1] = 3;
  global[2] = 1;
  global[3] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, -1.0, 0.0, 0.0);
    add_that_node(2, 0.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  global[0] = 0;
  global[1] = 2;
  global[2] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    local[2] = global[2];
    RSS(ref_cell_add(ref_grid_edg(ref_grid), local, &cell), "add edg");
  }

  global[0] = 3;
  global[1] = 0;
  global[2] = 10;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1])) {
    RSS(ref_node_local(ref_node, global[0], &(local[0])), "loc");
    RSS(ref_node_local(ref_node, global[1], &(local[1])), "loc");
    local[2] = global[2];
    RSS(ref_cell_add(ref_grid_edg(ref_grid), local, &cell), "add edg");
  }

  return REF_SUCCESS;
}

/*
 3 - 1
 |   | \
 4 - 0 - 2
 */

REF_STATUS ref_fixture_tri_qua_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 5;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

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
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 0.0, 1.0, 0.0);
    add_that_node(2, 1.0, 0.0, 0.0);

    RSS(ref_cell_add(ref_grid_tri(ref_grid), local, &cell), "add tri");
  }

  global[0] = 0;
  global[1] = 4;
  global[2] = 3;
  global[3] = 1;
  global[4] = 101;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, -1.0, 0.0, 0.0);
    add_that_node(2, -1.0, 1.0, 0.0);
    add_that_node(3, 0.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_qua(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_twod_cubic_edge(REF_GRID *ref_grid_ptr,
                                       REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 4;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  global[0] = 0;
  global[1] = 2;
  global[2] = 3;
  global[3] = 1;
  global[4] = 30;
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
    add_that_node(2, 1.0 / 3.0, 0.0, 0.0);
    add_that_node(3, 2.0 / 3.0, 0.0, 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

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

  ref_grid_twod(ref_grid) = REF_FALSE;

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

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_tet_cap_grid(REF_GRID *ref_grid_ptr,
                                        REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 7;

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
/*
                               inode11----------inode10
                                 /.\               /|
                                / . \             / |
                               /  .  \           /  |
                              /   .   \         /   |
                             9    .    \       10    6
                            /     7     \     /     |
                           /      .      \   /      |
                          /       .       \ /       |
                       inode8-8----------inode9     |
                         |      inode7.....|...5..inode6
                         |       . \       |       /
                         |      .   \      |      /
                         |     .     \     |     /
                         2    1       \    4    3
                         |   .         \   |   /
                         |  .           \  |  /
                         | .             \ | /
                         |.               \|/
                       inode4------0-----inode5

                               inode7-----11-----inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             9    .           10    6
                            /     7           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4-8----------inode5     |
                         |      inode3.....|...5..inode2
                         |       .         |       /
                         |      .          |      /
                         |     .           |     /
                         2    1            4    3
                         |   .             |   /
                         |  .              |  /
                         | .               | /  z y
                         |.                |/   |/
                       inode0------0-----inode1 -> x

*/

REF_STATUS ref_fixture_hanging_hex_pri_grid(REF_GRID *ref_grid_ptr,
                                            REF_MPI ref_mpi) {
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

  /*
  nodes[0] = 4;
  nodes[1] = 7;
  nodes[2] = 6;
  nodes[3] = 5;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");
  */

  nodes[0] = 1;
  nodes[1] = 5;
  nodes[2] = 6;
  nodes[3] = 2;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 0;
  nodes[1] = 3;
  nodes[2] = 7;
  nodes[3] = 4;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 0;
  nodes[1] = 4;
  nodes[2] = 5;
  nodes[3] = 1;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 2;
  nodes[1] = 6;
  nodes[2] = 7;
  nodes[3] = 3;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  RSS(ref_node_add(ref_node, 8, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 2.0;

  RSS(ref_node_add(ref_node, 9, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 0.0;
  ref_node_xyz(ref_node, 2, node) = 2.0;

  RSS(ref_node_add(ref_node, 10, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 1.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 2.0;

  RSS(ref_node_add(ref_node, 11, &node), "add node");
  ref_node_xyz(ref_node, 0, node) = 0.0;
  ref_node_xyz(ref_node, 1, node) = 1.0;
  ref_node_xyz(ref_node, 2, node) = 2.0;

  nodes[0] = 4;
  nodes[1] = 5;
  nodes[2] = 7;
  nodes[3] = 8;
  nodes[4] = 9;
  nodes[5] = 11;
  RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 5;
  nodes[1] = 6;
  nodes[2] = 7;
  nodes[3] = 9;
  nodes[4] = 10;
  nodes[5] = 11;
  RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 5;
  nodes[1] = 9;
  nodes[2] = 10;
  nodes[3] = 6;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 4;
  nodes[1] = 7;
  nodes[2] = 11;
  nodes[3] = 8;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 4;
  nodes[1] = 8;
  nodes[2] = 9;
  nodes[3] = 5;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 6;
  nodes[1] = 10;
  nodes[2] = 11;
  nodes[3] = 7;
  nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 8;
  nodes[1] = 11;
  nodes[2] = 9;
  nodes[3] = 10;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add quad");

  nodes[0] = 9;
  nodes[1] = 11;
  nodes[2] = 10;
  nodes[3] = 10;
  RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add quad");

  return REF_SUCCESS;
}
REF_STATUS ref_fixture_hex_brick_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_INT l = 4, m = 4, n = 4;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  RSS(ref_fixture_hex_brick_args_grid(ref_grid_ptr, ref_mpi, x0, x1, y0, y1, z0,
                                      z1, l, m, n),
      "args");
  return REF_SUCCESS;
}
REF_STATUS ref_fixture_hex_brick_args_grid(REF_GRID *ref_grid_ptr,
                                           REF_MPI ref_mpi, REF_DBL x0,
                                           REF_DBL x1, REF_DBL y0, REF_DBL y1,
                                           REF_DBL z0, REF_DBL z1, REF_INT l,
                                           REF_INT m, REF_INT n) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], cell;
  REF_INT quad[5];

  REF_INT i, j, k;
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
        RSS(ref_cell_add(ref_grid_hex(ref_grid), hex, &cell), "hex");
      }

  quad[4] = 1;
  i = 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[0];
      quad[1] = hex[3];
      quad[2] = hex[7];
      quad[3] = hex[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  quad[4] = 2;
  i = l - 1;
  for (k = 1; k < n; k++)
    for (j = 1; j < m; j++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[2];
      quad[1] = hex[1];
      quad[2] = hex[5];
      quad[3] = hex[6];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  quad[4] = 3;
  j = 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[1];
      quad[1] = hex[0];
      quad[2] = hex[4];
      quad[3] = hex[5];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  quad[4] = 4;
  j = m - 1;
  for (k = 1; k < n; k++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[3];
      quad[1] = hex[2];
      quad[2] = hex[6];
      quad[3] = hex[7];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  quad[4] = 5;
  k = 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[0];
      quad[1] = hex[1];
      quad[2] = hex[2];
      quad[3] = hex[3];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  quad[4] = 6;
  k = n - 1;
  for (j = 1; j < m; j++)
    for (i = 1; i < l; i++) {
      ijk2hex(i, j, k, l, m, n, hex);
      quad[0] = hex[5];
      quad[1] = hex[4];
      quad[2] = hex[7];
      quad[3] = hex[6];
      RSS(ref_cell_add(ref_grid_qua(ref_grid), quad, &cell), "qua");
    }

  RSS(ref_node_initialize_n_global(ref_node, l * m * n), "init glob");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_tet_brick_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi) {
  REF_INT l = 4, m = 4, n = 4;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  RSS(ref_fixture_tet_brick_args_grid(ref_grid_ptr, ref_mpi, x0, x1, y0, y1, z0,
                                      z1, l, m, n),
      "args");
  return REF_SUCCESS;
}
REF_STATUS ref_fixture_tet_brick_args_grid(REF_GRID *ref_grid_ptr,
                                           REF_MPI ref_mpi, REF_DBL x0,
                                           REF_DBL x1, REF_DBL y0, REF_DBL y1,
                                           REF_DBL z0, REF_DBL z1, REF_INT l,
                                           REF_INT m, REF_INT n) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], tet[4], cell;
  REF_INT quad[5], tri[4];

  REF_INT i, j, k;
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

REF_STATUS ref_fixture_twod_brick_grid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                       REF_INT dim) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, cell;
  REF_INT edg[3], tri[4];

  REF_INT m, n;
  REF_INT i, j;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z = 0.0;

  REF_DBL dx, dy;

  m = dim;
  n = dim;

  dx = (x1 - x0) / ((REF_DBL)(n - 1));
  dy = (y1 - y0) / ((REF_DBL)(m - 1));

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  /*
  y   3 --- 2
  ^   |     |
  |   |     |
  |   0 --- 1
  |
  +----> x
  */

#define ij2node(i, j, m, n) ((i) + (j) * (m))

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++) {
      global = ij2node(i, j, m, n);
      RSS(ref_node_add(ref_node, global, &node), "node");
      ref_node_xyz(ref_node, 0, node) = x0 + dx * (REF_DBL)i;
      ref_node_xyz(ref_node, 1, node) = y0 + dy * (REF_DBL)j;
      ref_node_xyz(ref_node, 2, node) = z;
    }
  }

  edg[2] = 1;
  j = 0;
  for (i = 0; i < m - 1; i++) {
    edg[0] = ij2node(i, j, m, n);
    edg[1] = ij2node(i + 1, j, m, n);
    RSS(ref_cell_add(ref_grid_edg(ref_grid), edg, &cell), "qua");
  }

  edg[2] = 2;
  i = m - 1;
  for (j = 0; j < n - 1; j++) {
    edg[0] = ij2node(i, j, m, n);
    edg[1] = ij2node(i, j + 1, m, n);
    RSS(ref_cell_add(ref_grid_edg(ref_grid), edg, &cell), "edg");
  }

  edg[2] = 3;
  j = n - 1;
  for (i = m - 2; i >= 0; i--) {
    edg[0] = ij2node(i + 1, j, m, n);
    edg[1] = ij2node(i, j, m, n);
    RSS(ref_cell_add(ref_grid_edg(ref_grid), edg, &cell), "edg");
  }

  edg[2] = 4;
  i = 0;
  for (j = n - 2; j >= 0; j--) {
    edg[0] = ij2node(i, j + 1, m, n);
    edg[1] = ij2node(i, j, m, n);
    RSS(ref_cell_add(ref_grid_edg(ref_grid), edg, &cell), "edg");
  }

  tri[3] = 1;
  j = 1;
  for (j = 0; j < n - 1; j++) {
    for (i = 0; i < m - 1; i++) {
      tri[0] = ij2node(i, j, m, n);
      tri[1] = ij2node(i + 1, j, m, n);
      tri[2] = ij2node(i + 1, j + 1, m, n);
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "tri");
      tri[0] = ij2node(i, j, m, n);
      tri[1] = ij2node(i + 1, j + 1, m, n);
      tri[2] = ij2node(i, j + 1, m, n);
      RSS(ref_cell_add(ref_grid_tri(ref_grid), tri, &cell), "tri");
    }
  }

  RSS(ref_node_initialize_n_global(ref_node, m * n), "init glob");

  return REF_SUCCESS;
}

/*
      2
   9 8 7 6
4 10     5 3
  11     4
   0 1 2 3
      1

      21
 7 20    22 8
  19      23
 18        12
  17      13
 6 16    14 5
      15
 */

REF_STATUS ref_fixture_twod_square_circle(REF_GRID *ref_grid_ptr,
                                          REF_MPI ref_mpi) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 24;
  REF_DBL x0 = 0.5, y0 = 0.5, r = 0.25, t;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  global[0] = 0;
  global[1] = 1;
  global[2] = 2;
  global[3] = 3;
  local[4] = 1;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 0.0, 0.0);
    add_that_node(1, 1.0 / 3.0, 0.0, 0.0);
    add_that_node(2, 2.0 / 3.0, 0.0, 0.0);
    add_that_node(3, 1.0, 0.0, 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 3;
  global[1] = 4;
  global[2] = 5;
  global[3] = 6;
  local[4] = 2;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 1.0, 0.0, 0.0);
    add_that_node(1, 1.0, 1.0 / 3.0, 0.0);
    add_that_node(2, 1.0, 2.0 / 3.0, 0.0);
    add_that_node(3, 1.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 6;
  global[1] = 7;
  global[2] = 8;
  global[3] = 9;
  local[4] = 3;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 1.0, 1.0, 0.0);
    add_that_node(1, 2.0 / 3.0, 1.0, 0.0);
    add_that_node(2, 1.0 / 3.0, 1.0, 0.0);
    add_that_node(3, 0.0, 1.0, 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 9;
  global[1] = 10;
  global[2] = 11;
  global[3] = 0;
  local[4] = 4;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    add_that_node(0, 0.0, 1.0, 0.0);
    add_that_node(1, 0.0, 2.0 / 3.0, 0.0);
    add_that_node(2, 0.0, 1.0 / 3.0, 0.0);
    add_that_node(3, 0.0, 0.0, 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 12;
  global[1] = 13;
  global[2] = 14;
  global[3] = 15;
  local[4] = 5;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    t = ref_math_pi * (12.0 / 6.0);
    add_that_node(0, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (11.0 / 6.0);
    add_that_node(1, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (10.0 / 6.0);
    add_that_node(2, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (9.0 / 6.0);
    add_that_node(3, x0 + r * cos(t), y0 + r * sin(t), 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 15;
  global[1] = 16;
  global[2] = 17;
  global[3] = 18;
  local[4] = 6;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    t = ref_math_pi * (9.0 / 6.0);
    add_that_node(0, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (8.0 / 6.0);
    add_that_node(1, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (7.0 / 6.0);
    add_that_node(2, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (6.0 / 6.0);
    add_that_node(3, x0 + r * cos(t), y0 + r * sin(t), 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 18;
  global[1] = 19;
  global[2] = 20;
  global[3] = 21;
  local[4] = 7;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    t = ref_math_pi * (6.0 / 6.0);
    add_that_node(0, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (5.0 / 6.0);
    add_that_node(1, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (4.0 / 6.0);
    add_that_node(2, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (3.0 / 6.0);
    add_that_node(3, x0 + r * cos(t), y0 + r * sin(t), 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  global[0] = 21;
  global[1] = 22;
  global[2] = 23;
  global[3] = 12;
  local[4] = 8;
  if (ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[0]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[1]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[2]) ||
      ref_mpi_rank(ref_mpi) ==
          ref_part_implicit(nnodesg, ref_mpi_n(ref_mpi), global[3])) {
    t = ref_math_pi * (3.0 / 6.0);
    add_that_node(0, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (2.0 / 6.0);
    add_that_node(1, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (1.0 / 6.0);
    add_that_node(2, x0 + r * cos(t), y0 + r * sin(t), 0.0);
    t = ref_math_pi * (0.0 / 6.0);
    add_that_node(3, x0 + r * cos(t), y0 + r * sin(t), 0.0);

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), local, &cell), "add tri");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnodesg), "init glob");

  return REF_SUCCESS;
}
