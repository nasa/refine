
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

#include "ref_grid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_shard.h"
#include "ref_sort.h"
#include "ref_validation.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  { /* init */
    REF_GRID ref_grid;
    REIS(REF_NULL, ref_grid_free(NULL), "dont free NULL");

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* each 3d element */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;
    REF_INT group;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    node_per = 3;
    each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
      node_per += 1;
      if (7 == node_per) node_per = 8;
      RES(node_per, ref_cell_node_per(ref_cell), "cells in order");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* cell with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    node_per = 4;
    RSS(ref_grid_cell_with(ref_grid, node_per, &ref_cell), "with");
    REIS(node_per, ref_cell_node_per(ref_cell), "match");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* face with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    node_per = 4;
    RSS(ref_grid_face_with(ref_grid, node_per, &ref_cell), "with");
    REIS(node_per, ref_cell_node_per(ref_cell), "match");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* unique nodes of one tri */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT nnode, ncell, *g2l, *l2g;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_cell = ref_grid_tri(ref_grid);

    nodes[0] = 5;
    nodes[1] = 8;
    nodes[2] = 6;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_grid_tri_qua_id_nodes(ref_grid, 10, &nnode, &ncell, &g2l, &l2g),
        "no list");
    REIS(1, ncell, "mis count");
    REIS(3, nnode, "mis count");
    REIS(5, l2g[0], "not in list");
    REIS(6, l2g[1], "not in list");
    REIS(8, l2g[2], "not in list");

    ref_free(l2g);
    ref_free(g2l);

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* unique nodes of two tri */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT nnode, ncell, *g2l, *l2g;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_cell = ref_grid_tri(ref_grid);

    nodes[0] = 5;
    nodes[1] = 8;
    nodes[2] = 6;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    nodes[0] = 6;
    nodes[1] = 8;
    nodes[2] = 9;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_grid_tri_qua_id_nodes(ref_grid, 10, &nnode, &ncell, &g2l, &l2g),
        "no list");
    REIS(2, ncell, "mis count");
    REIS(4, nnode, "mis count");
    REIS(5, l2g[0], "not in list");
    REIS(6, l2g[1], "not in list");
    REIS(8, l2g[2], "not in list");
    REIS(9, l2g[3], "not in list");

    ref_free(l2g);
    ref_free(g2l);

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* unique nodes of two edge */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT nnode, ncell, *g2l, *l2g;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    ref_cell = ref_grid_edg(ref_grid);

    nodes[0] = 5;
    nodes[1] = 8;
    nodes[2] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    nodes[0] = 3;
    nodes[1] = 8;
    nodes[2] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, 10, &nnode, &ncell, &g2l,
                               &l2g),
        "no list");
    REIS(2, ncell, "cell mis count");
    REIS(3, nnode, "node mis count");
    REIS(3, l2g[0], "not in list");
    REIS(5, l2g[1], "not in list");
    REIS(8, l2g[2], "not in list");
    REIS(0, g2l[3], "not in list");
    REIS(1, g2l[5], "not in list");
    REIS(2, g2l[8], "not in list");

    ref_free(l2g);
    ref_free(g2l);

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* face id flag range */
    REF_GRID ref_grid;
    REF_INT min_faceid, max_faceid;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_grid_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");
    REIS(10, min_faceid, "min");
    REIS(10, max_faceid, "max");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* orient inward tri */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_grid_tet(ref_grid), nodes, &cell), "add tet");

    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = 0;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &cell), "add outward tri");

    RSS(ref_grid_inward_boundary_orientation(ref_grid), "get out");

    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell, nodes), "get tri");
    REIS(0, nodes[0], "n0");
    REIS(1, nodes[1], "n1");
    REIS(2, nodes[2], "n2");
    REIS(10, nodes[3], "id");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* orient inward qua */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &cell), "add tri");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 4;
    nodes[3] = 3;
    nodes[4] = 20;
    RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &cell), "add outward quad");

    RSS(ref_grid_inward_boundary_orientation(ref_grid), "get out");

    RSS(ref_cell_nodes(ref_grid_qua(ref_grid), cell, nodes), "get qua");
    REIS(3, nodes[0], "n0");
    REIS(4, nodes[1], "n1");
    REIS(1, nodes[2], "n2");
    REIS(0, nodes[3], "n2");
    REIS(20, nodes[4], "id");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* around tet pri */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT nnode, list[7];

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

    RSS(ref_grid_node_list_around(ref_grid, 3, 7, &nnode, list), "no list");
    REIS(6, nnode, "mis count");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* single tet enclosing */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "fix");

    xyz[0] = 0.2;
    xyz[1] = 0.1;
    xyz[2] = 0.3;
    tet = 0;
    RSS(ref_grid_enclosing_tet(ref_grid, xyz, &tet, bary), "enclose");

    REIS(0, tet, "tet");
    RWDS(0.4, bary[0], -1, "b0");
    RWDS(0.2, bary[1], -1, "b1");
    RWDS(0.1, bary[2], -1, "b2");
    RWDS(0.3, bary[3], -1, "b3");

    xyz[0] = 0.2;
    xyz[1] = 0.1;
    xyz[2] = 0.3;
    tet = REF_EMPTY;
    RSS(ref_grid_enclosing_tet(ref_grid, xyz, &tet, bary), "enclose");

    REIS(0, tet, "tet");
    RWDS(0.4, bary[0], -1, "b0");
    RWDS(0.2, bary[1], -1, "b1");
    RWDS(0.1, bary[2], -1, "b2");
    RWDS(0.3, bary[3], -1, "b3");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* walk to find enclosing tet */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "fix");

    xyz[0] = 0.5;
    xyz[1] = 0.5;
    xyz[2] = 0.5;
    tet = REF_EMPTY;
    RSS(ref_grid_enclosing_tet(ref_grid, xyz, &tet, bary), "enclose");

    REIS(79, tet, "tet");
    RWDS(0.0, bary[0], -1, "b0");
    RWDS(0.5, bary[1], -1, "b1");
    RWDS(0.0, bary[2], -1, "b2");
    RWDS(0.5, bary[3], -1, "b3");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  { /* walk to find enclosing tet falls outside*/
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "fix");

    xyz[0] = 0.5;
    xyz[1] = 0.5;
    xyz[2] = -0.01;
    tet = REF_EMPTY;
    RSS(ref_grid_enclosing_tet(ref_grid, xyz, &tet, bary), "enclose");

    REIS(25, tet, "tet");
    RWDS(0.53, bary[0], -1, "b0");
    RWDS(0.50, bary[1], -1, "b1");
    RWDS(0.00, bary[2], -1, "b2");
    RWDS(-0.03, bary[3], -1, "b3");

    RSS(ref_grid_free(ref_grid), "cleanup");
  }

  if (!ref_mpi_para(ref_mpi)) {
    REF_GRID twod_grid, ref_grid;
    REF_BOOL twod_brick_wound_correctly = REF_FALSE;
    RSS(ref_fixture_twod_brick_grid(&twod_grid, ref_mpi, 4), "2d");
    RSS(ref_grid_extrude_twod(&ref_grid, twod_grid, 2), "extrude");

    RSS(ref_validation_all(ref_grid), "pri val");

    if (twod_brick_wound_correctly) {
      REF_GRID tet_grid;
      RSS(ref_grid_deep_copy(&tet_grid, ref_grid), "deep copy");
      RSS(ref_shard_prism_into_tet(tet_grid, 0, REF_EMPTY), "shard to tet");
      RSS(ref_validation_all(ref_grid), "pri val");
      RSS(ref_grid_free(tet_grid), "cleanup");
    }

    RSS(ref_grid_free(ref_grid), "cleanup");
    RSS(ref_grid_free(twod_grid), "cleanup");
  }

  if (!ref_mpi_para(ref_mpi)) {
    REF_GRID twod_grid, ref_grid;
    REF_BOOL twod_brick_wound_correctly = REF_FALSE;
    RSS(ref_fixture_twod_brick_grid(&twod_grid, ref_mpi, 4), "2d");
    RSS(ref_grid_extrude_twod(&ref_grid, twod_grid, 5), "extrude");

    RSS(ref_validation_all(ref_grid), "pri val");

    if (twod_brick_wound_correctly) {
      REF_GRID tet_grid;
      RSS(ref_grid_deep_copy(&tet_grid, ref_grid), "deep copy");
      RSS(ref_shard_prism_into_tet(tet_grid, 0, REF_EMPTY), "shard to tet");
      RSS(ref_validation_all(ref_grid), "pri val");
      RSS(ref_grid_free(tet_grid), "cleanup");
    }

    RSS(ref_grid_free(ref_grid), "cleanup");
    RSS(ref_grid_free(twod_grid), "cleanup");
  }

  if (!ref_mpi_para(ref_mpi)) { /* count total local cell */
    REF_GRID ref_grid;
    REF_INT ncell;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_grid_ncell(ref_grid, &ncell), "total cells");
    REIS(3, ncell, "ncell");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* map contiguous cell index tet */
    REF_GRID ref_grid;
    REF_INT contiguous_cell, cell_group, cell;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");

    contiguous_cell = REF_EMPTY;
    REIS(REF_NOT_FOUND,
         ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                        &cell),
         "map to ref_cell");
    REIS(REF_EMPTY, cell_group, "cell_group");
    REIS(REF_EMPTY, cell, "cell");

    contiguous_cell = 0;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_EDG, cell_group, "cell_group");
    REIS(0, cell, "cell");

    contiguous_cell = 1;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_TRI, cell_group, "cell_group");
    REIS(0, cell, "cell");

    contiguous_cell = 2;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_TET, cell_group, "cell_group");
    REIS(0, cell, "cell");

    contiguous_cell = 3;
    REIS(REF_NOT_FOUND,
         ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                        &cell),
         "map to ref_cell");
    REIS(REF_EMPTY, cell_group, "cell_group");
    REIS(REF_EMPTY, cell, "cell");

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* map contiguous cell index tri2 */
    REF_GRID ref_grid;
    REF_INT contiguous_cell, cell_group, cell;
    RSS(ref_fixture_tri2_grid(&ref_grid, ref_mpi), "set up tet");

    contiguous_cell = 0;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_EDG, cell_group, "cell_group");
    REIS(0, cell, "cell");

    contiguous_cell = 1;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_EDG, cell_group, "cell_group");
    REIS(1, cell, "cell");

    contiguous_cell = 2;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_TRI, cell_group, "cell_group");
    REIS(0, cell, "cell");

    contiguous_cell = 3;
    RSS(ref_grid_contiguous_group_cell(ref_grid, contiguous_cell, &cell_group,
                                       &cell),
        "map to ref_cell");
    REIS(REF_CELL_TRI, cell_group, "cell_group");
    REIS(1, cell, "cell");

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* contiguous cell global tet */
    REF_GRID ref_grid;
    REF_INT ncell;
    REF_LONG *global;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_grid_ncell(ref_grid, &ncell), "total cells");
    REIS(3, ncell, "ref_fixture_tet_grid changed");
    ref_malloc(global, ncell, REF_LONG);
    RSS(ref_grid_contiguous_cell_global(ref_grid, global), "global");
    REIS(0, global[0], "global[0]");
    REIS(1, global[1], "global[1]");
    REIS(2, global[2], "global[2]");
    ref_free(global);
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* contiguous cell global tri2 */
    REF_GRID ref_grid;
    REF_INT ncell;
    REF_LONG *global;
    RSS(ref_fixture_tri2_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_grid_ncell(ref_grid, &ncell), "total cells");
    REIS(4, ncell, "ref_fixture_tri2_grid changed");
    ref_malloc(global, ncell, REF_LONG);
    RSS(ref_grid_contiguous_cell_global(ref_grid, global), "global");
    REIS(0, global[0], "global[0]");
    REIS(1, global[1], "global[1]");
    REIS(2, global[2], "global[2]");
    REIS(3, global[3], "global[3]");
    ref_free(global);
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
