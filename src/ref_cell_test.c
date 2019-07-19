
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

#include "ref_adj.h"
#include "ref_cell.h"

#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

#include "ref_fixture.h"
#include "ref_grid.h"

#include "ref_malloc.h"
#include "ref_math.h"

static REF_STATUS ref_edg(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 2, REF_TRUE);
}
static REF_STATUS ref_tri(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 3, REF_TRUE);
}
static REF_STATUS ref_qua(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 4, REF_TRUE);
}

static REF_STATUS ref_tet(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 4, REF_FALSE);
}
static REF_STATUS ref_pyr(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 5, REF_FALSE);
}
static REF_STATUS ref_pri(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 6, REF_FALSE);
}
static REF_STATUS ref_hex(REF_CELL *ref_cell_ptr) {
  return ref_cell_create(ref_cell_ptr, 8, REF_FALSE);
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create mpi");

  REIS(REF_NULL, ref_cell_free(NULL), "dont free NULL");

  { /* deep copy empty */
    REF_CELL ref_cell;
    REF_CELL original;

    RSS(ref_tet(&original), "create");
    RSS(ref_cell_deep_copy(&ref_cell, original), "deep copy");

    RSS(ref_cell_free(original), "cleanup");
    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* add */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tet(&ref_cell), "create");
    RES(0, ref_cell_n(ref_cell), "init zero cells");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    RES(0, cell, "first cell is zero");
    RES(1, ref_cell_n(ref_cell), "first cell incements n");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    RES(1, cell, "second cell is one");
    RES(2, ref_cell_n(ref_cell), "second cell incements n");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* add many global */
    REF_CELL ref_cell;
    REF_NODE ref_node;
    REF_GLOB nodes[4];
    REF_INT parts[4];
    REF_INT retrieved[4];

    RSS(ref_node_create(&ref_node, ref_mpi), "create node");

    RSS(ref_tet(&ref_cell), "create");
    RES(0, ref_cell_n(ref_cell), "init zero cells");

    nodes[0] = 0;
    nodes[1] = 10;
    nodes[2] = 20;
    nodes[3] = 30;
    parts[0] = 0;
    parts[1] = 1;
    parts[2] = 2;
    parts[3] = 3;
    RSS(ref_cell_add_many_global(ref_cell, ref_node, 1, nodes, parts,
                                 REF_EMPTY),
        "add many");

    RSS(ref_cell_nodes(ref_cell, 0, retrieved), "cell should exist");
    RES(0, retrieved[0], "node 0");
    RES(1, retrieved[1], "node 1");
    RES(2, retrieved[2], "node 2");
    RES(3, retrieved[3], "node 3");

    RES(20, ref_node_global(ref_node, 2), "global mapped");
    RES(2, ref_node_part(ref_node, 2), "part rec");

    RSS(ref_cell_free(ref_cell), "cleanup");
    RSS(ref_node_free(ref_node), "cleanup");
  }

  { /* ref_cell_nodes bounds */
    REF_CELL ref_cell;
    REF_INT max;
    REF_INT cell, nodes[4];
    RSS(ref_tet(&ref_cell), "create");
    max = ref_cell_max(ref_cell);
    cell = -1;
    REIS(REF_INVALID, ref_cell_nodes(ref_cell, cell, nodes), "invalid cell");
    cell = max;
    REIS(REF_INVALID, ref_cell_nodes(ref_cell, cell, nodes), "invalid cell");
  }

  { /* remove */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT item;
    REF_INT cell;

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    nodes[0] = 0;
    nodes[1] = 4;
    nodes[2] = 5;
    nodes[3] = 6;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    REIS(REF_INVALID, ref_cell_remove(ref_cell, 3), "remove cell missing cell");
    RES(2, ref_cell_n(ref_cell), "still there");

    RSS(ref_cell_remove(ref_cell, 0), "remove cell");
    RES(1, ref_cell_n(ref_cell), "reduced count")
    RAS(!ref_cell_valid(ref_cell, 0), "cell is still here");

    each_ref_adj_node_item_with_ref((ref_cell)->ref_adj, 0, item, cell)
        RAS(!(cell == 0), "removed cell still in adj");

    each_ref_adj_node_item_with_ref((ref_cell)->ref_adj, 2, item, cell)
        RAS(!(cell == 0), "removed cell still in adj");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* compact */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;
    REF_INT *o2n, *n2o;

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 4;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    nodes[0] = 3;
    nodes[1] = 4;
    nodes[2] = 5;
    nodes[3] = 6;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_remove(ref_cell, 1), "remove cell");

    RSS(ref_cell_compact(ref_cell, &o2n, &n2o), "compact");

    REIS(0, o2n[0], "o2n");
    REIS(REF_EMPTY, o2n[1], "o2n");
    REIS(1, o2n[2], "o2n");

    REIS(0, n2o[0], "n2o");
    REIS(2, n2o[1], "n2o");

    ref_free(n2o);
    ref_free(o2n);

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* remove tri*/
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 1;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_remove(ref_cell, cell), "remove cell");

    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 1), "id");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* replace whole tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = 3;
    nodes[3] = 20;
    RSS(ref_cell_replace_whole(ref_cell, cell, nodes), "replace cell");

    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 0), "old node");
    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 3), "new node");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* replace node of tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT retrieved[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_replace_node(ref_cell, 2, 3), "replace node");

    RSS(ref_cell_nodes(ref_cell, cell, retrieved), "cell should exist");
    REIS(0, retrieved[0], "node 0");
    REIS(1, retrieved[1], "node 1");
    REIS(3, retrieved[2], "node 2");
    REIS(10, retrieved[3], "id");

    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 2), "old node");
    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 3), "new node");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* replace node of 2 tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = 3;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_replace_node(ref_cell, 2, 4), "replace node");

    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 2), "old node");
    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 4), "new node");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* force realloc twice */
    REF_CELL ref_cell;
    REF_INT max, i;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell), "create");

    max = ref_cell_max(ref_cell);
    for (i = 0; i < max + 1; i++) {
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      nodes[3] = 3;
      RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
      RES(i, cell, "expected ordered new cells first block");
    }
    RAS(ref_cell_max(ref_cell) > max, "realloc max");

    max = ref_cell_max(ref_cell);
    for (i = ref_cell_n(ref_cell); i < max + 1; i++) {
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      nodes[3] = 3;
      RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
      RES(i, cell, "expected ordered new cells second block");
    }
    RAS(ref_cell_max(ref_cell) > max, "realloc max");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* reuse without reallocation */
    REF_CELL ref_cell;
    REF_INT max, i;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell), "create");

    max = ref_cell_max(ref_cell);
    for (i = 0; i < max + 10; i++) {
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      nodes[3] = 3;
      RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
      RSS(ref_cell_remove(ref_cell, cell), "rm cell");
    }
    REIS(max, ref_cell_max(ref_cell), "do not realloc max");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* get nodes */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT retrieved[4];

    RSS(ref_tet(&ref_cell), "create new");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    REIS(REF_INVALID, ref_cell_nodes(ref_cell, 0, nodes), "missing cell");
    REIS(REF_INVALID, ref_cell_nodes(ref_cell, -1, nodes),
         "-1 cell should fail");
    REIS(REF_INVALID, ref_cell_nodes(ref_cell, 1000000000, nodes),
         "large cell");

    nodes[0] = 10;
    nodes[1] = 20;
    nodes[2] = 30;
    nodes[3] = 40;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_nodes(ref_cell, 0, retrieved), "cell should exist");
    RES(10, retrieved[0], "node 0");
    RES(20, retrieved[1], "node 1");
    RES(30, retrieved[2], "node 2");
    RES(40, retrieved[3], "node 3");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* valid? */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    RSS(ref_tet(&ref_cell), "create new");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RAS(!ref_cell_valid(ref_cell, -1), "invlid -1");
    RAS(!ref_cell_valid(ref_cell, -1), "invlid 0");
    RAS(!ref_cell_valid(ref_cell, -1), "invlid 1");
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    RAS(ref_cell_valid(ref_cell, 0), "valid 0");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* orient */
    REF_INT nnode, node0, nodes[4];
    nnode = 4;
    node0 = 0;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_orient_node0(nnode, node0, nodes), "orient");
    REIS(0, nodes[0], "nodes[0]");
    REIS(1, nodes[1], "nodes[1]");
    REIS(2, nodes[2], "nodes[2]");
    REIS(3, nodes[3], "nodes[3]");
    node0 = 1;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_orient_node0(nnode, node0, nodes), "orient");
    REIS(1, nodes[0], "nodes[0]");
    REIS(0, nodes[1], "nodes[1]");
    REIS(3, nodes[2], "nodes[2]");
    REIS(2, nodes[3], "nodes[3]");
    node0 = 2;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_orient_node0(nnode, node0, nodes), "orient");
    REIS(2, nodes[0], "nodes[0]");
    REIS(3, nodes[1], "nodes[1]");
    REIS(0, nodes[2], "nodes[2]");
    REIS(1, nodes[3], "nodes[3]");
    node0 = 3;
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_orient_node0(nnode, node0, nodes), "orient");
    REIS(3, nodes[0], "nodes[0]");
    REIS(2, nodes[1], "nodes[1]");
    REIS(1, nodes[2], "nodes[2]");
    REIS(0, nodes[3], "nodes[3]");
  }

  { /* adj */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell), "create new");

    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 0), "first node");
    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 3), "last node");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 0), "first node");
    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 3), "last node");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* adj with id */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];

    RSS(ref_tri(&ref_cell), "create new");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 0), "first node");
    RAS(!ref_adj_empty(ref_cell_adj(ref_cell), 2), "last node");
    RAS(ref_adj_empty(ref_cell_adj(ref_cell), 3), "id");

    RSS(ref_cell_free(ref_cell), "free");
  }

  { /* loop cells */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell;

    RSS(ref_tet(&ref_cell), "create new");

    ncell = 0;
    each_ref_cell_valid_cell(ref_cell, cell) ncell++;

    RES(0, ncell, "start zero cells");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    nodes[0] = -1;
    nodes[1] = -1;
    nodes[2] = -1;
    nodes[3] = -1;

    ncell = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) ncell++;

    RES(1, ncell, "found cell");
    RES(0, nodes[0], "got node 0");
    RES(1, nodes[1], "got node 1");
    RES(2, nodes[2], "got node 2");
    RES(3, nodes[3], "got node 3");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* edge_per */
    REF_CELL ref_cell;

    RSS(ref_tet(&ref_cell), "create");
    REIS(6, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_pyr(&ref_cell), "create");
    REIS(8, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_pri(&ref_cell), "create");
    REIS(9, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_hex(&ref_cell), "create");
    REIS(12, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_edg(&ref_cell), "create");
    REIS(1, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_tri(&ref_cell), "create");
    REIS(3, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_qua(&ref_cell), "create");
    REIS(4, ref_cell_edge_per(ref_cell), "edge_per");
    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* face_per */
    REF_CELL ref_cell;

    RSS(ref_tet(&ref_cell), "create");
    RES(4, ref_cell_face_per(ref_cell), "face_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_pyr(&ref_cell), "create");
    RES(5, ref_cell_face_per(ref_cell), "face_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_pri(&ref_cell), "create");
    RES(5, ref_cell_face_per(ref_cell), "face_per");
    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_hex(&ref_cell), "create");
    RES(6, ref_cell_face_per(ref_cell), "face_per");
    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* tet face */
    REF_CELL ref_cell;
    RSS(ref_tet(&ref_cell), "create");
    REIS(1, ref_cell_f2n_gen(ref_cell, 0, 0), "tri face nodes");
    REIS(3, ref_cell_f2n_gen(ref_cell, 1, 0), "tri face nodes");
    REIS(2, ref_cell_f2n_gen(ref_cell, 2, 0), "tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell, 0, 0), ref_cell_f2n_gen(ref_cell, 3, 0),
         "quad not tri");
    REIS(0, ref_cell_f2n_gen(ref_cell, 0, 1), "tri face nodes");
    REIS(2, ref_cell_f2n_gen(ref_cell, 1, 1), "tri face nodes");
    REIS(3, ref_cell_f2n_gen(ref_cell, 2, 1), "tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell, 0, 1), ref_cell_f2n_gen(ref_cell, 3, 1),
         "quad not tri");
    REIS(0, ref_cell_f2n_gen(ref_cell, 0, 2), "tri face nodes");
    REIS(3, ref_cell_f2n_gen(ref_cell, 1, 2), "tri face nodes");
    REIS(1, ref_cell_f2n_gen(ref_cell, 2, 2), "tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell, 0, 2), ref_cell_f2n_gen(ref_cell, 3, 2),
         "quad not tri");
    REIS(0, ref_cell_f2n_gen(ref_cell, 0, 3), "tri face nodes");
    REIS(1, ref_cell_f2n_gen(ref_cell, 1, 3), "tri face nodes");
    REIS(2, ref_cell_f2n_gen(ref_cell, 2, 3), "tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell, 0, 3), ref_cell_f2n_gen(ref_cell, 3, 3),
         "quad not tri");
    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* tet face nodes*/
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT cell_face;
    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 10;
    nodes[1] = 20;
    nodes[2] = 30;
    nodes[3] = 40;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    cell = 0;
    cell_face = 3;
    REIS(10, ref_cell_f2n(ref_cell, 0, cell_face, cell), "tri face nodes");
    REIS(20, ref_cell_f2n(ref_cell, 1, cell_face, cell), "tri face nodes");
    REIS(30, ref_cell_f2n(ref_cell, 2, cell_face, cell), "tri face nodes");
    REIS(ref_cell_f2n(ref_cell, 0, cell_face, cell),
         ref_cell_f2n(ref_cell, 3, cell_face, cell), "tri face nodes");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* hex edge faces */
    REF_CELL ref_cell;
    REF_INT edge, face0, face1;
    RSS(ref_hex(&ref_cell), "create");

    edge = 0;
    RSS(ref_cell_gen_edge_face(ref_cell, edge, &face0, &face1), "edge face");
    REIS(0, face0, "face0");
    REIS(4, face1, "face1");

    edge = 1;
    RSS(ref_cell_gen_edge_face(ref_cell, edge, &face0, &face1), "edge face");
    REIS(3, face0, "face0");
    REIS(4, face1, "face1");

    edge = 2;
    RSS(ref_cell_gen_edge_face(ref_cell, edge, &face0, &face1), "edge face");
    REIS(0, face0, "face0");
    REIS(3, face1, "face1");

    edge = 6;
    RSS(ref_cell_gen_edge_face(ref_cell, edge, &face0, &face1), "edge face");
    REIS(1, face0, "face0");
    REIS(2, face1, "face1");

    edge = 11;
    RSS(ref_cell_gen_edge_face(ref_cell, edge, &face0, &face1), "edge face");
    REIS(2, face0, "face0");
    REIS(5, face1, "face1");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* tet with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT found;

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_with(ref_cell, nodes, &found), "not found");
    REIS(cell, found, "not same");

    nodes[0] = 5;
    REIS(REF_NOT_FOUND, ref_cell_with(ref_cell, nodes, &found), "found");
    REIS(REF_EMPTY, found, "not same")

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* list of one tet with two nodes */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell, list[5];

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_list_with2(ref_cell, 0, 1, 5, &ncell, list), "no list");
    REIS(1, ncell, "mis count");
    REIS(0, list[0], "not in list");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* list of two tet with two nodes */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell, list[5];

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 1;
    nodes[3] = 4;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    REIS(REF_INCREASE_LIMIT,
         ref_cell_list_with2(ref_cell, 0, 1, 1, &ncell, list),
         "one should be too small");

    RSS(ref_cell_list_with2(ref_cell, 0, 1, 5, &ncell, list), "no list");
    REIS(2, ncell, "mis count");
    REIS(1, list[0], "not in list");
    REIS(0, list[1], "not in list");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* degree of tets with two nodes */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT degree;

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 1;
    nodes[3] = 4;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_degree_with2(ref_cell, 3, 4, &degree), "degree 0");
    REIS(0, degree, "zero edge");
    RSS(ref_cell_degree_with2(ref_cell, 0, 3, &degree), "degree 1");
    REIS(1, degree, "one edge");
    RSS(ref_cell_degree_with2(ref_cell, 0, 1, &degree), "degree 2");
    REIS(2, degree, "two edge");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* list of node */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT nnode, list[6];

    RSS(ref_tet(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_node_list_around(ref_cell, 0, 6, &nnode, list), "no list");
    REIS(3, nnode, "mis count");
    REIS(1, list[0], "not in list");
    REIS(2, list[1], "not in list");
    REIS(3, list[2], "not in list");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 1;
    nodes[3] = 4;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_node_list_around(ref_cell, 0, 6, &nnode, list), "no list");
    REIS(4, nnode, "mis count");
    REIS(2, list[0], "not in list");
    REIS(1, list[1], "not in list");
    REIS(4, list[2], "not in list");
    REIS(3, list[3], "not in list");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* list of faceid */
    REF_CELL ref_cell;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT nfaceid, faceids[3];

    RSS(ref_tri(&ref_cell), "create");

    RSS(ref_cell_id_list_around(ref_cell, 0, 3, &nfaceid, faceids), "no list");
    REIS(0, nfaceid, "mis count");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add tri");

    RSS(ref_cell_id_list_around(ref_cell, 0, 3, &nfaceid, faceids), "list1");
    REIS(1, nfaceid, "mis count");
    REIS(10, faceids[0], "not in list");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 20;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_id_list_around(ref_cell, 0, 3, &nfaceid, faceids), "list2");
    REIS(2, nfaceid, "mis count");
    REIS(20, faceids[0], "not in list");
    REIS(10, faceids[1], "not in list");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 30;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_id_list_around(ref_cell, 0, 3, &nfaceid, faceids), "list3");
    REIS(3, nfaceid, "mis count");
    REIS(30, faceids[0], "not in list");
    REIS(20, faceids[1], "not in list");
    REIS(10, faceids[2], "not in list");

    nodes[0] = 0;
    nodes[1] = 2;
    nodes[2] = 3;
    nodes[3] = 40;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    REIS(REF_INCREASE_LIMIT,
         ref_cell_id_list_around(ref_cell, 0, 3, &nfaceid, faceids), "list4");
    REIS(3, nfaceid, "mis count");
    REIS(40, faceids[0], "not in list");
    REIS(30, faceids[1], "not in list");
    REIS(20, faceids[2], "not in list");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* tet tris */
    REF_CELL ref_cell;
    REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT ntri, cell;

    RSS(ref_tri(&ref_cell), "create");

    tet_nodes[0] = 0;
    tet_nodes[1] = 1;
    tet_nodes[2] = 2;
    tet_nodes[3] = 3;

    RSS(ref_cell_ntri_with_tet_nodes(ref_cell, tet_nodes, &ntri), "tet tris");
    REIS(0, ntri, "should be none");

    tri_nodes[0] = 0;
    tri_nodes[1] = 1;
    tri_nodes[2] = 2;
    tri_nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, tri_nodes, &cell), "add tri");

    RSS(ref_cell_ntri_with_tet_nodes(ref_cell, tet_nodes, &ntri), "tet tris");
    REIS(1, ntri, "should be one");

    tri_nodes[0] = 0;
    tri_nodes[1] = 3;
    tri_nodes[2] = 1;
    tri_nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, tri_nodes, &cell), "add tri");

    RSS(ref_cell_ntri_with_tet_nodes(ref_cell, tet_nodes, &ntri), "tet tris");
    REIS(2, ntri, "should be two");

    tri_nodes[0] = 1;
    tri_nodes[1] = 3;
    tri_nodes[2] = 2;
    tri_nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, tri_nodes, &cell), "add tri");

    RSS(ref_cell_ntri_with_tet_nodes(ref_cell, tet_nodes, &ntri), "tet tris");
    REIS(3, ntri, "should be three");

    tri_nodes[0] = 0;
    tri_nodes[1] = 2;
    tri_nodes[2] = 3;
    tri_nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, tri_nodes, &cell), "add tri");

    RSS(ref_cell_ntri_with_tet_nodes(ref_cell, tet_nodes, &ntri), "tet tris");
    REIS(4, ntri, "should be all");

    RSS(ref_cell_free(ref_cell), "cleanup tri");
  }

  { /* tri with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT found;

    RSS(ref_tri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    RSS(ref_cell_with(ref_cell, nodes, &found), "not found");
    REIS(cell, found, "not same");

    nodes[3] = 5;
    RSS(ref_cell_with(ref_cell, nodes, &found), "not found");
    REIS(cell, found, "not same");

    nodes[0] = 5;
    REIS(REF_NOT_FOUND, ref_cell_with(ref_cell, nodes, &found), "found");
    REIS(REF_EMPTY, found, "not same")

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* pri cell has face */
    REF_CELL ref_cell;
    REF_INT cell;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT face_nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT found0, found1;

    RSS(ref_pri(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    face_nodes[0] = 0;
    face_nodes[1] = 1;
    face_nodes[2] = 5;
    face_nodes[3] = 0;
    RSS(ref_cell_with_face(ref_cell, face_nodes, &found0, &found1), "with");
    REIS(REF_EMPTY, found0, "false positive");
    REIS(REF_EMPTY, found1, "false positive");

    face_nodes[0] = 1;
    face_nodes[1] = 0;
    face_nodes[2] = 3;
    face_nodes[3] = 4;
    RSS(ref_cell_with_face(ref_cell, face_nodes, &found0, &found1), "with");
    REIS(0, found0, "false negative");
    REIS(REF_EMPTY, found1, "false positive");

    face_nodes[0] = 0;
    face_nodes[1] = 1;
    face_nodes[2] = 2;
    face_nodes[3] = 0;
    RSS(ref_cell_with_face(ref_cell, face_nodes, &found0, &found1), "with");
    REIS(0, found0, "false negative");
    REIS(REF_EMPTY, found1, "false positive");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* has side */
    REF_CELL ref_cell;
    REF_INT cell, nodes[8];
    REF_INT node0, node1;
    REF_BOOL has_side;

    RSS(ref_hex(&ref_cell), "create");

    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 4;
    nodes[4] = 4;
    nodes[5] = 5;
    nodes[6] = 6;
    nodes[7] = 7;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    node0 = 0;
    node1 = 1;
    has_side = REF_FALSE;
    RSS(ref_cell_has_side(ref_cell, node0, node1, &has_side), "side");
    REIS(REF_TRUE, has_side, "side expected");

    node0 = 0;
    node1 = 2;
    RSS(ref_cell_has_side(ref_cell, node0, node1, &has_side), "side");
    REIS(REF_FALSE, has_side, "diagonal, not side");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* triangle has face */
    REF_CELL ref_cell;
    REF_INT cell0, cell1, nodes[8];

    RSS(ref_tri(&ref_cell), "create");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;

    RSS(ref_cell_add(ref_cell, nodes, &cell0), "add cell");

    nodes[3] = 0;
    cell0 = 55;
    cell1 = 55;
    RSS(ref_cell_with_face(ref_cell, nodes, &cell0, &cell1), "has");
    REIS(0, cell0, "wrong cell");
    REIS(REF_EMPTY, cell1, "false pos");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 0;
    nodes[3] = 1;
    cell0 = 55;
    cell1 = 55;
    RSS(ref_cell_with_face(ref_cell, nodes, &cell0, &cell1), "has");
    REIS(0, cell0, "wrong cell");
    REIS(REF_EMPTY, cell1, "false pos");

    RSS(ref_cell_free(ref_cell), "cleanup");
  }

  { /* cell part */
    REF_CELL ref_cell;
    REF_NODE ref_node;
    REF_INT cell, nodes[3];
    REF_INT global, node;
    REF_INT part;

    RSS(ref_node_create(&ref_node, ref_mpi), "create node");
    global = 0;
    RSS(ref_node_add(ref_node, global, &node), "first add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &node), "first add");

    RSS(ref_edg(&ref_cell), "create");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 10;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    ref_node_part(ref_node, 0) = 0;
    ref_node_part(ref_node, 1) = 1;
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "get part");
    REIS(0, part, "wrong part");

    ref_node_part(ref_node, 0) = 3;
    ref_node_part(ref_node, 1) = 2;
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "get part");
    REIS(3, part, "wrong part");

    RSS(ref_cell_free(ref_cell), "cleanup");

    RSS(ref_edg(&ref_cell), "create");
    nodes[0] = 1;
    nodes[1] = 0;
    nodes[2] = 20;
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");

    ref_node_part(ref_node, 0) = 1;
    ref_node_part(ref_node, 1) = 0;
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "get part");
    REIS(1, part, "wrong part");

    ref_node_part(ref_node, 0) = 2;
    ref_node_part(ref_node, 1) = 3;
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "get part");
    REIS(2, part, "wrong part");

    RSS(ref_cell_free(ref_cell), "cleanup");
    RSS(ref_node_free(ref_node), "cleanup");
  }

  { /* gem local? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    if (ref_cell_valid(ref_grid_tet(ref_grid), 0)) {
      node0 = 0;
      node1 = 1;
      RSS(ref_cell_local_gem(ref_grid_tet(ref_grid), ref_grid_node(ref_grid),
                             node0, node1, &allowed),
          "split");

      REIS(!ref_mpi_para(ref_grid_mpi(ref_grid)), allowed,
           "local split allowed?");
    }

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (4 <= ref_mpi_n(ref_mpi)) { /* cell ghost */
    REF_CELL ref_cell;
    REF_NODE ref_node;
    REF_INT global, nodes[4], cell;
    REF_LONG *data;

    RSS(ref_tet(&ref_cell), "create");
    RSS(ref_node_create(&ref_node, ref_mpi), "create node");

    if (0 <= ref_mpi_rank(ref_mpi) && ref_mpi_rank(ref_mpi) <= 3) {
      global = 0;
      RSS(ref_node_add(ref_node, global, &(nodes[0])), "node");
      ref_node_part(ref_node, nodes[0]) = 0;
      global = 1;
      RSS(ref_node_add(ref_node, global, &(nodes[1])), "node");
      ref_node_part(ref_node, nodes[1]) = 1;
      global = 2;
      RSS(ref_node_add(ref_node, global, &(nodes[2])), "node");
      ref_node_part(ref_node, nodes[2]) = 2;
      global = 3;
      RSS(ref_node_add(ref_node, global, &(nodes[3])), "node");
      ref_node_part(ref_node, nodes[3]) = 3;
      RSS(ref_cell_add(ref_cell, nodes, &cell), "add cell");
    }

    ref_malloc(data, 1, REF_LONG);
    data[0] = REF_EMPTY;
    if (ref_mpi_rank(ref_mpi) == 0) {
      data[0] = 10;
    }
    RSS(ref_cell_ghost_long(ref_cell, ref_node, data), "ghost");
    if (ref_cell_n(ref_cell) == 1) {
      REIS(10, data[0], "set ghost");
    }
    ref_free(data);
    ref_node_free(ref_node);
    ref_cell_free(ref_cell);
  }

  RSS(ref_mpi_free(ref_mpi), "cleanup");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
