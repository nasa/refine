
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

#include "ref_twod.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

REF_STATUS ref_twod_opposite_node(REF_CELL pri, REF_INT node,
                                  REF_INT *opposite) {
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;

  *opposite = REF_EMPTY;

  cell = ref_cell_first_with(pri, node);

  RSS(ref_cell_nodes(pri, cell, nodes), "get first prism about node");

  if (node == nodes[0]) *opposite = nodes[3];
  if (node == nodes[3]) *opposite = nodes[0];

  if (node == nodes[1]) *opposite = nodes[4];
  if (node == nodes[4]) *opposite = nodes[1];

  if (node == nodes[2]) *opposite = nodes[5];
  if (node == nodes[5]) *opposite = nodes[2];

  return ((REF_EMPTY == (*opposite)) ? REF_NOT_FOUND : REF_SUCCESS);
}

REF_STATUS ref_twod_opposite_edge(REF_CELL pri, REF_INT node0, REF_INT node1,
                                  REF_INT *node2, REF_INT *node3) {
  REF_INT ncell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell_to_split[2];

  *node2 = REF_EMPTY;
  *node3 = REF_EMPTY;

  RSS(ref_cell_list_with2(pri, node0, node1, 2, &ncell, cell_to_split),
      "more than two");

  RAS(ncell > 0, "no cells found");

  RSS(ref_cell_nodes(pri, cell_to_split[0], nodes), "nodes");

  if (node0 == nodes[0] && node1 == nodes[1]) {
    *node2 = nodes[3];
    *node3 = nodes[4];
    return REF_SUCCESS;
  }
  if (node0 == nodes[1] && node1 == nodes[0]) {
    *node2 = nodes[4];
    *node3 = nodes[3];
    return REF_SUCCESS;
  }

  if (node0 == nodes[1] && node1 == nodes[2]) {
    *node2 = nodes[4];
    *node3 = nodes[5];
    return REF_SUCCESS;
  }
  if (node0 == nodes[2] && node1 == nodes[1]) {
    *node2 = nodes[5];
    *node3 = nodes[4];
    return REF_SUCCESS;
  }

  if (node0 == nodes[2] && node1 == nodes[0]) {
    *node2 = nodes[5];
    *node3 = nodes[3];
    return REF_SUCCESS;
  }
  if (node0 == nodes[0] && node1 == nodes[2]) {
    *node2 = nodes[3];
    *node3 = nodes[5];
    return REF_SUCCESS;
  }

  if (node0 == nodes[3] && node1 == nodes[4]) {
    *node2 = nodes[0];
    *node3 = nodes[1];
    return REF_SUCCESS;
  }
  if (node0 == nodes[4] && node1 == nodes[3]) {
    *node2 = nodes[1];
    *node3 = nodes[0];
    return REF_SUCCESS;
  }

  if (node0 == nodes[4] && node1 == nodes[5]) {
    *node2 = nodes[1];
    *node3 = nodes[2];
    return REF_SUCCESS;
  }
  if (node0 == nodes[5] && node1 == nodes[4]) {
    *node2 = nodes[2];
    *node3 = nodes[1];
    return REF_SUCCESS;
  }

  if (node0 == nodes[5] && node1 == nodes[3]) {
    *node2 = nodes[2];
    *node3 = nodes[0];
    return REF_SUCCESS;
  }
  if (node0 == nodes[3] && node1 == nodes[5]) {
    *node2 = nodes[0];
    *node3 = nodes[2];
    return REF_SUCCESS;
  }

  return REF_FAILURE;
}

REF_STATUS ref_twod_tri_pri_tri(REF_CELL tri, REF_CELL pri, REF_INT cell,
                                REF_INT *pri_cell, REF_INT *tri_cell) {
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face[4];
  REF_INT pri_cell2;

  RSS(ref_cell_nodes(tri, cell, tri_nodes), "grab tri");
  face[0] = tri_nodes[0];
  face[1] = tri_nodes[1];
  face[2] = tri_nodes[2];
  face[3] = face[0];
  RSS(ref_cell_with_face(pri, face, pri_cell, &pri_cell2), "pri");

  if (REF_EMPTY == *pri_cell) THROW("prism missing");
  if (REF_EMPTY != pri_cell2) {
    ref_cell_inspect(pri);
    printf("prisms %d %d\n", *pri_cell, pri_cell2);
    printf("face %d %d %d\n", face[0], face[1], face[2]);
    THROW("mulitple prisms found");
  }
  RSS(ref_cell_nodes(pri, *pri_cell, pri_nodes), "grab pri");
  if (tri_nodes[0] == pri_nodes[0] || tri_nodes[0] == pri_nodes[1] ||
      tri_nodes[0] == pri_nodes[2]) {
    tri_nodes[0] = pri_nodes[3];
    tri_nodes[1] = pri_nodes[5];
    tri_nodes[2] = pri_nodes[4];
  } else {
    tri_nodes[0] = pri_nodes[0];
    tri_nodes[1] = pri_nodes[1];
    tri_nodes[2] = pri_nodes[2];
  }

  RSS(ref_cell_with(tri, tri_nodes, tri_cell), "tri");

  return REF_SUCCESS;
}
