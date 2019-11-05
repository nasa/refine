
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

#include <stdio.h>
#include <stdlib.h>

#include "ref_dist.h"

#include "ref_cell.h"
#include "ref_grid.h"

#include "ref_edge.h"

static REF_BOOL ref_dist_exclude(REF_INT node0, REF_INT node1, REF_INT *nodes) {
  if (node0 == nodes[0] && node1 == nodes[1]) return REF_TRUE;
  if (node0 == nodes[1] && node1 == nodes[2]) return REF_TRUE;
  if (node0 == nodes[2] && node1 == nodes[0]) return REF_TRUE;

  if (node1 == nodes[0] && node0 == nodes[1]) return REF_TRUE;
  if (node1 == nodes[1] && node0 == nodes[2]) return REF_TRUE;
  if (node1 == nodes[2] && node0 == nodes[0]) return REF_TRUE;

  return REF_FALSE;
}

REF_STATUS ref_dist_collisions(REF_GRID ref_grid, REF_BOOL report,
                               REF_INT *n_collisions) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL vol0, vol1;
  REF_DBL side0, side1, side2;

  *n_collisions = 0;
  RSS(ref_edge_create(&ref_edge, ref_grid), "create edge");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (ref_dist_exclude(ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), nodes))
        continue;
      tet_nodes[0] = nodes[0];
      tet_nodes[1] = nodes[1];
      tet_nodes[2] = nodes[2];
      tet_nodes[3] = ref_edge_e2n(ref_edge, 0, edge);
      RSS(ref_node_tet_vol(ref_node, tet_nodes, &vol0), "vol0");
      tet_nodes[3] = ref_edge_e2n(ref_edge, 1, edge);
      RSS(ref_node_tet_vol(ref_node, tet_nodes, &vol1), "vol1");
      if ((vol0 >= 0.0 && vol1 <= 0.0) ||
          (vol0 <= 0.0 && vol1 >= 0.0)) { /* canidate, above and below */
        tet_nodes[2] = ref_edge_e2n(ref_edge, 0, edge);
        tet_nodes[3] = ref_edge_e2n(ref_edge, 1, edge);
        tet_nodes[0] = nodes[1];
        tet_nodes[1] = nodes[2];
        RSS(ref_node_tet_vol(ref_node, tet_nodes, &side0), "side0");
        tet_nodes[0] = nodes[2];
        tet_nodes[1] = nodes[0];
        RSS(ref_node_tet_vol(ref_node, tet_nodes, &side1), "side1");
        tet_nodes[0] = nodes[0];
        tet_nodes[1] = nodes[1];
        RSS(ref_node_tet_vol(ref_node, tet_nodes, &side2), "side2");
        if ((side0 > 0.0 && side1 > 0.0 && side2 > 0.0) ||
            (side0 < 0.0 && side1 < 0.0 && side2 < 0.0)) { /* inside */
          (*n_collisions) += 1;
          if (report) {
            RSS(ref_node_location(ref_node, ref_edge_e2n(ref_edge, 0, edge)),
                "edge0");
            RSS(ref_node_location(ref_node, ref_edge_e2n(ref_edge, 1, edge)),
                "edge1");
          }
        }
      }
    }
  }

  RSS(ref_edge_free(ref_edge), "free");
  return REF_SUCCESS;
}
