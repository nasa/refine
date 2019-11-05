
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
static REF_STATUS ref_dist_pierce(REF_NODE ref_node, REF_INT node0,
                                  REF_INT node1, REF_INT *nodes,
                                  REF_BOOL *pierce) {
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL vol0, vol1;
  REF_DBL side0, side1, side2;
  *pierce = REF_FALSE;
  tet_nodes[0] = nodes[0];
  tet_nodes[1] = nodes[1];
  tet_nodes[2] = nodes[2];
  tet_nodes[3] = node0;
  RSS(ref_node_tet_vol(ref_node, tet_nodes, &vol0), "vol0");
  tet_nodes[3] = node1;
  RSS(ref_node_tet_vol(ref_node, tet_nodes, &vol1), "vol1");
  if ((vol0 > 0.0 && vol1 > 0.0) || (vol0 < 0.0 && vol1 < 0.0)) {
    /* segment above or below triangle */
    return REF_SUCCESS;
  }
  tet_nodes[2] = node0;
  tet_nodes[3] = node1;
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
    *pierce = REF_TRUE;
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_dist_bounding_sphere3(REF_NODE ref_node, REF_INT *nodes,
                                            REF_DBL *center, REF_DBL *radius) {
  REF_INT i;
  for (i = 0; i < 3; i++)
    center[i] = (1.0 / 3.0) * (ref_node_xyz(ref_node, i, nodes[0]) +
                               ref_node_xyz(ref_node, i, nodes[1]) +
                               ref_node_xyz(ref_node, i, nodes[2]));
  *radius = 0.0;
  for (i = 0; i < 3; i++)
    *radius = MAX(
        *radius, sqrt(pow(ref_node_xyz(ref_node, 0, nodes[i]) - center[0], 2) +
                      pow(ref_node_xyz(ref_node, 1, nodes[i]) - center[1], 2) +
                      pow(ref_node_xyz(ref_node, 2, nodes[i]) - center[2], 2)));
  return REF_SUCCESS;
}

static REF_STATUS ref_dist_bounding_sphere2(REF_NODE ref_node, REF_INT node0,
                                            REF_INT node1, REF_DBL *center,
                                            REF_DBL *radius) {
  REF_INT i;
  for (i = 0; i < 3; i++)
    center[i] = 0.5 * (ref_node_xyz(ref_node, i, node0) +
                       ref_node_xyz(ref_node, i, node1));
  *radius = sqrt(pow(ref_node_xyz(ref_node, 0, node0) - center[0], 2) +
                 pow(ref_node_xyz(ref_node, 1, node0) - center[1], 2) +
                 pow(ref_node_xyz(ref_node, 2, node0) - center[2], 2));
  return REF_SUCCESS;
}

REF_STATUS ref_dist_collisions(REF_GRID ref_grid, REF_BOOL report,
                               REF_INT *n_collisions) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_EDGE ref_edge;
  REF_INT item;
  REF_INT edge;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius, scale = 1.01;
  REF_BOOL pierce;
  REF_SEARCH ref_search;
  REF_LIST ref_list;

  *n_collisions = 0;

  RSS(ref_edge_create(&ref_edge, ref_grid), "create edge");
  RSS(ref_search_create(&ref_search, ref_cell_n(ref_cell)), "create search");
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_dist_bounding_sphere3(ref_node, nodes, center, &radius), "b");
    RSS(ref_search_insert(ref_search, cell, center, scale * radius), "insert");
  }

  RSS(ref_list_create(&ref_list), "create list");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_dist_bounding_sphere2(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                  ref_edge_e2n(ref_edge, 1, edge), center,
                                  &radius),
        "b");
    RSS(ref_search_touching(ref_search, ref_list, center, scale * radius),
        "touch");
    each_ref_list_item(ref_list, item) {
      cell = ref_list_value(ref_list, item);
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      if (ref_dist_exclude(ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), nodes))
        continue;
      RSS(ref_dist_pierce(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), nodes, &pierce),
          "hits?");
      if (pierce) {
        (*n_collisions) += 1;
        if (report) {
          printf("%5d face with %f %f %f vertex\n",
                 nodes[ref_cell_id_index(ref_cell)],
                 ref_node_xyz(ref_node, 0, nodes[0]),
                 ref_node_xyz(ref_node, 1, nodes[0]),
                 ref_node_xyz(ref_node, 2, nodes[0]));
        }
      }
    }
    RSS(ref_list_erase(ref_list), "erase");
  }

  RSS(ref_list_free(ref_list), "free");
  RSS(ref_search_free(ref_search), "free");
  RSS(ref_edge_free(ref_edge), "free");
  return REF_SUCCESS;
}
