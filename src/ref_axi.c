
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

#include "ref_cell.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_node.h"
#include "ref_sort.h"

REF_STATUS ref_axi_wedge(REF_GRID ref_grid) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_DBL pole_tol;
  REF_INT *o2n, node;
  REF_INT nhalf;

  REF_INT cell, nodes[8], new_nodes[8], pyr_nodes[8], pri_nodes[8];
  REF_INT nunique, unique[8];
  REF_INT new_cell;

  REF_DBL radius, wedge_angle;

  ref_node = ref_grid_node(ref_grid);

  ref_malloc_init(o2n, ref_node_n(ref_node), REF_INT, REF_EMPTY);

  pole_tol = 1.0e-6;
  wedge_angle = ref_math_in_radians(5.0);

  nhalf = ref_node_n(ref_node) / 2;
  each_ref_node_valid_node(ref_node, node) {
    if (ABS(ref_node_xyz(ref_node, 2, node)) < pole_tol &&
        ABS(ref_node_xyz(ref_node, 1, node)) > 0.5) {
      if (node < nhalf) {
        o2n[node] = node + nhalf;
      } else {
        o2n[node] = node - nhalf;
      }
      RSS(ref_node_remove(ref_node, node), "remove");
    } else {
      o2n[node] = node;
    }
    if (ABS(ref_node_xyz(ref_node, 1, node)) > 0.5) {
      radius = ref_node_xyz(ref_node, 2, node);
      ref_node_xyz(ref_node, 1, node) = radius * sin(wedge_angle);
      ref_node_xyz(ref_node, 2, node) = radius * cos(wedge_angle);
    }
  }

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < 3; node++) new_nodes[node] = o2n[nodes[node]];
    new_nodes[3] = nodes[3];
    RSS(ref_cell_replace_whole(ref_cell, cell, new_nodes), "renum");
  }

  ref_cell = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < 4; node++) new_nodes[node] = o2n[nodes[node]];
    new_nodes[4] = nodes[4];
    RSS(ref_sort_unique_int(4, new_nodes, &nunique, unique), "uniq");
    if (4 > nunique) RSS(ref_cell_remove(ref_cell, cell), "rm qua");
    if (3 == nunique) {
      if (new_nodes[2] == new_nodes[3]) {
        new_nodes[3] = new_nodes[4];
      }
      if (new_nodes[1] == new_nodes[2]) {
        new_nodes[2] = new_nodes[3];
        new_nodes[3] = new_nodes[4];
      }
      if (new_nodes[0] == new_nodes[1]) {
        new_nodes[1] = new_nodes[2];
        new_nodes[2] = new_nodes[3];
        new_nodes[3] = new_nodes[4];
      }
      if (new_nodes[3] == new_nodes[0]) {
        new_nodes[0] = new_nodes[1];
        new_nodes[1] = new_nodes[2];
        new_nodes[2] = new_nodes[3];
        new_nodes[3] = new_nodes[4];
      }
      RSS(ref_cell_add(ref_grid_tri(ref_grid), new_nodes, &new_cell),
          "new cell");
    }
    if (4 == nunique) {
      for (node = 0; node < 4; node++) new_nodes[node] = o2n[nodes[node]];
      new_nodes[4] = nodes[4];
      RSS(ref_cell_replace_whole(ref_cell, cell, new_nodes), "renum");
    }
  }

  ref_cell = ref_grid_pri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < 6; node++) new_nodes[node] = o2n[nodes[node]];
    RSS(ref_sort_unique_int(6, new_nodes, &nunique, unique), "uniq");
    if (6 > nunique) RSS(ref_cell_remove(ref_cell, cell), "rm qua");
    if (5 == nunique) {
      if (new_nodes[0] == new_nodes[3]) {
        pyr_nodes[0] = new_nodes[1];
        pyr_nodes[1] = new_nodes[2];
        pyr_nodes[2] = new_nodes[0];
        pyr_nodes[3] = new_nodes[4];
        pyr_nodes[4] = new_nodes[5];
      }
      if (new_nodes[1] == new_nodes[4]) {
        pyr_nodes[0] = new_nodes[2];
        pyr_nodes[1] = new_nodes[0];
        pyr_nodes[2] = new_nodes[1];
        pyr_nodes[3] = new_nodes[5];
        pyr_nodes[4] = new_nodes[3];
      }
      if (new_nodes[2] == new_nodes[5]) {
        pyr_nodes[0] = new_nodes[0];
        pyr_nodes[1] = new_nodes[1];
        pyr_nodes[2] = new_nodes[2];
        pyr_nodes[3] = new_nodes[3];
        pyr_nodes[4] = new_nodes[4];
      }
      RSS(ref_cell_add(ref_grid_pyr(ref_grid), pyr_nodes, &new_cell),
          "new cell");
    }
    if (4 == nunique) {
      if (new_nodes[1] != new_nodes[4]) new_nodes[3] = new_nodes[4];
      if (new_nodes[2] != new_nodes[5]) new_nodes[3] = new_nodes[5];
      RSS(ref_cell_add(ref_grid_tet(ref_grid), new_nodes, &new_cell),
          "new cell");
    }
  }

  ref_cell = ref_grid_hex(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < 8; node++) new_nodes[node] = o2n[nodes[node]];
    RSS(ref_sort_unique_int(8, new_nodes, &nunique, unique), "uniq");
    if (8 == nunique) continue;
    if (8 > nunique) RSS(ref_cell_remove(ref_cell, cell), "rm qua");
    REIS(6, nunique, "expected 6 unique hex nodes on pole");
    if (6 == nunique) {
      pri_nodes[0] = REF_EMPTY;

      /* 0: 0 4 5 1 */
      if (new_nodes[0] == new_nodes[1] && new_nodes[4] == new_nodes[5]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[2];
        pri_nodes[2] = new_nodes[3];
        pri_nodes[3] = new_nodes[4];
        pri_nodes[4] = new_nodes[6];
        pri_nodes[5] = new_nodes[7];
      }
      if (new_nodes[0] == new_nodes[4] && new_nodes[1] == new_nodes[5]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[3];
        pri_nodes[2] = new_nodes[7];
        pri_nodes[3] = new_nodes[1];
        pri_nodes[4] = new_nodes[2];
        pri_nodes[5] = new_nodes[6];
      }

      /* 1: 1 5 6 2 */
      if (new_nodes[1] == new_nodes[2] && new_nodes[5] == new_nodes[6]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[1];
        pri_nodes[2] = new_nodes[3];
        pri_nodes[3] = new_nodes[4];
        pri_nodes[4] = new_nodes[5];
        pri_nodes[5] = new_nodes[7];
      }
      if (new_nodes[1] == new_nodes[5] && new_nodes[2] == new_nodes[6]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[4];
        pri_nodes[2] = new_nodes[1];
        pri_nodes[3] = new_nodes[0];
        pri_nodes[4] = new_nodes[7];
        pri_nodes[5] = new_nodes[2];
      }

      /* 2: 2 6 7 3 */
      if (new_nodes[3] == new_nodes[7] && new_nodes[2] == new_nodes[6]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[3];
        pri_nodes[2] = new_nodes[4];
        pri_nodes[3] = new_nodes[1];
        pri_nodes[4] = new_nodes[2];
        pri_nodes[5] = new_nodes[5];
      }
      if (new_nodes[3] == new_nodes[2] && new_nodes[7] == new_nodes[6]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[1];
        pri_nodes[2] = new_nodes[2];
        pri_nodes[3] = new_nodes[4];
        pri_nodes[4] = new_nodes[5];
        pri_nodes[5] = new_nodes[6];
      }

      /* 3: 0 3 7 4 */
      if (new_nodes[0] == new_nodes[4] && new_nodes[3] == new_nodes[7]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[5];
        pri_nodes[2] = new_nodes[1];
        pri_nodes[3] = new_nodes[3];
        pri_nodes[4] = new_nodes[6];
        pri_nodes[5] = new_nodes[1];
      }
      if (new_nodes[0] == new_nodes[3] && new_nodes[4] == new_nodes[7]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[1];
        pri_nodes[2] = new_nodes[2];
        pri_nodes[3] = new_nodes[4];
        pri_nodes[4] = new_nodes[5];
        pri_nodes[5] = new_nodes[6];
      }

      /* 4: 0 1 2 3 */
      if (new_nodes[0] == new_nodes[1] && new_nodes[3] == new_nodes[2]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[4];
        pri_nodes[2] = new_nodes[5];
        pri_nodes[3] = new_nodes[3];
        pri_nodes[4] = new_nodes[7];
        pri_nodes[5] = new_nodes[6];
      }
      if (new_nodes[0] == new_nodes[3] && new_nodes[1] == new_nodes[2]) {
        pri_nodes[0] = new_nodes[0];
        pri_nodes[1] = new_nodes[7];
        pri_nodes[2] = new_nodes[4];
        pri_nodes[3] = new_nodes[1];
        pri_nodes[4] = new_nodes[6];
        pri_nodes[5] = new_nodes[5];
      }

      /*
      printf(" new %d %d %d %d    %d %d %d %d\n", new_nodes[0], new_nodes[1],
             new_nodes[2], new_nodes[3], new_nodes[4], new_nodes[5],
             new_nodes[6], new_nodes[7]);
             */
      RSS(ref_cell_add(ref_grid_pri(ref_grid), pri_nodes, &new_cell),
          "new cell");
    }
  }

  ref_free(o2n);

  return REF_SUCCESS;
}
