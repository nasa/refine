
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

#include "ref_comprow.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_edge.h"
#include "ref_grid.h"
#include "ref_malloc.h"
#include "ref_node.h"

REF_STATUS ref_comprow_create(REF_COMPROW *ref_comprow_ptr, REF_GRID ref_grid) {
  REF_COMPROW ref_comprow;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge, node;
  REF_INT n0, n1;

  RSS(ref_edge_create(&ref_edge, ref_grid), "make edges");

  ref_malloc(*ref_comprow_ptr, 1, REF_COMPROW_STRUCT);

  ref_comprow = *ref_comprow_ptr;

  ref_comprow_max(ref_comprow) = ref_node_max(ref_node);
  ref_comprow_nnz(ref_comprow) = 0;
  ref_malloc_init(ref_comprow->first, 1 + ref_node_max(ref_node), REF_INT, 0);
  /* count nodes that the edges touch */
  each_ref_edge(ref_edge, edge) {
    ref_comprow->first[ref_edge_e2n(ref_edge, 0, edge)]++;
    ref_comprow->first[ref_edge_e2n(ref_edge, 1, edge)]++;
  }
  /* add the diagonal */
  for (node = 0; node < ref_node_max(ref_node); node++)
    if (ref_comprow->first[node] > 0) {
      ref_comprow->first[node]++;
    }
  /* cumulative sum to set first to be the last entry on row */
  for (node = 0; node < ref_node_max(ref_node); node++)
    ref_comprow->first[node + 1] += ref_comprow->first[node];

  ref_comprow_nnz(ref_comprow) = ref_comprow->first[ref_node_max(ref_node)];

  ref_malloc_init(ref_comprow->col, ref_comprow_nnz(ref_comprow), REF_INT,
                  REF_EMPTY);

  /* add diagonal */
  for (node = ref_node_max(ref_node); node > 0; node--)
    if (ref_comprow->first[node] > ref_comprow->first[node - 1]) {
      ref_comprow->first[node]--;
      ref_comprow->col[ref_comprow->first[node]] = node;
    }
  node = 0;
  if (ref_comprow->first[node] > 0) {
    ref_comprow->first[node]--;
    ref_comprow->col[ref_comprow->first[node]] = node;
  }
  /* add off-diagonals */
  each_ref_edge(ref_edge, edge) {
    n0 = ref_edge_e2n(ref_edge, 0, edge);
    n1 = ref_edge_e2n(ref_edge, 1, edge);
    ref_comprow->first[n0]--;
    ref_comprow->col[ref_comprow->first[n0]] = n1;
    n0 = ref_edge_e2n(ref_edge, 1, edge);
    n1 = ref_edge_e2n(ref_edge, 0, edge);
    ref_comprow->first[n0]--;
    ref_comprow->col[ref_comprow->first[n0]] = n1;
  }

  RSS(ref_edge_free(ref_edge), "free");

  return REF_SUCCESS;
}

REF_STATUS ref_comprow_free(REF_COMPROW ref_comprow) {
  if (NULL == (void *)ref_comprow) return REF_NULL;

  ref_free(ref_comprow->col);
  ref_free(ref_comprow->first);

  ref_free(ref_comprow);

  return REF_SUCCESS;
}

REF_STATUS ref_comprow_inspect(REF_COMPROW ref_comprow) {
  REF_INT node, entry;
  for (node = 0; node < ref_comprow_max(ref_comprow); node++)
    if (ref_comprow->first[node + 1] > ref_comprow->first[node]) {
      printf(" row %d (%d) :", node,
             ref_comprow->first[node + 1] - ref_comprow->first[node]);
      for (entry = ref_comprow->first[node];
           entry < ref_comprow->first[node + 1]; entry++) {
        printf(" %d", ref_comprow->col[entry]);
      }
      printf("\n");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_comprow_entry(REF_COMPROW ref_comprow, REF_INT row, REF_INT col,
                             REF_INT *entry) {
  (*entry) = REF_EMPTY;
  if (row < 0 || ref_comprow_max(ref_comprow) <= row || col < 0 ||
      ref_comprow_max(ref_comprow) <= col)
    return REF_INVALID;
  if (ref_comprow->first[row + 1] == ref_comprow->first[row])
    return REF_NOT_FOUND;
  for ((*entry) = ref_comprow->first[row];
       (*entry) < ref_comprow->first[row + 1]; (*entry)++) {
    if (col == ref_comprow->col[*entry]) return REF_SUCCESS;
  }
  (*entry) = REF_EMPTY;
  return REF_NOT_FOUND;
}
