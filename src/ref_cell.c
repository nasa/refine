
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

#include "ref_cell.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"
#include "ref_sort.h"

static REF_STATUS ref_cell_initialize(REF_CELL ref_cell, REF_CELL_TYPE type) {
  ref_cell_type(ref_cell) = type;

  ref_cell_last_node_is_an_id(ref_cell) = REF_FALSE;
  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
    case REF_CELL_TRI:
    case REF_CELL_QUA:
    case REF_CELL_ED3:
      ref_cell_last_node_is_an_id(ref_cell) = REF_TRUE;
      break;
    case REF_CELL_TET:
    case REF_CELL_PYR:
    case REF_CELL_PRI:
    case REF_CELL_HEX:
      ref_cell_last_node_is_an_id(ref_cell) = REF_FALSE;
      break;
  }

  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
      ref_cell_node_per(ref_cell) = 2;
      break;
    case REF_CELL_ED3:
      ref_cell_node_per(ref_cell) = 4;
      break;
    case REF_CELL_TRI:
      ref_cell_node_per(ref_cell) = 3;
      break;
    case REF_CELL_QUA:
      ref_cell_node_per(ref_cell) = 4;
      break;
    case REF_CELL_TET:
      ref_cell_node_per(ref_cell) = 4;
      break;
    case REF_CELL_PYR:
      ref_cell_node_per(ref_cell) = 5;
      break;
    case REF_CELL_PRI:
      ref_cell_node_per(ref_cell) = 6;
      break;
    case REF_CELL_HEX:
      ref_cell_node_per(ref_cell) = 8;
      break;
    default:
      return REF_IMPLEMENT;
  }

  ref_cell_size_per(ref_cell) = ref_cell_node_per(ref_cell) +
                                (ref_cell_last_node_is_an_id(ref_cell) ? 1 : 0);

  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
      ref_cell_edge_per(ref_cell) = 1;
      break;
    case REF_CELL_ED3:
      ref_cell_edge_per(ref_cell) = 1;
      break;
    case REF_CELL_TRI:
      ref_cell_edge_per(ref_cell) = 3;
      break;
    case REF_CELL_QUA:
      ref_cell_edge_per(ref_cell) = 4;
      break;
    case REF_CELL_TET:
      ref_cell_edge_per(ref_cell) = 6;
      break;
    case REF_CELL_PYR:
      ref_cell_edge_per(ref_cell) = 8;
      break;
    case REF_CELL_PRI:
      ref_cell_edge_per(ref_cell) = 9;
      break;
    case REF_CELL_HEX:
      ref_cell_edge_per(ref_cell) = 12;
      break;
    default:
      ref_cell_edge_per(ref_cell) = 0;
      break;
  }

  ref_cell->e2n = NULL;
  if (ref_cell_edge_per(ref_cell) > 0)
    ref_malloc(ref_cell->e2n, 2 * ref_cell_edge_per(ref_cell), REF_INT);

  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      break;
    case REF_CELL_ED3:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      break;
    case REF_CELL_TRI:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 0;
      break;
    case REF_CELL_QUA:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 3) = 3;
      ref_cell_e2n_gen(ref_cell, 1, 3) = 0;
      break;
    case REF_CELL_TET:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 3) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 3) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 4) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 4) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 5) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 5) = 3;
      break;
    case REF_CELL_PYR:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 3) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 3) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 4) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 4) = 4;
      ref_cell_e2n_gen(ref_cell, 0, 5) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 5) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 6) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 6) = 4;
      ref_cell_e2n_gen(ref_cell, 0, 7) = 3;
      ref_cell_e2n_gen(ref_cell, 1, 7) = 4;
      break;

    case REF_CELL_PRI:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 3) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 3) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 4) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 4) = 4;
      ref_cell_e2n_gen(ref_cell, 0, 5) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 5) = 5;
      ref_cell_e2n_gen(ref_cell, 0, 6) = 3;
      ref_cell_e2n_gen(ref_cell, 1, 6) = 4;
      ref_cell_e2n_gen(ref_cell, 0, 7) = 3;
      ref_cell_e2n_gen(ref_cell, 1, 7) = 5;
      ref_cell_e2n_gen(ref_cell, 0, 8) = 4;
      ref_cell_e2n_gen(ref_cell, 1, 8) = 5;
      break;
    case REF_CELL_HEX:
      ref_cell_e2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_e2n_gen(ref_cell, 0, 1) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 1) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_e2n_gen(ref_cell, 1, 2) = 4;
      ref_cell_e2n_gen(ref_cell, 0, 3) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 3) = 2;
      ref_cell_e2n_gen(ref_cell, 0, 4) = 1;
      ref_cell_e2n_gen(ref_cell, 1, 4) = 5;
      ref_cell_e2n_gen(ref_cell, 0, 5) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 5) = 3;
      ref_cell_e2n_gen(ref_cell, 0, 6) = 2;
      ref_cell_e2n_gen(ref_cell, 1, 6) = 6;
      ref_cell_e2n_gen(ref_cell, 0, 7) = 3;
      ref_cell_e2n_gen(ref_cell, 1, 7) = 7;
      ref_cell_e2n_gen(ref_cell, 0, 8) = 4;
      ref_cell_e2n_gen(ref_cell, 1, 8) = 5;
      ref_cell_e2n_gen(ref_cell, 0, 9) = 4;
      ref_cell_e2n_gen(ref_cell, 1, 9) = 7;
      ref_cell_e2n_gen(ref_cell, 0, 10) = 5;
      ref_cell_e2n_gen(ref_cell, 1, 10) = 6;
      ref_cell_e2n_gen(ref_cell, 0, 11) = 6;
      ref_cell_e2n_gen(ref_cell, 1, 11) = 7;
      break;
  }

  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
      ref_cell_face_per(ref_cell) = 0;
      break;
    case REF_CELL_ED3:
      ref_cell_face_per(ref_cell) = 0;
      break;
    case REF_CELL_TRI:
      ref_cell_face_per(ref_cell) = 1;
      break;
    case REF_CELL_QUA:
      ref_cell_face_per(ref_cell) = 1;
      break;
    case REF_CELL_TET:
      ref_cell_face_per(ref_cell) = 4;
      break;
    case REF_CELL_PYR:
      ref_cell_face_per(ref_cell) = 5;
      break;
    case REF_CELL_PRI:
      ref_cell_face_per(ref_cell) = 5;
      break;
    case REF_CELL_HEX:
      ref_cell_face_per(ref_cell) = 6;
      break;
  }

  ref_cell->f2n = NULL;
  if (ref_cell_face_per(ref_cell) > 0)
    ref_malloc(ref_cell->f2n, 4 * ref_cell_face_per(ref_cell), REF_INT);

  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
    case REF_CELL_ED3:
      break;
    case REF_CELL_TRI:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 0) = 0;
      break;
    case REF_CELL_QUA:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 0) = 3;
      break;
    case REF_CELL_TET:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 1;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 3;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 0) = ref_cell_f2n_gen(ref_cell, 0, 0);
      ref_cell_f2n_gen(ref_cell, 0, 1) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 1) = 2;
      ref_cell_f2n_gen(ref_cell, 2, 1) = 3;
      ref_cell_f2n_gen(ref_cell, 3, 1) = ref_cell_f2n_gen(ref_cell, 0, 1);
      ref_cell_f2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 2) = 3;
      ref_cell_f2n_gen(ref_cell, 2, 2) = 1;
      ref_cell_f2n_gen(ref_cell, 3, 2) = ref_cell_f2n_gen(ref_cell, 0, 2);
      ref_cell_f2n_gen(ref_cell, 0, 3) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 3) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 3) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 3) = ref_cell_f2n_gen(ref_cell, 0, 3);
      break;
    case REF_CELL_PYR:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 0) = ref_cell_f2n_gen(ref_cell, 0, 0);
      ref_cell_f2n_gen(ref_cell, 0, 1) = 1;
      ref_cell_f2n_gen(ref_cell, 1, 1) = 4;
      ref_cell_f2n_gen(ref_cell, 2, 1) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 1) = ref_cell_f2n_gen(ref_cell, 0, 1);
      ref_cell_f2n_gen(ref_cell, 0, 2) = 2;
      ref_cell_f2n_gen(ref_cell, 1, 2) = 4;
      ref_cell_f2n_gen(ref_cell, 2, 2) = 3;
      ref_cell_f2n_gen(ref_cell, 3, 2) = ref_cell_f2n_gen(ref_cell, 0, 2);
      ref_cell_f2n_gen(ref_cell, 0, 3) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 3) = 2;
      ref_cell_f2n_gen(ref_cell, 2, 3) = 3;
      ref_cell_f2n_gen(ref_cell, 3, 3) = ref_cell_f2n_gen(ref_cell, 0, 3);
      ref_cell_f2n_gen(ref_cell, 0, 4) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 4) = 3;
      ref_cell_f2n_gen(ref_cell, 2, 4) = 4;
      ref_cell_f2n_gen(ref_cell, 3, 4) = 1;
      break;
    case REF_CELL_PRI:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 3;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 4;
      ref_cell_f2n_gen(ref_cell, 3, 0) = 1;

      ref_cell_f2n_gen(ref_cell, 0, 1) = 1;
      ref_cell_f2n_gen(ref_cell, 1, 1) = 4;
      ref_cell_f2n_gen(ref_cell, 2, 1) = 5;
      ref_cell_f2n_gen(ref_cell, 3, 1) = 2;

      ref_cell_f2n_gen(ref_cell, 0, 2) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 2) = 2;
      ref_cell_f2n_gen(ref_cell, 2, 2) = 5;
      ref_cell_f2n_gen(ref_cell, 3, 2) = 3;

      ref_cell_f2n_gen(ref_cell, 0, 3) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 3) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 3) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 3) = ref_cell_f2n_gen(ref_cell, 0, 3);

      ref_cell_f2n_gen(ref_cell, 0, 4) = 3;
      ref_cell_f2n_gen(ref_cell, 1, 4) = 5;
      ref_cell_f2n_gen(ref_cell, 2, 4) = 4;
      ref_cell_f2n_gen(ref_cell, 3, 4) = ref_cell_f2n_gen(ref_cell, 0, 4);
      break;
    case REF_CELL_HEX:
      ref_cell_f2n_gen(ref_cell, 0, 0) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 0) = 4;
      ref_cell_f2n_gen(ref_cell, 2, 0) = 5;
      ref_cell_f2n_gen(ref_cell, 3, 0) = 1;

      ref_cell_f2n_gen(ref_cell, 0, 1) = 1;
      ref_cell_f2n_gen(ref_cell, 1, 1) = 5;
      ref_cell_f2n_gen(ref_cell, 2, 1) = 6;
      ref_cell_f2n_gen(ref_cell, 3, 1) = 2;

      ref_cell_f2n_gen(ref_cell, 0, 2) = 2;
      ref_cell_f2n_gen(ref_cell, 1, 2) = 6;
      ref_cell_f2n_gen(ref_cell, 2, 2) = 7;
      ref_cell_f2n_gen(ref_cell, 3, 2) = 3;

      ref_cell_f2n_gen(ref_cell, 0, 3) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 3) = 3;
      ref_cell_f2n_gen(ref_cell, 2, 3) = 7;
      ref_cell_f2n_gen(ref_cell, 3, 3) = 4;

      ref_cell_f2n_gen(ref_cell, 0, 4) = 0;
      ref_cell_f2n_gen(ref_cell, 1, 4) = 1;
      ref_cell_f2n_gen(ref_cell, 2, 4) = 2;
      ref_cell_f2n_gen(ref_cell, 3, 4) = 3;

      ref_cell_f2n_gen(ref_cell, 0, 5) = 4;
      ref_cell_f2n_gen(ref_cell, 1, 5) = 7;
      ref_cell_f2n_gen(ref_cell, 2, 5) = 6;
      ref_cell_f2n_gen(ref_cell, 3, 5) = 5;

      break;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_create(REF_CELL *ref_cell_ptr, REF_CELL_TYPE type) {
  REF_INT cell;
  REF_INT max;
  REF_CELL ref_cell;

  (*ref_cell_ptr) = NULL;

  ref_malloc(*ref_cell_ptr, 1, REF_CELL_STRUCT);

  ref_cell = (*ref_cell_ptr);

  RSS(ref_cell_initialize(ref_cell, type), "init");

  max = 100;

  ref_cell_n(ref_cell) = 0;
  ref_cell_max(ref_cell) = max;

  ref_malloc(ref_cell->c2n,
             ref_cell_max(ref_cell) * ref_cell_size_per(ref_cell), REF_INT);

  for (cell = 0; cell < max; cell++) {
    ref_cell_c2n(ref_cell, 0, cell) = REF_EMPTY;
    ref_cell_c2n(ref_cell, 1, cell) = cell + 1;
  }
  ref_cell_c2n(ref_cell, 1, max - 1) = REF_EMPTY;
  ref_cell_blank(ref_cell) = 0;

  RSS(ref_adj_create(&(ref_cell->ref_adj)), "create ref_adj for ref_cell");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_free(REF_CELL ref_cell) {
  if (NULL == (void *)ref_cell) return REF_NULL;
  ref_adj_free(ref_cell->ref_adj);
  ref_free(ref_cell->c2n);
  ref_free(ref_cell->f2n);
  ref_free(ref_cell->e2n);
  ref_free(ref_cell);
  return REF_SUCCESS;
}

REF_STATUS ref_cell_deep_copy(REF_CELL *ref_cell_ptr, REF_CELL original) {
  REF_INT node, cell;
  REF_INT max;
  REF_CELL ref_cell;
  ref_malloc(*ref_cell_ptr, 1, REF_CELL_STRUCT);

  ref_cell = (*ref_cell_ptr);

  RSS(ref_cell_initialize(ref_cell, ref_cell_type(original)), "init");

  max = ref_cell_max(original);
  ref_cell_n(ref_cell) = ref_cell_n(original);
  ref_cell_max(ref_cell) = max;

  ref_malloc(ref_cell->c2n,
             ref_cell_max(ref_cell) * ref_cell_size_per(ref_cell), REF_INT);
  for (cell = 0; cell < max; cell++)
    for (node = 0; node < ref_cell_size_per(ref_cell); node++)
      ref_cell_c2n(ref_cell, node, cell) = ref_cell_c2n(original, node, cell);

  ref_cell_blank(ref_cell) = ref_cell_blank(original);

  RSS(ref_adj_deep_copy(&(ref_cell->ref_adj), original->ref_adj),
      "deep copy ref_adj for ref_cell");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_pack(REF_CELL ref_cell, REF_INT *o2n) {
  REF_INT node, cell, compact;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  compact = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      ref_cell_c2n(ref_cell, node, compact) = o2n[nodes[node]];
    if (ref_cell_last_node_is_an_id(ref_cell))
      ref_cell_c2n(ref_cell, ref_cell_node_per(ref_cell), compact) =
          nodes[ref_cell_node_per(ref_cell)];
    compact++;
  }
  REIS(compact, ref_cell_n(ref_cell), "count is off");

  if (ref_cell_n(ref_cell) < ref_cell_max(ref_cell)) {
    for (cell = ref_cell_n(ref_cell); cell < ref_cell_max(ref_cell); cell++) {
      ref_cell_c2n(ref_cell, 0, cell) = REF_EMPTY;
      ref_cell_c2n(ref_cell, 1, cell) = cell + 1;
    }
    ref_cell_c2n(ref_cell, 1, ref_cell_max(ref_cell) - 1) = REF_EMPTY;
    ref_cell_blank(ref_cell) = ref_cell_n(ref_cell);
  } else {
    ref_cell_blank(ref_cell) = REF_EMPTY;
  }

  {
    REF_INT *key, *order, *c2n;
    ref_malloc(key, ref_cell_n(ref_cell), REF_INT);
    ref_malloc(order, ref_cell_n(ref_cell), REF_INT);
    ref_malloc(c2n, ref_cell_size_per(ref_cell) * ref_cell_n(ref_cell),
               REF_INT);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      key[cell] = nodes[0];
      for (node = 1; node < ref_cell_node_per(ref_cell); node++) {
        if (nodes[node] < key[cell]) key[cell] = nodes[node];
      }
    }
    RSS(ref_sort_heap_int(ref_cell_n(ref_cell), key, order), "sort smallest");
    for (cell = 0; cell < ref_cell_n(ref_cell); cell++) {
      for (node = 0; node < ref_cell_size_per(ref_cell); node++) {
        c2n[node + cell * ref_cell_size_per(ref_cell)] =
            ref_cell_c2n(ref_cell, node, order[cell]);
      }
    }
    for (cell = 0; cell < ref_cell_n(ref_cell); cell++) {
      for (node = 0; node < ref_cell_size_per(ref_cell); node++) {
        ref_cell_c2n(ref_cell, node, cell) =
            c2n[node + cell * ref_cell_size_per(ref_cell)];
      }
    }
    ref_free(c2n);
    ref_free(order);
    ref_free(key);
  }

  RSS(ref_adj_free(ref_cell_adj(ref_cell)), "free adj");
  RSS(ref_adj_create(&(ref_cell->ref_adj)), "fresh ref_adj for ref_cell");

  for (cell = 0; cell < ref_cell_n(ref_cell); cell++)
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      RSS(ref_adj_add(ref_cell_adj(ref_cell),
                      ref_cell_c2n(ref_cell, node, cell), cell),
          "register cell");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_meshb_keyword(REF_CELL ref_cell, REF_INT *keyword) {
  switch (ref_cell_type(ref_cell)) {
    case REF_CELL_EDG:
      *keyword = 5;
      break;
    case REF_CELL_TRI:
      *keyword = 6;
      break;
    case REF_CELL_QUA:
      *keyword = 7;
      break;
    case REF_CELL_ED3:
      *keyword = 25;
      break;
    case REF_CELL_TET:
      *keyword = 8;
      break;
    case REF_CELL_PYR:
      *keyword = 49;
      break;
    case REF_CELL_PRI:
      *keyword = 9;
      break;
    case REF_CELL_HEX:
      *keyword = 10;
      break;
    default:
      return REF_IMPLEMENT;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cell_inspect(REF_CELL ref_cell) {
  REF_INT cell, node;
  printf("ref_cell = %p\n", (void *)ref_cell);
  printf(" size_per = %d\n", ref_cell_size_per(ref_cell));
  printf(" node_per = %d\n", ref_cell_node_per(ref_cell));
  printf(" n = %d\n", ref_cell_n(ref_cell));
  printf(" max = %d\n", ref_cell_max(ref_cell));
  printf(" blank = %d\n", ref_cell->blank);
  for (cell = 0; cell < ref_cell_max(ref_cell); cell++) {
    if (ref_cell_valid(ref_cell, cell)) {
      printf(" %d:", cell);
      for (node = 0; node < ref_cell_size_per(ref_cell); node++)
        printf(" %d", ref_cell_c2n(ref_cell, node, cell));
      printf("\n");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_tattle(REF_CELL ref_cell, REF_INT cell) {
  REF_INT node;
  printf("cell %d:", cell);
  for (node = 0; node < ref_cell_size_per(ref_cell); node++) {
    printf(" %d", ref_cell_c2n(ref_cell, node, cell));
  }
  printf("\n");
  return REF_SUCCESS;
}

REF_STATUS ref_cell_add(REF_CELL ref_cell, REF_INT *nodes, REF_INT *new_cell) {
  REF_INT node, cell;
  REF_INT orig, chunk;
  REF_INT max_limit = REF_INT_MAX / 4;

  (*new_cell) = REF_EMPTY;

  if (REF_EMPTY == ref_cell_blank(ref_cell)) {
    RAS(ref_cell_max(ref_cell) != max_limit,
        "the number of cells is too large for integers, cannot grow");
    orig = ref_cell_max(ref_cell);
    /* geometric growth for efficiency */
    chunk = MAX(5000, (REF_INT)(1.5 * (REF_DBL)orig));

    /* try to keep under 32-bit limit */
    RAS(max_limit - orig > 0, "chunk limit at max");
    chunk = MIN(chunk, max_limit - orig);

    ref_cell_max(ref_cell) = orig + chunk;

    ref_realloc(ref_cell->c2n,
                ref_cell_size_per(ref_cell) * ref_cell_max(ref_cell), REF_INT);

    for (cell = orig; cell < ref_cell_max(ref_cell); cell++) {
      ref_cell_c2n(ref_cell, 0, cell) = REF_EMPTY;
      ref_cell_c2n(ref_cell, 1, cell) = cell + 1;
    }
    ref_cell_c2n(ref_cell, 1, ref_cell_max(ref_cell) - 1) = REF_EMPTY;
    ref_cell_blank(ref_cell) = orig;
  }

  cell = ref_cell_blank(ref_cell);
  ref_cell_blank(ref_cell) = ref_cell_c2n(ref_cell, 1, cell);
  for (node = 0; node < ref_cell_size_per(ref_cell); node++)
    ref_cell_c2n(ref_cell, node, cell) = nodes[node];

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    RSS(ref_adj_add(ref_cell->ref_adj, nodes[node], cell), "register cell");

  ref_cell_n(ref_cell)++;

  (*new_cell) = cell;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_add_many_global(REF_CELL ref_cell, REF_NODE ref_node,
                                    REF_INT n, REF_GLOB *c2n, REF_INT *part,
                                    REF_INT exclude_part_id) {
  REF_GLOB *global;
  REF_INT nnode;
  REF_INT node, cell;
  REF_INT local, local_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;

  ref_malloc(global, ref_cell_node_per(ref_cell) * n, REF_GLOB);

  nnode = 0;
  for (cell = 0; cell < n; cell++)
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (exclude_part_id != part[node + ref_cell_size_per(ref_cell) * cell]) {
        global[nnode] = c2n[node + ref_cell_size_per(ref_cell) * cell];
        nnode++;
      }
  RSS(ref_node_add_many(ref_node, nnode, global), "many nodes");

  ref_free(global);

  /* set parts */
  for (cell = 0; cell < n; cell++) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      RSS(ref_node_local(
              ref_node, c2n[node + ref_cell_size_per(ref_cell) * cell], &local),
          "local");
      ref_node_part(ref_node, local) =
          part[node + ref_cell_size_per(ref_cell) * cell];
      local_nodes[node] = local;
    }
    if (ref_cell_last_node_is_an_id(ref_cell))
      local_nodes[ref_cell_size_per(ref_cell) - 1] =
          (REF_INT)c2n[(ref_cell_size_per(ref_cell) - 1) +
                       ref_cell_size_per(ref_cell) * cell];

    RXS(ref_cell_with(ref_cell, local_nodes, &new_cell), REF_NOT_FOUND,
        "with failed");

    if (REF_EMPTY == new_cell)
      RSS(ref_cell_add(ref_cell, local_nodes, &new_cell), "add cell");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_remove(REF_CELL ref_cell, REF_INT cell) {
  REF_INT node;
  if (!ref_cell_valid(ref_cell, cell)) return REF_INVALID;
  ref_cell_n(ref_cell)--;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    RSS(ref_adj_remove(ref_cell->ref_adj, ref_cell_c2n(ref_cell, node, cell),
                       cell),
        "unregister cell");

  ref_cell_c2n(ref_cell, 0, cell) = REF_EMPTY;
  ref_cell_c2n(ref_cell, 1, cell) = ref_cell_blank(ref_cell);
  ref_cell_blank(ref_cell) = cell;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_replace_whole(REF_CELL ref_cell, REF_INT cell,
                                  REF_INT *nodes) {
  REF_INT node;
  if (!ref_cell_valid(ref_cell, cell)) return REF_FAILURE;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
    RSS(ref_adj_remove(ref_cell->ref_adj, ref_cell_c2n(ref_cell, node, cell),
                       cell),
        "unregister cell");
    ref_cell_c2n(ref_cell, node, cell) = nodes[node];
    RSS(ref_adj_add(ref_cell->ref_adj, nodes[node], cell),
        "register cell with id");
  }

  if (ref_cell_last_node_is_an_id(ref_cell)) {
    node = ref_cell_size_per(ref_cell) - 1;
    ref_cell_c2n(ref_cell, node, cell) = nodes[node];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_replace_node(REF_CELL ref_cell, REF_INT old_node,
                                 REF_INT new_node) {
  REF_ADJ ref_adj = ref_cell_adj(ref_cell);
  REF_INT node;
  REF_INT item, cell;

  if (old_node == new_node) return REF_SUCCESS;

  item = ref_adj_first(ref_adj, old_node);
  while (ref_adj_valid(item)) {
    cell = ref_adj_item_ref(ref_adj, item);

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (old_node == ref_cell_c2n(ref_cell, node, cell)) {
        RSS(ref_adj_remove(ref_cell->ref_adj,
                           ref_cell_c2n(ref_cell, node, cell), cell),
            "unregister cell");
        ref_cell_c2n(ref_cell, node, cell) = new_node;
        RSS(ref_adj_add(ref_cell->ref_adj, ref_cell_c2n(ref_cell, node, cell),
                        cell),
            "register cell with id");
      }

    item = ref_adj_first(ref_adj, old_node);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_compact(REF_CELL ref_cell, REF_INT **o2n_ptr,
                            REF_INT **n2o_ptr) {
  REF_INT cell;
  REF_INT ncell;
  REF_INT *o2n, *n2o;

  ref_malloc_init(*o2n_ptr, ref_cell_max(ref_cell), REF_INT, REF_EMPTY);
  o2n = *o2n_ptr;
  ref_malloc(*n2o_ptr, ref_cell_n(ref_cell), REF_INT);
  n2o = *n2o_ptr;

  ncell = 0;

  each_ref_cell_valid_cell(ref_cell, cell) {
    o2n[cell] = ncell;
    ncell++;
  }

  RES(ncell, ref_cell_n(ref_cell), "ncell miscount");

  each_ref_cell_valid_cell(ref_cell, cell) n2o[o2n[cell]] = cell;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_nodes(REF_CELL ref_cell, REF_INT cell, REF_INT *nodes) {
  REF_INT node;

  if (!ref_cell_valid(ref_cell, cell)) return REF_INVALID;

  for (node = 0; node < ref_cell_size_per(ref_cell); node++)
    nodes[node] = ref_cell_c2n(ref_cell, node, cell);

  return REF_SUCCESS;
}

REF_STATUS ref_cell_part_cell_node(REF_CELL ref_cell, REF_NODE ref_node,
                                   REF_INT cell, REF_INT *cell_node) {
  REF_GLOB global, smallest_global;
  REF_INT node, smallest_global_node;
  *cell_node = REF_EMPTY;
  if (cell < 0 || cell > ref_cell_max(ref_cell)) return REF_INVALID;
  if (REF_EMPTY == ref_cell_c2n(ref_cell, 0, cell)) return REF_INVALID;
  /* set first node as smallest */
  node = 0;
  global = ref_node_global(ref_node, ref_cell_c2n(ref_cell, node, cell));
  smallest_global = global;
  smallest_global_node = node;
  /* search other nodes for smaller global */
  for (node = 1; node < ref_cell_node_per(ref_cell); node++) {
    global = ref_node_global(ref_node, ref_cell_c2n(ref_cell, node, cell));
    if (global < smallest_global) {
      smallest_global = global;
      smallest_global_node = node;
    }
  }
  RUS(REF_EMPTY, smallest_global_node, "node is empty?");
  *cell_node = smallest_global_node;
  return REF_SUCCESS;
}

REF_STATUS ref_cell_part(REF_CELL ref_cell, REF_NODE ref_node, REF_INT cell,
                         REF_INT *output_part) {
  REF_INT cell_node;
  *output_part = REF_EMPTY;
  RSS(ref_cell_part_cell_node(ref_cell, ref_node, cell, &cell_node),
      "part_cell_node");
  RUS(REF_EMPTY, cell_node, "cell_node is empty?");
  *output_part =
      ref_node_part(ref_node, ref_cell_c2n(ref_cell, cell_node, cell));
  return REF_SUCCESS;
}

REF_STATUS ref_cell_all_local(REF_CELL ref_cell, REF_NODE ref_node,
                              REF_INT cell, REF_BOOL *all_local_nodes) {
  REF_INT cell_node;

  *all_local_nodes = REF_TRUE;
  RAS(ref_cell_valid(ref_cell, cell), "invalid cell");

  each_ref_cell_cell_node(ref_cell, cell_node) {
    *all_local_nodes =
        *all_local_nodes &&
        ref_node_owned(ref_node, ref_cell_c2n(ref_cell, cell_node, cell));
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_local_gem(REF_CELL ref_cell, REF_NODE ref_node,
                              REF_INT node0, REF_INT node1, REF_BOOL *local) {
  REF_INT item, cell, search_node, test_node;

  *local = REF_FALSE;

  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    for (search_node = 0; search_node < ref_cell_node_per(ref_cell);
         search_node++) {
      if (node1 == ref_cell_c2n(ref_cell, search_node, cell)) {
        for (test_node = 0; test_node < ref_cell_node_per(ref_cell);
             test_node++) {
          if (!ref_node_owned(ref_node,
                              ref_cell_c2n(ref_cell, test_node, cell))) {
            return REF_SUCCESS;
          }
        }
      }
    }
  }

  *local = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_ncell(REF_CELL ref_cell, REF_NODE ref_node,
                          REF_LONG *ncell) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, part;

  *ncell = 0;
  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "cell part");
    if (ref_mpi_rank(ref_mpi) == part) {
      (*ncell)++;
    }
  }

  RSS(ref_mpi_allsum(ref_mpi, ncell, 1, REF_LONG_TYPE), "sum");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_orient_node0(REF_INT nnode, REF_INT node0, REF_INT *nodes) {
  REF_INT temp;
  REIS(4, nnode, "only implemented for tets");
  if (node0 == nodes[0]) return REF_SUCCESS;
  if (node0 == nodes[1]) {
    temp = nodes[0];
    nodes[0] = nodes[1];
    nodes[1] = temp;
    temp = nodes[2];
    nodes[2] = nodes[3];
    nodes[3] = temp;
    return REF_SUCCESS;
  }
  if (node0 == nodes[2]) {
    temp = nodes[0];
    nodes[0] = nodes[2];
    nodes[2] = temp;
    temp = nodes[3];
    nodes[3] = nodes[1];
    nodes[1] = temp;
    return REF_SUCCESS;
  }
  if (node0 == nodes[3]) {
    temp = nodes[0];
    nodes[0] = nodes[3];
    nodes[3] = temp;
    temp = nodes[2];
    nodes[2] = nodes[1];
    nodes[1] = temp;
    return REF_SUCCESS;
  }

  RSS(REF_NOT_FOUND, "node0 not found in nodes");
  return REF_SUCCESS;
}

REF_STATUS ref_cell_has_side(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                             REF_BOOL *has_side) {
  REF_INT item, cell;
  REF_INT cell_edge;

  *has_side = REF_FALSE;

  each_ref_adj_node_item_with_ref(ref_cell_adj(ref_cell), node0, item, cell)
      each_ref_cell_cell_edge(
          ref_cell,
          cell_edge) if ((node0 == ref_cell_e2n(ref_cell, 0, cell_edge, cell) &&
                          node1 ==
                              ref_cell_e2n(ref_cell, 1, cell_edge, cell)) ||
                         (node0 == ref_cell_e2n(ref_cell, 1, cell_edge, cell) &&
                          node1 ==
                              ref_cell_e2n(ref_cell, 0, cell_edge, cell))) {
    *has_side = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_side_has_id(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                                REF_INT id, REF_BOOL *has_id) {
  REF_INT item, cell;
  REF_INT cell_edge;

  *has_id = REF_FALSE;

  if (!ref_cell_last_node_is_an_id(ref_cell)) return REF_SUCCESS;

  each_ref_adj_node_item_with_ref(ref_cell_adj(ref_cell), node0, item, cell)
      each_ref_cell_cell_edge(
          ref_cell,
          cell_edge) if ((node0 == ref_cell_e2n(ref_cell, 0, cell_edge, cell) &&
                          node1 ==
                              ref_cell_e2n(ref_cell, 1, cell_edge, cell)) ||
                         (node0 == ref_cell_e2n(ref_cell, 1, cell_edge, cell) &&
                          node1 ==
                              ref_cell_e2n(
                                  ref_cell, 0, cell_edge,
                                  cell))) if (id ==
                                              ref_cell_c2n(
                                                  ref_cell,
                                                  ref_cell_node_per(ref_cell),
                                                  cell)) {
    *has_id = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_id_range(REF_CELL ref_cell, REF_MPI ref_mpi,
                             REF_INT *min_id, REF_INT *max_id) {
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_id = REF_INT_MAX;
  *max_id = REF_INT_MIN;

  RAS(ref_cell_last_node_is_an_id(ref_cell), "cell does not have ids");

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_id = MIN(*min_id, nodes[ref_cell_id_index(ref_cell)]);
    *max_id = MAX(*max_id, nodes[ref_cell_id_index(ref_cell)]);
  }

  if (NULL != ref_mpi && ref_mpi_para(ref_mpi)) {
    REF_INT global;

    RSS(ref_mpi_min(ref_mpi, min_id, &global, REF_INT_TYPE), "mpi min face");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi min face");
    *min_id = global;

    RSS(ref_mpi_max(ref_mpi, max_id, &global, REF_INT_TYPE), "mpi max face");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi max face");
    *max_id = global;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_with_face(REF_CELL ref_cell, REF_INT *face_nodes,
                              REF_INT *cell0, REF_INT *cell1) {
  REF_INT item, node, same, cell_face, cell;
  REF_INT ntarget, target[REF_CELL_MAX_SIZE_PER];
  REF_INT ncandidate, candidate[REF_CELL_MAX_SIZE_PER];
  REF_INT orig[REF_CELL_MAX_SIZE_PER];

  (*cell0) = REF_EMPTY;
  (*cell1) = REF_EMPTY;

  RSS(ref_sort_unique_int(4, face_nodes, &ntarget, target), "t uniq");

  each_ref_cell_having_node(ref_cell, face_nodes[0], item, cell) {
    each_ref_cell_cell_face(ref_cell, cell_face) {
      for (node = 0; node < 4; node++) {
        orig[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
      }
      RSS(ref_sort_unique_int(4, orig, &ncandidate, candidate), "c uniq");

      if (ntarget == ncandidate) {
        same = 0;
        for (node = 0; node < ntarget; node++) {
          if (target[node] == candidate[node]) same++;
        }

        if (ntarget == same) {
          if (REF_EMPTY == *cell0) {
            (*cell0) = cell;
          } else {
            if (REF_EMPTY != *cell1)
              return REF_INVALID; /* more than 2 cells with face */
            (*cell1) = cell;
          }
        }
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_ntri_with_tet_nodes(REF_CELL ref_cell, REF_INT *nodes,
                                        REF_INT *ntri) {
  REF_INT face[3], cell;

  *ntri = 0;

  face[0] = nodes[1];
  face[1] = nodes[2];
  face[2] = nodes[3];
  RXS(ref_cell_with(ref_cell, face, &cell), REF_NOT_FOUND, "t0");
  if (REF_EMPTY != cell) (*ntri)++;

  face[0] = nodes[0];
  face[1] = nodes[2];
  face[2] = nodes[3];
  RXS(ref_cell_with(ref_cell, face, &cell), REF_NOT_FOUND, "t0");
  if (REF_EMPTY != cell) (*ntri)++;

  face[0] = nodes[0];
  face[1] = nodes[3];
  face[2] = nodes[1];
  RXS(ref_cell_with(ref_cell, face, &cell), REF_NOT_FOUND, "t0");
  if (REF_EMPTY != cell) (*ntri)++;

  face[0] = nodes[0];
  face[1] = nodes[1];
  face[2] = nodes[2];
  RXS(ref_cell_with(ref_cell, face, &cell), REF_NOT_FOUND, "t0");
  if (REF_EMPTY != cell) (*ntri)++;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_with(REF_CELL ref_cell, REF_INT *nodes, REF_INT *cell) {
  REF_INT item, ref, node, same;
  REF_INT ntarget, target[REF_CELL_MAX_SIZE_PER];
  REF_INT ncandidate, candidate[REF_CELL_MAX_SIZE_PER];
  REF_INT orig[REF_CELL_MAX_SIZE_PER];

  (*cell) = REF_EMPTY;

  RSS(ref_sort_unique_int(ref_cell_node_per(ref_cell), nodes, &ntarget, target),
      "canonical");

  each_ref_adj_node_item_with_ref(ref_cell_adj(ref_cell), nodes[0], item, ref) {
    RSS(ref_cell_nodes(ref_cell, ref, orig), "get orig");
    RSS(ref_sort_unique_int(ref_cell_node_per(ref_cell), orig, &ncandidate,
                            candidate),
        "canonical");
    if (ntarget == ncandidate) {
      same = 0;
      for (node = 0; node < ntarget; node++)
        if (target[node] == candidate[node]) same++;
      if (ntarget == same) {
        (*cell) = ref;
        return REF_SUCCESS;
      }
    }
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_cell_degree_with2(REF_CELL ref_cell, REF_INT node0,
                                 REF_INT node1, REF_INT *degree) {
  REF_INT item, cell_node, cell;

  *degree = 0;
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    (*degree)++;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_list_with2(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                               REF_INT max_cell, REF_INT *ncell,
                               REF_INT *cell_list) {
  REF_INT cell, item, cell_node;

  *ncell = 0;
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    if (*ncell >= max_cell) return REF_INCREASE_LIMIT;
    cell_list[*ncell] = cell;
    (*ncell)++;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_node_list_around(REF_CELL ref_cell, REF_INT node,
                                     REF_INT max_node, REF_INT *nnode,
                                     REF_INT *node_list) {
  REF_INT cell, item, cell_node, haves;
  REF_BOOL already_have_it;

  *nnode = 0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (node == ref_cell_c2n(ref_cell, cell_node, cell)) continue;
      already_have_it = REF_FALSE;
      for (haves = 0; haves < *nnode; haves++)
        if (node_list[haves] == ref_cell_c2n(ref_cell, cell_node, cell)) {
          already_have_it = REF_TRUE;
          break;
        }
      if (!already_have_it) {
        if (*nnode >= max_node) {
          RSS(REF_INCREASE_LIMIT, "max_node too small");
        }
        node_list[*nnode] = ref_cell_c2n(ref_cell, cell_node, cell);
        (*nnode)++;
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_id_list_around(REF_CELL ref_cell, REF_INT node,
                                   REF_INT max_ids, REF_INT *n_ids,
                                   REF_INT *ids) {
  REF_INT item, cell, id, deg;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL already_have_it;

  *n_ids = 0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    id = nodes[ref_cell_node_per(ref_cell)];
    already_have_it = REF_FALSE;
    for (deg = 0; deg < *n_ids; deg++)
      if (id == ids[deg]) already_have_it = REF_TRUE;
    if (!already_have_it) {
      if (*n_ids >= max_ids) {
        return REF_INCREASE_LIMIT;
      }
      ids[*n_ids] = id;
      (*n_ids)++;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_gen_edge_face(REF_CELL ref_cell, REF_INT edge,
                                  REF_INT *face0, REF_INT *face1) {
  REF_INT face, node0, node1;
  REF_BOOL have_node0;
  REF_BOOL have_node1;

  *face0 = REF_EMPTY;
  *face1 = REF_EMPTY;

  node0 = ref_cell_e2n_gen(ref_cell, 0, edge);
  node1 = ref_cell_e2n_gen(ref_cell, 1, edge);

  for (face = 0; face < ref_cell_face_per(ref_cell); face++) {
    have_node0 = (node0 == ref_cell_f2n_gen(ref_cell, 0, face) ||
                  node0 == ref_cell_f2n_gen(ref_cell, 1, face) ||
                  node0 == ref_cell_f2n_gen(ref_cell, 2, face) ||
                  node0 == ref_cell_f2n_gen(ref_cell, 3, face));
    have_node1 = (node1 == ref_cell_f2n_gen(ref_cell, 0, face) ||
                  node1 == ref_cell_f2n_gen(ref_cell, 1, face) ||
                  node1 == ref_cell_f2n_gen(ref_cell, 2, face) ||
                  node1 == ref_cell_f2n_gen(ref_cell, 3, face));
    if (have_node0 && have_node1) {
      if ((*face0) == REF_EMPTY) {
        (*face0) = face;
      } else {
        RAS(REF_EMPTY == (*face1), "face1 set twice");
        (*face1) = face;
      }
    }
  }

  RAS(REF_EMPTY != (*face0), "face0 not set");
  RAS(REF_EMPTY != (*face1), "face1 not set");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_ghost_long(REF_CELL ref_cell, REF_NODE ref_node,
                               REF_LONG *data) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT cell;
  REF_INT part;
  REF_INT *a_next, *a_cell;
  REF_GLOB *a_nodes, *b_nodes;
  REF_LONG *a_data, *b_data;

  REF_INT cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_GLOB global;
  REF_INT local;
  REF_INT request;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "cell part");
    if (part != ref_mpi_rank(ref_mpi)) a_size[part]++;
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_nodes, ref_cell_node_per(ref_cell) * a_total, REF_GLOB);
  ref_malloc(a_data, a_total, REF_LONG);
  ref_malloc(a_cell, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_nodes, ref_cell_node_per(ref_cell) * b_total, REF_GLOB);
  ref_malloc(b_data, b_total, REF_LONG);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "cell part");
    if (part != ref_mpi_rank(ref_mpi)) {
      a_cell[a_next[part]] = cell;
      for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell);
           cell_node++) {
        a_nodes[cell_node + ref_cell_node_per(ref_cell) * a_next[part]] =
            ref_node_global(ref_node, nodes[cell_node]);
      }
      (a_next[part])++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_nodes, a_size, b_nodes, b_size,
                        ref_cell_node_per(ref_cell), REF_GLOB_TYPE),
      "alltoallv requested nodes");

  for (request = 0; request < b_total; request++) {
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      global = b_nodes[cell_node + ref_cell_node_per(ref_cell) * request];
      RSS(ref_node_local(ref_node, global, &(local)), "local");
      nodes[cell_node] = local;
    }
    RSS(ref_cell_with(ref_cell, nodes, &cell), "find cell");
    b_data[request] = data[cell];
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_data, b_size, a_data, a_size, 1,
                        REF_LONG_TYPE),
      "alltoallv return data");

  for (request = 0; request < a_total; request++) {
    data[a_cell[request]] = a_data[request];
  }

  ref_free(a_next);

  ref_free(b_data);
  ref_free(b_nodes);

  ref_free(a_cell);

  ref_free(a_data);
  ref_free(a_nodes);

  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_cell_global(REF_CELL ref_cell, REF_NODE ref_node,
                           REF_LONG **global) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT ncell, cell, part, *counts, proc;
  REF_LONG offset;

  ref_malloc(*global, ref_cell_max(ref_cell), REF_LONG);

  ncell = 0;
  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "cell part");
    if (ref_mpi_rank(ref_mpi) == part) {
      ncell++;
    }
  }

  ref_malloc(counts, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &ncell, counts, REF_INT_TYPE), "gather size");
  offset = 0;
  for (proc = 0; proc < ref_mpi_rank(ref_mpi); proc++) offset += counts[proc];
  ref_free(counts);

  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "cell part");
    if (ref_mpi_rank(ref_mpi) == part) {
      (*global)[cell] = offset;
      offset++;
    }
  }

  RSS(ref_cell_ghost_long(ref_cell, ref_node, *global), "ghost");

  return REF_SUCCESS;
}

REF_STATUS ref_cell_tec_fill(REF_CELL ref_cell, const char *filename) {
  REF_INT cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine cell fill\"\n");
  fprintf(file, "variables = \"node\" \"cell\"\n");

  fprintf(file, "zone t=\"fill\", i=%d, datapacking=%s\n",
          ref_cell_node_per(ref_cell) * ref_cell_n(ref_cell), "point");

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    each_ref_cell_cell_node(ref_cell, cell_node) {
      fprintf(file, " %d %d\n", ref_cell_c2n(ref_cell, cell_node, cell), cell);
    }
  }

  fclose(file);

  return REF_SUCCESS;
}
