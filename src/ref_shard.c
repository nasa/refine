
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

#include "ref_shard.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_export.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

REF_STATUS ref_shard_create(REF_SHARD *ref_shard_ptr, REF_GRID ref_grid) {
  REF_SHARD ref_shard;
  REF_INT face;

  ref_malloc(*ref_shard_ptr, 1, REF_SHARD_STRUCT);

  ref_shard = *ref_shard_ptr;

  ref_shard_grid(ref_shard) = ref_grid;

  RSS(ref_face_create(&(ref_shard_face(ref_shard)), ref_shard_grid(ref_shard)),
      "create face");

  ref_malloc(ref_shard->mark, ref_face_n(ref_shard_face(ref_shard)), REF_INT);

  for (face = 0; face < ref_face_n(ref_shard_face(ref_shard)); face++)
    ref_shard_mark(ref_shard, face) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_shard_free(REF_SHARD ref_shard) {
  if (NULL == (void *)ref_shard) return REF_NULL;

  ref_free(ref_shard->mark);
  RSS(ref_face_free(ref_shard_face(ref_shard)), "free face");

  ref_free(ref_shard);

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_to_split(REF_SHARD ref_shard, REF_INT node0,
                                   REF_INT node1) {
  REF_INT face;

  RSS(ref_face_spanning(ref_shard_face(ref_shard), node0, node1, &face),
      "missing face");

  if (node0 == ref_face_f2n(ref_shard_face(ref_shard), 2, face) ||
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 2, face)) {
    if (3 == ref_shard_mark(ref_shard, face))
      RSS(REF_FAILURE, "2-3 mark mismatch");
    ref_shard_mark(ref_shard, face) = 2;
    return REF_SUCCESS;
  }

  if (node0 == ref_face_f2n(ref_shard_face(ref_shard), 3, face) ||
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 3, face)) {
    if (2 == ref_shard_mark(ref_shard, face))
      RSS(REF_FAILURE, "3-2 mark mismatch");
    ref_shard_mark(ref_shard, face) = 3;
    return REF_SUCCESS;
  }

  return REF_FAILURE;
}

REF_STATUS ref_shard_marked(REF_SHARD ref_shard, REF_INT node0, REF_INT node1,
                            REF_BOOL *marked) {
  REF_INT face;

  *marked = REF_FALSE;

  RSS(ref_face_spanning(ref_shard_face(ref_shard), node0, node1, &face),
      "missing face");

  if (0 == ref_shard_mark(ref_shard, face)) return REF_SUCCESS;

  if (2 == ref_shard_mark(ref_shard, face) &&
      node0 == ref_face_f2n(ref_shard_face(ref_shard), 2, face) &&
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 0, face))
    *marked = REF_TRUE;

  if (2 == ref_shard_mark(ref_shard, face) &&
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 2, face) &&
      node0 == ref_face_f2n(ref_shard_face(ref_shard), 0, face))
    *marked = REF_TRUE;

  if (3 == ref_shard_mark(ref_shard, face) &&
      node0 == ref_face_f2n(ref_shard_face(ref_shard), 3, face) &&
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 1, face))
    *marked = REF_TRUE;

  if (3 == ref_shard_mark(ref_shard, face) &&
      node1 == ref_face_f2n(ref_shard_face(ref_shard), 3, face) &&
      node0 == ref_face_f2n(ref_shard_face(ref_shard), 1, face))
    *marked = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_n(REF_SHARD ref_shard, REF_INT *face_marks,
                            REF_INT *hex_marks) {
  REF_INT face;
  REF_INT cell, hex_nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL marked14, marked05;
  REF_BOOL marked16, marked25;
  REF_BOOL marked02, marked13;

  (*face_marks) = 0;

  for (face = 0; face < ref_face_n(ref_shard_face(ref_shard)); face++)
    if (0 != ref_shard_mark(ref_shard, face)) (*face_marks)++;

  (*hex_marks) = 0;
  each_ref_cell_valid_cell_with_nodes(ref_grid_hex(ref_shard_grid(ref_shard)),
                                      cell, hex_nodes) {
    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[4], &marked14),
        "1-4");
    RSS(ref_shard_marked(ref_shard, hex_nodes[0], hex_nodes[5], &marked05),
        "0-5");

    if (marked14 || marked05) (*hex_marks)++;

    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[6], &marked16),
        "1-6");
    RSS(ref_shard_marked(ref_shard, hex_nodes[2], hex_nodes[5], &marked25),
        "2-5");

    if (marked16 || marked25) (*hex_marks)++;

    RSS(ref_shard_marked(ref_shard, hex_nodes[0], hex_nodes[2], &marked02),
        "0-2");
    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[3], &marked13),
        "1-3");

    if (marked02 || marked13) (*hex_marks)++;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_cell_edge_split(REF_SHARD ref_shard, REF_INT cell,
                                          REF_INT cell_edge) {
  REF_CELL ref_cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_cell = ref_grid_hex(ref_shard_grid(ref_shard));

  RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");

  switch (cell_edge) {
    case 0:
    case 11:
      RSS(ref_shard_mark_to_split(ref_shard, nodes[1], nodes[6]), "mark1");
      RSS(ref_shard_mark_to_split(ref_shard, nodes[0], nodes[7]), "mark2");
      break;
    case 5:
    case 8:
      RSS(ref_shard_mark_to_split(ref_shard, nodes[2], nodes[5]), "mark1");
      RSS(ref_shard_mark_to_split(ref_shard, nodes[3], nodes[4]), "mark2");
      break;

      /*
    case 2:  case 6:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[4], nodes[6] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[0], nodes[2] ), "mark2" );
      break;
    case 4:  case 7:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[7], nodes[5] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[3], nodes[1] ), "mark2" );
      break;
      */

    default:
      printf("%s: %d: %s: cell edge %d not implemented, skipping.\n", __FILE__,
             __LINE__, __func__, cell_edge);
      /*
      RSB( REF_IMPLEMENT, "can not handle cell edge",
           printf("cell edge %d\n",cell_edge););
      */
      break;
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_shard_pair(REF_SHARD ref_shard, REF_BOOL *again,
                                 REF_INT a0, REF_INT a1, REF_INT b0,
                                 REF_INT b1) {
  REF_BOOL a_marked, b_marked;

  RSS(ref_shard_marked(ref_shard, a0, a1, &a_marked), "marked? a0-a1");
  RSS(ref_shard_marked(ref_shard, b0, b1, &b_marked), "marked? b0-b1");
  if (a_marked != b_marked) {
    *again = REF_TRUE;
    RSS(ref_shard_mark_to_split(ref_shard, a0, a1), "mark a0-a1");
    RSS(ref_shard_mark_to_split(ref_shard, b0, b1), "mark b0-b1");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_relax(REF_SHARD ref_shard) {
  REF_CELL ref_cell;
  REF_BOOL again;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_cell = ref_grid_hex(ref_shard_grid(ref_shard));

  again = REF_TRUE;

  while (again) {
    again = REF_FALSE;

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_shard_pair(ref_shard, &again, nodes[1], nodes[4], nodes[2],
                         nodes[7]),
          "not consist");
      RSS(ref_shard_pair(ref_shard, &again, nodes[5], nodes[0], nodes[6],
                         nodes[3]),
          "not consist");

      RSS(ref_shard_pair(ref_shard, &again, nodes[1], nodes[6], nodes[0],
                         nodes[7]),
          "not consist");
      RSS(ref_shard_pair(ref_shard, &again, nodes[5], nodes[2], nodes[4],
                         nodes[3]),
          "not consist");

      RSS(ref_shard_pair(ref_shard, &again, nodes[0], nodes[2], nodes[4],
                         nodes[6]),
          "not consist");
      RSS(ref_shard_pair(ref_shard, &again, nodes[1], nodes[3], nodes[5],
                         nodes[7]),
          "not consist");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_split(REF_SHARD ref_shard) {
  REF_GRID ref_grid;
  REF_CELL hex, pri, tri, qua;
  REF_INT cell, hex_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pri_nodes[6], new_cell;
  REF_INT tri_nodes[4], qua_nodes[5];
  REF_BOOL marked;

  ref_grid = ref_shard_grid(ref_shard);
  hex = ref_grid_hex(ref_grid);
  pri = ref_grid_pri(ref_grid);
  tri = ref_grid_tri(ref_grid);
  qua = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes(qua, cell, qua_nodes) {
    RSS(ref_shard_marked(ref_shard, qua_nodes[0], qua_nodes[2], &marked),
        "0-2");
    if (marked) {
      RSS(ref_cell_remove(qua, cell), "remove qua");
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[2];
      tri_nodes[3] = qua_nodes[4]; /* bound id */
      RSS(ref_cell_add(tri, tri_nodes, &new_cell), "add tri");
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[2];
      tri_nodes[2] = qua_nodes[3];
      tri_nodes[3] = qua_nodes[4]; /* bound id */
      RSS(ref_cell_add(tri, tri_nodes, &new_cell), "add tri");
    }
    RSS(ref_shard_marked(ref_shard, qua_nodes[1], qua_nodes[3], &marked),
        "1-3");
    if (marked) {
      RSS(ref_cell_remove(qua, cell), "remove qua");
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[3];
      tri_nodes[3] = qua_nodes[4]; /* bound id */
      RSS(ref_cell_add(tri, tri_nodes, &new_cell), "add tri");
      tri_nodes[0] = qua_nodes[1];
      tri_nodes[1] = qua_nodes[2];
      tri_nodes[2] = qua_nodes[3];
      tri_nodes[3] = qua_nodes[4]; /* bound id */
      RSS(ref_cell_add(tri, tri_nodes, &new_cell), "add tri");
    }
  }

  each_ref_cell_valid_cell_with_nodes(hex, cell, hex_nodes) {
    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[4], &marked),
        "1-4");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[0];
      pri_nodes[2] = hex_nodes[4];
      pri_nodes[3] = hex_nodes[2];
      pri_nodes[4] = hex_nodes[3];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[4];
      pri_nodes[2] = hex_nodes[5];
      pri_nodes[3] = hex_nodes[2];
      pri_nodes[4] = hex_nodes[7];
      pri_nodes[5] = hex_nodes[6];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex_pri 2");
      continue;
    }
    RSS(ref_shard_marked(ref_shard, hex_nodes[0], hex_nodes[5], &marked),
        "0-5");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[0];
      pri_nodes[1] = hex_nodes[5];
      pri_nodes[2] = hex_nodes[1];
      pri_nodes[3] = hex_nodes[3];
      pri_nodes[4] = hex_nodes[6];
      pri_nodes[5] = hex_nodes[2];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[4];
      pri_nodes[1] = hex_nodes[5];
      pri_nodes[2] = hex_nodes[0];
      pri_nodes[3] = hex_nodes[7];
      pri_nodes[4] = hex_nodes[6];
      pri_nodes[5] = hex_nodes[3];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex_pri 2");
      continue;
    }

    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[6], &marked),
        "1-6");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[5];
      pri_nodes[2] = hex_nodes[6];
      pri_nodes[3] = hex_nodes[0];
      pri_nodes[4] = hex_nodes[4];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[6];
      pri_nodes[2] = hex_nodes[2];
      pri_nodes[3] = hex_nodes[0];
      pri_nodes[4] = hex_nodes[7];
      pri_nodes[5] = hex_nodes[3];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex_pri 2");
      continue;
    }
    RSS(ref_shard_marked(ref_shard, hex_nodes[2], hex_nodes[5], &marked),
        "2-5");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[5];
      pri_nodes[2] = hex_nodes[2];
      pri_nodes[3] = hex_nodes[0];
      pri_nodes[4] = hex_nodes[4];
      pri_nodes[5] = hex_nodes[3];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[2];
      pri_nodes[1] = hex_nodes[5];
      pri_nodes[2] = hex_nodes[6];
      pri_nodes[3] = hex_nodes[3];
      pri_nodes[4] = hex_nodes[4];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 2");
      continue;
    }

    RSS(ref_shard_marked(ref_shard, hex_nodes[0], hex_nodes[2], &marked),
        "0-2");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[0];
      pri_nodes[1] = hex_nodes[1];
      pri_nodes[2] = hex_nodes[2];
      pri_nodes[3] = hex_nodes[4];
      pri_nodes[4] = hex_nodes[5];
      pri_nodes[5] = hex_nodes[6];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[0];
      pri_nodes[1] = hex_nodes[2];
      pri_nodes[2] = hex_nodes[3];
      pri_nodes[3] = hex_nodes[4];
      pri_nodes[4] = hex_nodes[6];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex_pri 2");
      continue;
    }

    RSS(ref_shard_marked(ref_shard, hex_nodes[1], hex_nodes[3], &marked),
        "1-3");
    if (marked) {
      RSS(ref_cell_remove(hex, cell), "remove hex");
      pri_nodes[0] = hex_nodes[0];
      pri_nodes[1] = hex_nodes[1];
      pri_nodes[2] = hex_nodes[3];
      pri_nodes[3] = hex_nodes[4];
      pri_nodes[4] = hex_nodes[5];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex pri 1");
      pri_nodes[0] = hex_nodes[1];
      pri_nodes[1] = hex_nodes[2];
      pri_nodes[2] = hex_nodes[3];
      pri_nodes[3] = hex_nodes[5];
      pri_nodes[4] = hex_nodes[6];
      pri_nodes[5] = hex_nodes[7];
      RSS(ref_cell_add(pri, pri_nodes, &new_cell), "add hex_pri 2");
      continue;
    }
  }

  return REF_SUCCESS;
}

#define check_tet_volume()                                                   \
  {                                                                          \
    REF_DBL vol;                                                             \
    RSS(ref_node_tet_vol(ref_node, tet_nodes, &vol), "tet vol");             \
    if (vol <= 0.0) {                                                        \
      printf("tet vol %e\n", vol);                                           \
      printf("minnode " REF_GLOB_FMT "\n", minnode);                         \
      printf("orig %d %d %d %d %d %d\n", orig[0], orig[1], orig[2], orig[3], \
             orig[4], orig[5]);                                              \
      printf("prism %d %d %d %d %d %d\n", pri_nodes[0], pri_nodes[1],        \
             pri_nodes[2], pri_nodes[3], pri_nodes[4], pri_nodes[5]);        \
      printf("tet %d %d %d %d\n", tet_nodes[0], tet_nodes[1], tet_nodes[2],  \
             tet_nodes[3]);                                                  \
    }                                                                        \
  }

static REF_STATUS ref_shard_cell_add_local(REF_NODE ref_node, REF_CELL ref_cell,
                                           REF_INT *nodes) {
  REF_BOOL has_local;
  REF_INT node, new_cell;

  has_local = REF_FALSE;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    has_local = has_local || (ref_mpi_rank(ref_node_mpi(ref_node)) ==
                              ref_node_part(ref_node, nodes[node]));

  if (has_local) RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add");

  return REF_SUCCESS;
}

REF_STATUS ref_shard_prism_into_tet(REF_GRID ref_grid, REF_INT keeping_n_layers,
                                    REF_INT of_faceid) {
  REF_INT cell, tri_mark;
  REF_GLOB minnode;

  REF_INT orig[REF_CELL_MAX_SIZE_PER];
  REF_GLOB global[REF_CELL_MAX_SIZE_PER];
  REF_INT pri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_GLOB pri_global[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_CELL pyr = ref_grid_pyr(ref_grid);
  REF_CELL tet = ref_grid_tet(ref_grid);

  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua_nodes[REF_CELL_MAX_SIZE_PER];
  REF_GLOB qua_global[REF_CELL_MAX_SIZE_PER];
  REF_CELL qua = ref_grid_qua(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);

  REF_NODE ref_node = ref_grid_node(ref_grid);

  REF_INT relaxation;
  REF_INT node;
  REF_INT *mark, *mark_copy;

  ref_malloc_init(mark, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  ref_malloc_init(mark_copy, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  /* mark nodes on prism tris */
  if (0 < keeping_n_layers)
    each_ref_cell_valid_cell_with_nodes(pri, cell, orig) {
      tri_nodes[0] = orig[0];
      tri_nodes[1] = orig[1];
      tri_nodes[2] = orig[2];
      RXS(ref_cell_with(tri, tri_nodes, &tri_mark), REF_NOT_FOUND, "with");
      if (REF_EMPTY != tri_mark &&
          (of_faceid == REF_EMPTY ||
           of_faceid == ref_cell_c2n(tri, 3, tri_mark))) {
        mark[tri_nodes[0]] = 0;
        mark[tri_nodes[1]] = 0;
        mark[tri_nodes[2]] = 0;
      }
      tri_nodes[0] = orig[3];
      tri_nodes[1] = orig[5];
      tri_nodes[2] = orig[4];
      RXS(ref_cell_with(tri, tri_nodes, &tri_mark), REF_NOT_FOUND, "with");
      if (REF_EMPTY != tri_mark &&
          (of_faceid == REF_EMPTY ||
           of_faceid == ref_cell_c2n(tri, 3, tri_mark))) {
        mark[tri_nodes[0]] = 0;
        mark[tri_nodes[1]] = 0;
        mark[tri_nodes[2]] = 0;
      }
    }
  RSS(ref_node_ghost_int(ref_node, mark, 1), "update ghost mark");

  for (relaxation = 0; relaxation < keeping_n_layers; relaxation++) {
    for (node = 0; node < ref_node_max(ref_node); node++)
      mark_copy[node] = mark[node];
    each_ref_cell_valid_cell_with_nodes(pri, cell, orig) {
      if (mark_copy[orig[0]] == REF_EMPTY && mark_copy[orig[3]] != REF_EMPTY)
        mark[orig[0]] = mark_copy[orig[3]] + 1;
      if (mark_copy[orig[3]] == REF_EMPTY && mark_copy[orig[0]] != REF_EMPTY)
        mark[orig[3]] = mark_copy[orig[0]] + 1;

      if (mark_copy[orig[1]] == REF_EMPTY && mark_copy[orig[4]] != REF_EMPTY)
        mark[orig[1]] = mark_copy[orig[4]] + 1;
      if (mark_copy[orig[4]] == REF_EMPTY && mark_copy[orig[1]] != REF_EMPTY)
        mark[orig[4]] = mark_copy[orig[1]] + 1;

      if (mark_copy[orig[2]] == REF_EMPTY && mark_copy[orig[5]] != REF_EMPTY)
        mark[orig[2]] = mark_copy[orig[5]] + 1;
      if (mark_copy[orig[5]] == REF_EMPTY && mark_copy[orig[2]] != REF_EMPTY)
        mark[orig[5]] = mark_copy[orig[2]] + 1;
    }
    RSS(ref_node_ghost_int(ref_node, mark, 1), "update ghost mark");
  }

  each_ref_cell_valid_cell_with_nodes(pri, cell, orig) {
    if (mark[orig[0]] != REF_EMPTY && mark[orig[3]] != REF_EMPTY) continue;

    RSS(ref_cell_remove(pri, cell), "remove pri");

    for (node = 0; node < ref_cell_node_per(pri); node++)
      global[node] = ref_node_global(ref_node, orig[node]);

    minnode = MIN(MIN(global[0], global[1]), MIN(global[2], global[3]));
    minnode = MIN(MIN(global[4], global[5]), minnode);

    pri_nodes[0] = orig[0];
    pri_nodes[1] = orig[1];
    pri_nodes[2] = orig[2];
    pri_nodes[3] = orig[3];
    pri_nodes[4] = orig[4];
    pri_nodes[5] = orig[5];

    if (global[1] == minnode) {
      pri_nodes[0] = orig[1];
      pri_nodes[1] = orig[2];
      pri_nodes[2] = orig[0];
      pri_nodes[3] = orig[4];
      pri_nodes[4] = orig[5];
      pri_nodes[5] = orig[3];
    }

    if (global[2] == minnode) {
      pri_nodes[0] = orig[2];
      pri_nodes[1] = orig[0];
      pri_nodes[2] = orig[1];
      pri_nodes[3] = orig[5];
      pri_nodes[4] = orig[3];
      pri_nodes[5] = orig[4];
    }

    if (global[3] == minnode) {
      pri_nodes[0] = orig[3];
      pri_nodes[1] = orig[5];
      pri_nodes[2] = orig[4];
      pri_nodes[3] = orig[0];
      pri_nodes[4] = orig[2];
      pri_nodes[5] = orig[1];
    }

    if (global[4] == minnode) {
      pri_nodes[0] = orig[4];
      pri_nodes[1] = orig[3];
      pri_nodes[2] = orig[5];
      pri_nodes[3] = orig[1];
      pri_nodes[4] = orig[0];
      pri_nodes[5] = orig[2];
    }

    if (global[5] == minnode) {
      pri_nodes[0] = orig[5];
      pri_nodes[1] = orig[4];
      pri_nodes[2] = orig[3];
      pri_nodes[3] = orig[2];
      pri_nodes[4] = orig[1];
      pri_nodes[5] = orig[0];
    }

    /* node 0 is now the smallest global index of prism */

    tet_nodes[0] = pri_nodes[0];
    tet_nodes[1] = pri_nodes[4];
    tet_nodes[2] = pri_nodes[5];
    tet_nodes[3] = pri_nodes[3];
    RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "add tet");
    check_tet_volume();

    for (node = 0; node < ref_cell_node_per(pri); node++)
      pri_global[node] = ref_node_global(ref_node, pri_nodes[node]);

    if ((pri_global[1] < pri_global[2] && pri_global[1] < pri_global[4]) ||
        (pri_global[5] < pri_global[2] && pri_global[5] < pri_global[4])) {
      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[4];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      check_tet_volume();

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[2];
      tet_nodes[3] = pri_nodes[5];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      check_tet_volume();
    } else {
      tet_nodes[0] = pri_nodes[2];
      tet_nodes[1] = pri_nodes[0];
      tet_nodes[2] = pri_nodes[4];
      tet_nodes[3] = pri_nodes[5];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      check_tet_volume();

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[2];
      tet_nodes[3] = pri_nodes[4];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      check_tet_volume();
    }
  }

  each_ref_cell_valid_cell_with_nodes(pyr, cell, orig) {
    if (mark[orig[0]] != REF_EMPTY && mark[orig[1]] != REF_EMPTY &&
        mark[orig[3]] != REF_EMPTY && mark[orig[4]] != REF_EMPTY)
      continue;

    for (node = 0; node < ref_cell_node_per(pyr); node++)
      global[node] = ref_node_global(ref_node, orig[node]);

    RSS(ref_cell_remove(pyr, cell), "remove qua");
    if ((global[0] < global[1] && global[0] < global[3]) ||
        (global[4] < global[1] &&
         global[4] < global[3])) { /* 0-4 diag split of quad */
                                   /* 4-1\
                                      |\| 2
                                      3-0/ */
      tet_nodes[0] = orig[0];
      tet_nodes[1] = orig[4];
      tet_nodes[2] = orig[1];
      tet_nodes[3] = orig[2];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      tet_nodes[0] = orig[0];
      tet_nodes[1] = orig[3];
      tet_nodes[2] = orig[4];
      tet_nodes[3] = orig[2];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
    } else { /* 3-1 diag split of quad */
             /* 4-1\
                |/| 2
                3-0/ */
      tet_nodes[0] = orig[0];
      tet_nodes[1] = orig[3];
      tet_nodes[2] = orig[1];
      tet_nodes[3] = orig[2];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
      tet_nodes[0] = orig[1];
      tet_nodes[1] = orig[3];
      tet_nodes[2] = orig[4];
      tet_nodes[3] = orig[2];
      RSS(ref_shard_cell_add_local(ref_node, tet, tet_nodes), "a tet");
    }
  }

  each_ref_cell_valid_cell_with_nodes(qua, cell, qua_nodes) {
    if (mark[qua_nodes[0]] != REF_EMPTY && mark[qua_nodes[1]] != REF_EMPTY &&
        mark[qua_nodes[2]] != REF_EMPTY && mark[qua_nodes[3]] != REF_EMPTY)
      continue;

    RSS(ref_cell_remove(qua, cell), "remove qua");
    tri_nodes[3] = qua_nodes[4]; /* patch id */

    for (node = 0; node < ref_cell_node_per(qua); node++)
      qua_global[node] = ref_node_global(ref_node, qua_nodes[node]);

    if ((qua_global[0] < qua_global[1] && qua_global[0] < qua_global[3]) ||
        (qua_global[2] < qua_global[1] &&
         qua_global[2] < qua_global[3])) { /* 0-2 diag split of quad */
                                           /* 2-1
                                              |\|
                                              3-0 */
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[2];
      tri_nodes[2] = qua_nodes[3];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[2];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
    } else { /* 3-1 diag split of quad */
             /* 2-1
                |/|
                3-0 */
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[3];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
      tri_nodes[0] = qua_nodes[2];
      tri_nodes[1] = qua_nodes[3];
      tri_nodes[2] = qua_nodes[1];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
    }
  }

  ref_free(mark_copy);
  ref_free(mark);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_adj_empty(ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_hex(ref_grid)), node)) {
      RSS(ref_node_remove_without_global(ref_node, node), "hang node");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_extract_tri(REF_GRID ref_grid, REF_CELL *ref_cell_ptr) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri;
  REF_CELL qua;
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua_nodes[REF_CELL_MAX_SIZE_PER];
  REF_GLOB qua_global[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, node;
  RSS(ref_cell_deep_copy(ref_cell_ptr, ref_grid_tri(ref_grid)),
      "deep tri copy");
  tri = *ref_cell_ptr;
  qua = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes(qua, cell, qua_nodes) {
    tri_nodes[3] = qua_nodes[4]; /* patch id */

    for (node = 0; node < ref_cell_node_per(qua); node++)
      qua_global[node] = ref_node_global(ref_node, qua_nodes[node]);

    if ((qua_global[0] < qua_global[1] && qua_global[0] < qua_global[3]) ||
        (qua_global[2] < qua_global[1] &&
         qua_global[2] < qua_global[3])) { /* 0-2 diag split of quad */
                                           /* 2-1
                                              |\|
                                              3-0 */
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[2];
      tri_nodes[2] = qua_nodes[3];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[2];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
    } else { /* 3-1 diag split of quad */
             /* 2-1
                |/|
                3-0 */
      tri_nodes[0] = qua_nodes[0];
      tri_nodes[1] = qua_nodes[1];
      tri_nodes[2] = qua_nodes[3];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
      tri_nodes[0] = qua_nodes[2];
      tri_nodes[1] = qua_nodes[3];
      tri_nodes[2] = qua_nodes[1];
      RSS(ref_shard_cell_add_local(ref_node, tri, tri_nodes), "a tri");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_extract_tet(REF_GRID ref_grid, REF_CELL *ref_cell_ptr) {
  RSS(ref_cell_deep_copy(ref_cell_ptr, ref_grid_tet(ref_grid)),
      "deep tri copy");
  return REF_SUCCESS;
}
