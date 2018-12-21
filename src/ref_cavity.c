
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

#include "ref_cavity.h"

#include "ref_list.h"
#include "ref_malloc.h"

#include "ref_twod.h"

#include "ref_adapt.h"
#include "ref_edge.h"
#include "ref_sort.h"

#include "ref_dict.h"

REF_STATUS ref_cavity_create(REF_CAVITY *ref_cavity_ptr) {
  REF_CAVITY ref_cavity;
  REF_INT face;

  ref_malloc(*ref_cavity_ptr, 1, REF_CAVITY_STRUCT);
  ref_cavity = (*ref_cavity_ptr);

  ref_cavity_state(ref_cavity) = REF_CAVITY_UNKNOWN;

  ref_cavity_grid(ref_cavity) = (REF_GRID)NULL;
  ref_cavity_node(ref_cavity) = REF_EMPTY;
  ref_cavity_n(ref_cavity) = 0;

  ref_cavity_max(ref_cavity) = 10;

  ref_malloc_init(ref_cavity->f2n, ref_cavity_max(ref_cavity) * 3, REF_INT, 0);
  for (face = 0; face < ref_cavity_max(ref_cavity); face++) {
    ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
    ref_cavity_f2n(ref_cavity, 1, face) = face + 1;
  }
  ref_cavity_f2n(ref_cavity, 1, ref_cavity_max(ref_cavity) - 1) = REF_EMPTY;
  ref_cavity_blank(ref_cavity) = 0;

  RSS(ref_list_create(&(ref_cavity->ref_list)), "add list");

  ref_cavity->debug = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_free(REF_CAVITY ref_cavity) {
  if (NULL == (void *)ref_cavity) return REF_NULL;
  ref_list_free(ref_cavity->ref_list);
  ref_free(ref_cavity->f2n);
  ref_free(ref_cavity);
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_inspect(REF_CAVITY ref_cavity) {
  REF_INT face, node;
  if (NULL == (void *)ref_cavity) return REF_NULL;
  printf("n = %d max = %d blank = %d\n", ref_cavity_n(ref_cavity),
         ref_cavity_max(ref_cavity), ref_cavity_blank(ref_cavity));
  for (face = 0; face < ref_cavity_max(ref_cavity); face++) {
    printf(" f2n[%d] = ", face);
    for (node = 0; node < 3; node++)
      printf(" %d ", ref_cavity_f2n(ref_cavity, node, face));
    printf("\n");
  }
  RSS(ref_list_inspect(ref_cavity_list(ref_cavity)), "insp");
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_insert(REF_CAVITY ref_cavity, REF_INT *nodes) {
  REF_INT node, face;
  REF_INT orig, chunk;
  REF_BOOL reversed;

  RXS(ref_cavity_find(ref_cavity, nodes, &face, &reversed), REF_NOT_FOUND,
      "find existing");

  if (REF_EMPTY != face) {
    if (reversed) { /* two faces with opposite orientation destroy each other */
      ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity, 1, face) = ref_cavity_blank(ref_cavity);
      ref_cavity_blank(ref_cavity) = face;
      ref_cavity_n(ref_cavity)--;
      return REF_SUCCESS;
    } else { /* can't happen, added same face twice */
      return REF_INVALID;
    }
  }

  /* if I need to grow my array of faces */
  if (REF_EMPTY == ref_cavity_blank(ref_cavity)) {
    orig = ref_cavity_max(ref_cavity);
    chunk = MAX(100, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_cavity_max(ref_cavity) = orig + chunk;

    ref_realloc(ref_cavity->f2n, 3 * ref_cavity_max(ref_cavity), REF_INT);

    for (face = orig; face < ref_cavity_max(ref_cavity); face++) {
      ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity, 1, face) = face + 1;
    }
    ref_cavity_f2n(ref_cavity, 1, (ref_cavity->max) - 1) = REF_EMPTY;
    ref_cavity_blank(ref_cavity) = orig;
  }

  face = ref_cavity_blank(ref_cavity);
  ref_cavity_blank(ref_cavity) = ref_cavity_f2n(ref_cavity, 1, face);
  for (node = 0; node < 3; node++)
    ref_cavity_f2n(ref_cavity, node, face) = nodes[node];

  ref_cavity_n(ref_cavity)++;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_find(REF_CAVITY ref_cavity, REF_INT *nodes,
                           REF_INT *found_face, REF_BOOL *reversed) {
  REF_INT face;

  *found_face = REF_EMPTY;

  each_ref_cavity_valid_face(ref_cavity, face) {
    if ((nodes[0] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[1] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[2] == ref_cavity_f2n(ref_cavity, 2, face)) ||
        (nodes[1] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[2] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[0] == ref_cavity_f2n(ref_cavity, 2, face)) ||
        (nodes[2] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[0] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[1] == ref_cavity_f2n(ref_cavity, 2, face))) {
      *found_face = face;
      *reversed = REF_FALSE;
      return REF_SUCCESS;
    }
    if ((nodes[2] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[1] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[0] == ref_cavity_f2n(ref_cavity, 2, face)) ||
        (nodes[1] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[0] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[2] == ref_cavity_f2n(ref_cavity, 2, face)) ||
        (nodes[0] == ref_cavity_f2n(ref_cavity, 0, face) &&
         nodes[2] == ref_cavity_f2n(ref_cavity, 1, face) &&
         nodes[1] == ref_cavity_f2n(ref_cavity, 2, face))) {
      *found_face = face;
      *reversed = REF_TRUE;
      return REF_SUCCESS;
    }
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_cavity_add_tet(REF_CAVITY ref_cavity, REF_INT tet) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT cell_face, node;
  REF_INT face_nodes[3];
  REF_INT already_have_it;

  RAS(ref_cell_valid(ref_cell, tet), "invalid tet");

  RSS(ref_list_contains(ref_cavity_list(ref_cavity), tet, &already_have_it),
      "have tet?");
  if (already_have_it) return REF_SUCCESS;

  RSS(ref_list_push(ref_cavity_list(ref_cavity), tet), "save tet");

  each_ref_cell_cell_face(ref_cell, cell_face) {
    each_ref_cavity_face_node(ref_cavity, node) {
      face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, tet);
      if (!ref_node_owned(ref_node, face_nodes[node])) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      }
    }
    RSS(ref_cavity_insert(ref_cavity, face_nodes), "tet side");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_rm_tet(REF_CAVITY ref_cavity, REF_INT tet) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_INT cell_face;
  REF_INT face_nodes[4];

  RSS(ref_list_delete(ref_cavity_list(ref_cavity), tet), "dump tet");

  each_ref_cell_cell_face(ref_cell, cell_face) {
    /* reverse face nodes orientation */
    face_nodes[0] = ref_cell_f2n(ref_cell, 1, cell_face, tet);
    face_nodes[1] = ref_cell_f2n(ref_cell, 0, cell_face, tet);
    face_nodes[2] = ref_cell_f2n(ref_cell, 2, cell_face, tet);
    RSS(ref_cavity_insert(ref_cavity, face_nodes), "tet side");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace_tet(REF_CAVITY ref_cavity) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face;
  REF_INT cell;
  REF_INT i;
  REF_DBL volume;

  if (ref_cavity_debug(ref_cavity))
    RSS(ref_cavity_tec(ref_cavity, "ref_cavity_tet_replace.tec"),
        "tec for enlarge_face fail");

  each_ref_cavity_valid_face(ref_cavity, face) {
    nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
    nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
    nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
    nodes[3] = ref_cavity_node(ref_cavity);
    if (nodes[3] == nodes[0] || nodes[3] == nodes[1] || nodes[3] == nodes[2])
      continue; /* attached face */
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add");
    RSS(ref_node_tet_vol(ref_node, nodes, &volume), "norm");
    if (volume <= ref_node_min_volume(ref_node))
      printf("%d %d %d %d %e\n", nodes[0], nodes[1], nodes[2], nodes[3],
             volume);
  }

  while (ref_list_n(ref_cavity_list(ref_cavity)) > 0) {
    RSS(ref_list_pop(ref_cavity_list(ref_cavity), &cell), "list");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "rm");
    RSS(ref_cell_remove(ref_cell, cell), "rm");
    for (i = 0; i < 4; i++)
      if (ref_adj_empty(ref_cell_adj(ref_cell), nodes[i]))
        RSS(ref_node_remove(ref_node, nodes[i]), "remove");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_empty(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                 REF_INT node) {
  ref_cavity_grid(ref_cavity) = ref_grid;
  ref_cavity_node(ref_cavity) = node;
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_ball(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                REF_INT node) {
  REF_INT item, cell;
  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "init form empty");

  each_ref_cell_having_node(ref_grid_tet(ref_grid), node, item, cell) {
    RSS(ref_cavity_add_tet(ref_cavity, cell), "insert");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_gem(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                               REF_INT node0, REF_INT node1, REF_INT node) {
  REF_INT cell, ncell;
  REF_INT cell_to_add[50];
  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "init form empty");

  RSS(ref_cell_list_with2(ref_grid_tet(ref_grid), node0, node1, 50, &ncell,
                          cell_to_add),
      "get list");
  for (cell = 0; cell < ncell; cell++) {
    RSS(ref_cavity_add_tet(ref_cavity, cell_to_add[cell]), "insert");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_edge_split(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                      REF_INT node0, REF_INT node1,
                                      REF_INT new_node) {
  REF_INT cell, ncell;
  REF_INT cell_to_add[50];
  REF_INT face;
  REF_INT node, nodes[3];
  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, new_node), "init form empty");

  RSS(ref_cell_list_with2(ref_grid_tet(ref_grid), node0, node1, 50, &ncell,
                          cell_to_add),
      "get list");
  for (cell = 0; cell < ncell; cell++) {
    RSS(ref_cavity_add_tet(ref_cavity, cell_to_add[cell]), "insert");
  }

  each_ref_cavity_valid_face(ref_cavity, face) {
    each_ref_cavity_face_node(ref_cavity, node) nodes[node] =
        ref_cavity_f2n(ref_cavity, node, face);
    if ((node0 == nodes[1] && node1 == nodes[2]) ||
        (node1 == nodes[1] && node0 == nodes[2])) {
      ref_cavity_f2n(ref_cavity, 1, face) = new_node;
      nodes[2] = new_node;
      RSS(ref_cavity_insert(ref_cavity, nodes), "insert edge 0");
      continue;
    }
    if ((node0 == nodes[0] && node1 == nodes[2]) ||
        (node1 == nodes[0] && node0 == nodes[2])) {
      ref_cavity_f2n(ref_cavity, 0, face) = new_node;
      nodes[2] = new_node;
      RSS(ref_cavity_insert(ref_cavity, nodes), "insert edge 1");
      continue;
    }
    if ((node0 == nodes[0] && node1 == nodes[1]) ||
        (node1 == nodes[0] && node0 == nodes[1])) {
      ref_cavity_f2n(ref_cavity, 0, face) = new_node;
      nodes[1] = new_node;
      RSS(ref_cavity_insert(ref_cavity, nodes), "insert edge 2");
      continue;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_visible(REF_CAVITY ref_cavity, REF_INT face,
                              REF_BOOL *visible) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL volume;

  *visible = REF_FALSE;

  nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
  nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
  nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
  nodes[3] = ref_cavity_node(ref_cavity);

  RSS(ref_node_tet_vol(ref_node, nodes, &volume), "norm");

  if (volume <= ref_node_min_volume(ref_node)) return REF_SUCCESS;

  *visible = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_visible(REF_CAVITY ref_cavity) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_INT face;
  REF_BOOL local;
  REF_BOOL visible;
  REF_BOOL keep_growing;

  RAS(ref_node_owned(ref_node, node), "cavity part must own node");

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge start %d tets %d faces\n",
           ref_list_n(ref_cavity_list(ref_cavity)), ref_cavity_n(ref_cavity));

  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  /* make sure all cell nodes to be replaced are owned */
  RSS(ref_cavity_local(ref_cavity, &local), "locality");
  if (!local) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  keep_growing = REF_TRUE;
  while (keep_growing) {
    keep_growing = REF_FALSE;
    each_ref_cavity_valid_face(ref_cavity, face) {
      /* skip a face attached to node */
      if (node == ref_cavity_f2n(ref_cavity, 0, face) ||
          node == ref_cavity_f2n(ref_cavity, 1, face))
        continue;
      if (node == ref_cavity_f2n(ref_cavity, 2, face)) continue;

      RSS(ref_cavity_visible(ref_cavity, face, &visible), "free");
      if (!visible) {
        RSS(ref_cavity_enlarge_face(ref_cavity, face), "enlarge face");
        if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) {
          if (ref_cavity_debug(ref_cavity))
            RSS(ref_cavity_tec(ref_cavity, "enlarge.tec"),
                "tec for enlarge_face fail");
          return REF_SUCCESS;
        }
        keep_growing = REF_TRUE;
      }
    }
  }

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge final %d tets %d faces\n",
           ref_list_n(ref_cavity_list(ref_cavity)), ref_cavity_n(ref_cavity));

  if (ref_cavity_debug(ref_cavity)) RSS(ref_cavity_topo(ref_cavity), "topo");

  ref_cavity_state(ref_cavity) = REF_CAVITY_VISIBLE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_shrink_visible(REF_CAVITY ref_cavity) {
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_INT face;
  REF_BOOL visible;
  REF_BOOL keep_growing;
  REF_STATUS status;

  keep_growing = REF_TRUE;
  while (keep_growing) {
    keep_growing = REF_FALSE;
    each_ref_cavity_valid_face(ref_cavity, face) {
      /* skip a face attached to node */
      if (node == ref_cavity_f2n(ref_cavity, 0, face) ||
          node == ref_cavity_f2n(ref_cavity, 1, face))
        continue;
      if (node == ref_cavity_f2n(ref_cavity, 2, face)) continue;

      RSS(ref_cavity_visible(ref_cavity, face, &visible), "free");
      if (!visible) {
        status = ref_cavity_shrink_face(ref_cavity, face);
        RXS(status, REF_INVALID, "shrink face");
        if (REF_SUCCESS == status) {
          keep_growing = REF_TRUE;
        } else {
          RSS(ref_cavity_tec(ref_cavity, "ref_cavity_debug_shrink.tec"), "tec");
          THROW("boundary, see debug");
        }
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_face(REF_CAVITY ref_cavity, REF_INT face) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT face_node, face_nodes[4];
  REF_BOOL have_cell0, have_cell1;
  REF_INT tet0, tet1;

  /* make sure all face nodes are owned */
  each_ref_cavity_face_node(ref_cavity, face_node) {
    if (!ref_node_owned(ref_node,
                        ref_cavity_f2n(ref_cavity, face_node, face))) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }
  }

  face_nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
  face_nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
  face_nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
  face_nodes[3] = face_nodes[0];
  RSB(ref_cell_with_face(ref_grid_tet(ref_grid), face_nodes, &tet0, &tet1),
      "found too many tets with face_nodes", {
        printf("%d face_nodes %d %d %d %d\n", face, face_nodes[0],
               face_nodes[1], face_nodes[2], face_nodes[3]);
        ref_cavity_tec(ref_cavity, "ref_cavity_error_too_many_tet.tec");
      });
  if (REF_EMPTY == tet0) THROW("cavity tets missing");
  if (REF_EMPTY == tet1) {
    REF_INT tri;
    RSS(ref_cell_with(ref_grid_tri(ref_grid), face_nodes, &tri),
        "verify boundary face");
    ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_contains(ref_cavity_list(ref_cavity), tet0, &have_cell0),
      "cell0");
  RSS(ref_list_contains(ref_cavity_list(ref_cavity), tet1, &have_cell1),
      "cell1");
  if (have_cell0 == have_cell1) THROW("cavity same state");
  if (have_cell0) RSS(ref_cavity_add_tet(ref_cavity, tet1), "add c1");
  if (have_cell1) RSS(ref_cavity_add_tet(ref_cavity, tet0), "add c0");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_shrink_face(REF_CAVITY ref_cavity, REF_INT face) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_INT face_nodes[4];
  REF_BOOL have_cell0, have_cell1;
  REF_INT tet0, tet1;

  face_nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
  face_nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
  face_nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
  face_nodes[3] = face_nodes[0];
  RSS(ref_cell_with_face(ref_cell, face_nodes, &tet0, &tet1),
      "unable to find tets with face");
  if (REF_EMPTY == tet0) THROW("cavity tets missing");
  /* boundary is allowed, use the interior tet */

  RSS(ref_list_contains(ref_cavity_list(ref_cavity), tet0, &have_cell0),
      "cell0");
  RSS(ref_list_contains(ref_cavity_list(ref_cavity), tet1, &have_cell1),
      "cell1");
  if (have_cell0 == have_cell1) THROW("cavity same state");
  if (!have_cell0) RSS(ref_cavity_rm_tet(ref_cavity, tet1), "add c1");
  if (!have_cell1) RSS(ref_cavity_rm_tet(ref_cavity, tet0), "add c0");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_tec(REF_CAVITY ref_cavity, const char *filename) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_DICT node_dict, face_dict;
  REF_CELL ref_cell;
  REF_INT face, face_node;
  REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local;
  REF_DBL xyz_phys[3];
  char *zonetype;
  FILE *f;

  ref_cell = ref_grid_tet(ref_grid);
  zonetype = "fetetrahedron";

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine cavity\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  each_ref_list_item(ref_cavity_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_list(ref_cavity), item);
    RSS(ref_dict_store(face_dict, cell, 0), "store");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  fprintf(
      f, "zone t=\"old\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      ref_dict_n(node_dict), ref_dict_n(face_dict), "point", zonetype);
  for (item = 0; item < ref_dict_n(node_dict); item++) {
    local = ref_dict_key(node_dict, item);
    xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local);
    xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local);
    xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local);
    fprintf(f, " %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
  }

  for (item = 0; item < ref_dict_n(face_dict); item++) {
    cell = ref_dict_key(face_dict, item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_location(node_dict, nodes[cell_node], &local), "ret");
      fprintf(f, " %d", local + 1);
    }
    fprintf(f, "\n");
  }

  RSS(ref_dict_free(face_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  RSS(ref_dict_store(node_dict, node, 0), "store");
  each_ref_cavity_valid_face(ref_cavity, face) {
    RSS(ref_dict_store(face_dict, face, 0), "store");
    each_ref_cavity_face_node(ref_cavity, face_node)
        RSS(ref_dict_store(node_dict,
                           ref_cavity_f2n(ref_cavity, face_node, face), 0),
            "store");
  }

  fprintf(
      f, "zone t=\"new\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      ref_dict_n(node_dict), ref_dict_n(face_dict), "point", zonetype);
  for (item = 0; item < ref_dict_n(node_dict); item++) {
    local = ref_dict_key(node_dict, item);
    xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local);
    xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local);
    xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local);
    fprintf(f, " %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
  }

  for (item = 0; item < ref_dict_n(face_dict); item++) {
    face = ref_dict_key(face_dict, item);
    RSS(ref_dict_location(node_dict, node, &local), "center node");
    fprintf(f, " %d", local + 1);
    each_ref_cavity_face_node(ref_cavity, face_node) {
      RSS(ref_dict_location(
              node_dict, ref_cavity_f2n(ref_cavity, face_node, face), &local),
          "ret");
      fprintf(f, " %d", local + 1);
    }
    fprintf(f, "\n");
  }

  RSS(ref_dict_free(face_dict), "free face");
  RSS(ref_dict_free(node_dict), "free nodes");

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_local(REF_CAVITY ref_cavity, REF_BOOL *local) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face, face_node;

  *local = REF_FALSE;

  each_ref_list_item(ref_cavity_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    each_ref_cell_cell_node(
        ref_cell, cell_node) if (!ref_node_owned(ref_node, nodes[cell_node])) {
      *local = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  each_ref_cavity_valid_face(ref_cavity, face) {
    each_ref_cavity_face_node(
        ref_cavity,
        face_node) if (!ref_node_owned(ref_node,
                                       ref_cavity_f2n(ref_cavity, face_node,
                                                      face))) {
      *local = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  *local = REF_TRUE;
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_change(REF_CAVITY ref_cavity, REF_DBL *min_del,
                             REF_DBL *min_add) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, min_quality, total_quality;
  REF_INT face, face_node, n, n_del, n_add;
  REF_BOOL skip;

  *min_del = -2.0;
  *min_add = -2.0;

  n = 0;
  min_quality = 1.0;
  total_quality = 0.0;
  each_ref_list_item(ref_cavity_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "new qual");
    n++;
    total_quality += quality;
    min_quality = MIN(min_quality, quality);
  }
  if (REF_FALSE && n > 0)
    printf("- min %12.8f avg %12.8f n %d\n", min_quality,
           total_quality / ((REF_DBL)n), n);
  *min_del = min_quality;
  n_del = n;

  n = 0;
  min_quality = 1.0;
  total_quality = 0.0;
  each_ref_cavity_valid_face(ref_cavity, face) {
    skip = REF_FALSE;
    /* skip a collapsed triangle that in on the boundary of cavity */
    each_ref_cavity_face_node(
        ref_cavity,
        face_node) if (node == ref_cavity_f2n(ref_cavity, face_node, face))
        skip = REF_TRUE;
    if (skip) continue;
    each_ref_cavity_face_node(ref_cavity, face_node) nodes[face_node] =
        ref_cavity_f2n(ref_cavity, face_node, face);
    nodes[3] = node;
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "new qual");
    n++;
    total_quality += quality;
    min_quality = MIN(min_quality, quality);
  }
  if (REF_FALSE && n > 0)
    printf("+ min %12.8f avg %12.8f n %d\n", min_quality,
           total_quality / ((REF_DBL)n), n);
  *min_add = min_quality;
  n_add = n;

  printf(" min %12.8f <- %12.8f diff %12.8f n %d <- %d\n", *min_add, *min_del,
         *min_add - *min_del, n_add, n_del);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_topo(REF_CAVITY ref_cavity) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_INT item, cell, face, face_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  each_ref_list_item(ref_cavity_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    printf("old %d %d %d %d\n", nodes[0], nodes[1], nodes[2], nodes[3]);
  }

  each_ref_cavity_valid_face(ref_cavity, face) {
    if (node == ref_cavity_f2n(ref_cavity, 0, face)) continue;
    if (node == ref_cavity_f2n(ref_cavity, 1, face)) continue;
    if (node == ref_cavity_f2n(ref_cavity, 2, face)) continue;
    printf("new ");
    for (face_node = 0; face_node < 3; face_node++)
      printf(" %d ", ref_cavity_f2n(ref_cavity, face_node, face));
    printf(" %d ", node);
    printf("\n");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, min_del, min_add;
  REF_CAVITY ref_cavity;
  REF_INT other, cell_edge;

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    if (quality < 0.1) {
      printf("cell %d qual %f\n", cell, quality);
      each_ref_cell_cell_edge(ref_cell, cell_edge) {
        RSS(ref_cavity_create(&ref_cavity), "create");
        other = 0;
        if (other == ref_cell_e2n_gen(ref_cell, 0, cell_edge)) other++;
        if (other == ref_cell_e2n_gen(ref_cell, 1, cell_edge)) other++;
        RSS(ref_cavity_form_gem(ref_cavity, ref_grid,
                                nodes[ref_cell_e2n_gen(ref_cell, 0, cell_edge)],
                                nodes[ref_cell_e2n_gen(ref_cell, 1, cell_edge)],
                                nodes[other]),
            "cavity gem");
        RSS(ref_cavity_enlarge_visible(ref_cavity), "enlarge viz");
        RSS(ref_cavity_change(ref_cavity, &min_del, &min_add), "change");
        if (min_add - min_del > 0.01 && min_add > 0.1) {
          printf("cavity accepted\n");
          RSS(ref_cavity_replace_tet(ref_cavity), "replace");
          RSS(ref_cavity_free(ref_cavity), "free");
          break;
        } else {
          RSS(ref_cavity_free(ref_cavity), "free");
        }
      }
    }
  }

  return REF_SUCCESS;
}
