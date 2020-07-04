
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

#include "ref_cavity.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_sort.h"
#include "ref_swap.h"

REF_STATUS ref_cavity_create(REF_CAVITY *ref_cavity_ptr) {
  REF_CAVITY ref_cavity;
  REF_INT seg, face;

  ref_malloc(*ref_cavity_ptr, 1, REF_CAVITY_STRUCT);
  ref_cavity = (*ref_cavity_ptr);

  ref_cavity_state(ref_cavity) = REF_CAVITY_UNKNOWN;

  ref_cavity_grid(ref_cavity) = (REF_GRID)NULL;
  ref_cavity_node(ref_cavity) = REF_EMPTY;
  ref_cavity_surf_node(ref_cavity) = REF_EMPTY;

  ref_cavity_nseg(ref_cavity) = 0;
  ref_cavity_maxseg(ref_cavity) = 10;

  ref_malloc_init(ref_cavity->s2n, ref_cavity_maxseg(ref_cavity) * 3, REF_INT,
                  0);
  for (seg = 0; seg < ref_cavity_maxseg(ref_cavity); seg++) {
    ref_cavity_s2n(ref_cavity, 0, seg) = REF_EMPTY;
    ref_cavity_s2n(ref_cavity, 1, seg) = seg + 1;
  }
  ref_cavity_s2n(ref_cavity, 1, ref_cavity_maxseg(ref_cavity) - 1) = REF_EMPTY;
  ref_cavity_blankseg(ref_cavity) = 0;

  ref_cavity_nface(ref_cavity) = 0;
  ref_cavity_maxface(ref_cavity) = 10;

  ref_malloc_init(ref_cavity->f2n, ref_cavity_maxface(ref_cavity) * 3, REF_INT,
                  0);
  for (face = 0; face < ref_cavity_maxface(ref_cavity); face++) {
    ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
    ref_cavity_f2n(ref_cavity, 1, face) = face + 1;
  }
  ref_cavity_f2n(ref_cavity, 1, ref_cavity_maxface(ref_cavity) - 1) = REF_EMPTY;
  ref_cavity_blankface(ref_cavity) = 0;

  RSS(ref_list_create(&(ref_cavity->tri_list)), "tri list");
  RSS(ref_list_create(&(ref_cavity->tet_list)), "tet list");

  ref_cavity->min_normdev = 0.5;

  ref_cavity->split_node0 = REF_EMPTY;
  ref_cavity->split_node1 = REF_EMPTY;
  ref_cavity->collapse_node0 = REF_EMPTY;
  ref_cavity->collapse_node1 = REF_EMPTY;

  ref_cavity->debug = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_free(REF_CAVITY ref_cavity) {
  if (NULL == (void *)ref_cavity) return REF_NULL;
  ref_list_free(ref_cavity->tet_list);
  ref_list_free(ref_cavity->tri_list);
  ref_free(ref_cavity->f2n);
  ref_free(ref_cavity->s2n);
  ref_free(ref_cavity);
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_inspect(REF_CAVITY ref_cavity) {
  REF_INT face, node;
  REF_INT item, cell, i, nodes[REF_CELL_MAX_SIZE_PER];

  if (NULL == (void *)ref_cavity) return REF_NULL;

  printf("node %d\nnface = %d maxface = %d blankface = %d\n",
         ref_cavity_node(ref_cavity), ref_cavity_nface(ref_cavity),
         ref_cavity_maxface(ref_cavity), ref_cavity_blankface(ref_cavity));
  for (face = 0; face < ref_cavity_maxface(ref_cavity); face++) {
    printf(" f2n[%d] = ", face);
    for (node = 0; node < 3; node++)
      printf(" %d ", ref_cavity_f2n(ref_cavity, node, face));
    printf("\n");
  }
  RSS(ref_list_inspect(ref_cavity_tet_list(ref_cavity)), "insp");
  printf("seg node %d\nnseg = %d maxseg = %d blankseg = %d\n",
         ref_cavity_seg_node(ref_cavity), ref_cavity_nseg(ref_cavity),
         ref_cavity_maxseg(ref_cavity), ref_cavity_blankseg(ref_cavity));
  for (face = 0; face < ref_cavity_maxseg(ref_cavity); face++) {
    printf(" s2n[%d] = ", face);
    for (node = 0; node < 3; node++)
      printf(" %d ", ref_cavity_s2n(ref_cavity, node, face));
    printf("\n");
  }
  RSS(ref_list_inspect(ref_cavity_tri_list(ref_cavity)), "insp");
  each_ref_list_item(ref_cavity_tri_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tri_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_grid_tri(ref_cavity_grid(ref_cavity)), cell, nodes),
        "cell");
    for (i = 0; i < 4; i++) printf(" %d", nodes[i]);
    printf("\n");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_remove_seg_face(REF_CAVITY ref_cavity,
                                      REF_INT *seg_nodes) {
  REF_INT face, face_nodes[3];
  REF_BOOL reversed;

  if (ref_list_n(ref_cavity_tet_list(ref_cavity)) == 0) return REF_SUCCESS;
  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  face_nodes[0] = seg_nodes[0];
  face_nodes[1] = seg_nodes[1];
  face_nodes[2] = ref_cavity_seg_node(ref_cavity);

  RSS(ref_cavity_find_face(ref_cavity, face_nodes, &face, &reversed),
      "find existing face for segment");
  RUS(REF_EMPTY, face, "face for a segment missing");
  RAS(reversed, "face for a segment removal should be backwards");

  ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
  ref_cavity_f2n(ref_cavity, 1, face) = ref_cavity_blankface(ref_cavity);
  ref_cavity_blankface(ref_cavity) = face;
  ref_cavity_nface(ref_cavity)--;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_seg_face(REF_CAVITY ref_cavity, REF_INT *seg_nodes) {
  REF_INT face_nodes[3];

  if (ref_list_n(ref_cavity_tet_list(ref_cavity)) == 0) return REF_SUCCESS;
  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;
  if (seg_nodes[0] == ref_cavity_seg_node(ref_cavity) ||
      seg_nodes[1] == ref_cavity_seg_node(ref_cavity))
    return REF_SUCCESS;

  face_nodes[0] = seg_nodes[0];
  face_nodes[1] = seg_nodes[1];
  face_nodes[2] = ref_cavity_seg_node(ref_cavity);

  RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "add face for seg");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_remove_seg_add_tets(REF_CAVITY ref_cavity,
                                          REF_INT *seg_nodes) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT ntri, tri_list[2];
  REF_INT item, cell_node, cell;
  REF_BOOL already_have_it, all_local;
  REF_INT cell_face, face_node, face_nodes[3], tri_cell;

  if (ref_list_n(ref_cavity_tet_list(ref_cavity)) == 0) return REF_SUCCESS;
  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  RSS(ref_cell_list_with2(ref_grid_tri(ref_grid), seg_nodes[0], seg_nodes[1], 2,
                          &ntri, tri_list),
      "tri with2");
  REIS(2, ntri, "cavity segment does not have two tri");

  each_ref_cell_having_node2(ref_cell, seg_nodes[0], seg_nodes[1], item,
                             cell_node, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    if (already_have_it) continue;
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local), "local cell");
    if (!all_local) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }
    each_ref_cell_cell_face(ref_cell, cell_face) {
      each_ref_cavity_face_node(ref_cavity, face_node) {
        face_nodes[face_node] =
            ref_cell_f2n(ref_cell, face_node, cell_face, cell);
      }
      RXS(ref_cell_with(ref_grid_tri(ref_grid), face_nodes, &tri_cell),
          REF_NOT_FOUND, "search for boundary tri");
      if (REF_EMPTY != tri_cell &&
          (tri_cell == tri_list[0] || tri_cell == tri_list[1])) {
        continue;
      }
      RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet face rm seg");
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_insert_seg(REF_CAVITY ref_cavity, REF_INT *nodes) {
  REF_INT node, seg;
  REF_INT orig, chunk;
  REF_BOOL reversed;

  RXS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed), REF_NOT_FOUND,
      "find existing");

  if (REF_EMPTY != seg) {
    if (reversed) { /* two segs with opposite orientation destroy each other */
      if (nodes[2] != ref_cavity_s2n(ref_cavity, 2, seg)) {
        /* mixing faceid would violate topology */
        ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
        return REF_SUCCESS;
      }
      ref_cavity_s2n(ref_cavity, 0, seg) = REF_EMPTY;
      ref_cavity_s2n(ref_cavity, 1, seg) = ref_cavity_blankseg(ref_cavity);
      ref_cavity_blankseg(ref_cavity) = seg;

      ref_cavity_nseg(ref_cavity)--;
      RSS(ref_cavity_remove_seg_face(ref_cavity, nodes), "rm seg face");
      RSS(ref_cavity_remove_seg_add_tets(ref_cavity, nodes), "rm seg tets");

      return REF_SUCCESS;
    } else { /* can't happen, added same seg twice */
      return REF_INVALID;
    }
  }

  /* if I need to grow my array of segs */
  if (REF_EMPTY == ref_cavity_blankseg(ref_cavity)) {
    orig = ref_cavity_maxseg(ref_cavity);
    chunk = MAX(100, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_cavity_maxseg(ref_cavity) = orig + chunk;

    ref_realloc(ref_cavity->s2n, 3 * ref_cavity_maxseg(ref_cavity), REF_INT);

    for (seg = orig; seg < ref_cavity_maxseg(ref_cavity); seg++) {
      ref_cavity_s2n(ref_cavity, 0, seg) = REF_EMPTY;
      ref_cavity_s2n(ref_cavity, 1, seg) = seg + 1;
    }
    ref_cavity_s2n(ref_cavity, 1, (ref_cavity->maxseg) - 1) = REF_EMPTY;
    ref_cavity_blankseg(ref_cavity) = orig;
  }

  seg = ref_cavity_blankseg(ref_cavity);
  ref_cavity_blankseg(ref_cavity) = ref_cavity_s2n(ref_cavity, 1, seg);
  for (node = 0; node < 2; node++)
    ref_cavity_s2n(ref_cavity, node, seg) = nodes[node];
  ref_cavity_s2n(ref_cavity, 2, seg) = nodes[2]; /* faceid */

  ref_cavity_nseg(ref_cavity)++;
  RSS(ref_cavity_add_seg_face(ref_cavity, nodes), "add seg face");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_find_seg(REF_CAVITY ref_cavity, REF_INT *nodes,
                               REF_INT *found_seg, REF_BOOL *reversed) {
  REF_INT seg;

  *found_seg = REF_EMPTY;

  each_ref_cavity_valid_seg(ref_cavity, seg) {
    if ((nodes[0] == ref_cavity_s2n(ref_cavity, 0, seg) &&
         nodes[1] == ref_cavity_s2n(ref_cavity, 1, seg))) {
      *found_seg = seg;
      *reversed = REF_FALSE;
      return REF_SUCCESS;
    }
    if ((nodes[1] == ref_cavity_s2n(ref_cavity, 0, seg) &&
         nodes[0] == ref_cavity_s2n(ref_cavity, 1, seg))) {
      *found_seg = seg;
      *reversed = REF_TRUE;
      return REF_SUCCESS;
    }
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_cavity_insert_face(REF_CAVITY ref_cavity, REF_INT *nodes) {
  REF_INT node, face;
  REF_INT orig, chunk;
  REF_BOOL reversed;

  RXS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed), REF_NOT_FOUND,
      "find existing");

  if (REF_EMPTY != face) {
    if (reversed) { /* two faces with opposite orientation destroy each other */
      ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity, 1, face) = ref_cavity_blankface(ref_cavity);
      ref_cavity_blankface(ref_cavity) = face;
      ref_cavity_nface(ref_cavity)--;
      return REF_SUCCESS;
    } else { /* can't happen, added same face twice */
      return REF_INVALID;
    }
  }

  /* if I need to grow my array of faces */
  if (REF_EMPTY == ref_cavity_blankface(ref_cavity)) {
    orig = ref_cavity_maxface(ref_cavity);
    chunk = MAX(100, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_cavity_maxface(ref_cavity) = orig + chunk;

    ref_realloc(ref_cavity->f2n, 3 * ref_cavity_maxface(ref_cavity), REF_INT);

    for (face = orig; face < ref_cavity_maxface(ref_cavity); face++) {
      ref_cavity_f2n(ref_cavity, 0, face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity, 1, face) = face + 1;
    }
    ref_cavity_f2n(ref_cavity, 1, (ref_cavity->maxface) - 1) = REF_EMPTY;
    ref_cavity_blankface(ref_cavity) = orig;
  }

  face = ref_cavity_blankface(ref_cavity);
  ref_cavity_blankface(ref_cavity) = ref_cavity_f2n(ref_cavity, 1, face);
  for (node = 0; node < 3; node++)
    ref_cavity_f2n(ref_cavity, node, face) = nodes[node];

  ref_cavity_nface(ref_cavity)++;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_find_face(REF_CAVITY ref_cavity, REF_INT *nodes,
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

REF_STATUS ref_cavity_find_face_with_side(REF_CAVITY ref_cavity, REF_INT node0,
                                          REF_INT node1, REF_INT *found_face) {
  REF_INT face;

  *found_face = REF_EMPTY;

  each_ref_cavity_valid_face(ref_cavity, face) {
    if ((node0 == ref_cavity_f2n(ref_cavity, 0, face) &&
         node1 == ref_cavity_f2n(ref_cavity, 1, face)) ||
        (node0 == ref_cavity_f2n(ref_cavity, 1, face) &&
         node1 == ref_cavity_f2n(ref_cavity, 2, face)) ||
        (node0 == ref_cavity_f2n(ref_cavity, 2, face) &&
         node1 == ref_cavity_f2n(ref_cavity, 0, face))) {
      if (REF_EMPTY != *found_face) { /* found face twice */
        ref_cavity_state(ref_cavity) = REF_CAVITY_INCONSISTENT;
        *found_face = REF_EMPTY;
        return REF_SUCCESS;
      }
      *found_face = face;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_verify_face_manifold(REF_CAVITY ref_cavity) {
  REF_INT face, found_face;

  if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity))
    return REF_SUCCESS;

  each_ref_cavity_valid_face(ref_cavity, face) {
    RSS(ref_cavity_find_face_with_side(
            ref_cavity, ref_cavity_f2n(ref_cavity, 1, face),
            ref_cavity_f2n(ref_cavity, 0, face), &found_face),
        "find side 01");
    if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity))
      return REF_SUCCESS;
    if (REF_EMPTY == found_face) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_INCONSISTENT;
      return REF_SUCCESS;
    }
    RSS(ref_cavity_find_face_with_side(
            ref_cavity, ref_cavity_f2n(ref_cavity, 2, face),
            ref_cavity_f2n(ref_cavity, 1, face), &found_face),
        "find side 12");
    if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity))
      return REF_SUCCESS;
    if (REF_EMPTY == found_face) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_INCONSISTENT;
      return REF_SUCCESS;
    }
    RSS(ref_cavity_find_face_with_side(
            ref_cavity, ref_cavity_f2n(ref_cavity, 0, face),
            ref_cavity_f2n(ref_cavity, 2, face), &found_face),
        "find side 20");
    if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity))
      return REF_SUCCESS;
    if (REF_EMPTY == found_face) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_INCONSISTENT;
      return REF_SUCCESS;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_verify_seg_manifold(REF_CAVITY ref_cavity) {
  REF_INT seg0, seg1, found_seg;

  if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity))
    return REF_SUCCESS;

  each_ref_cavity_valid_seg(ref_cavity, seg0) {
    found_seg = REF_EMPTY;
    each_ref_cavity_valid_seg(ref_cavity, seg1) {
      if (ref_cavity_s2n(ref_cavity, 1, seg0) ==
          ref_cavity_s2n(ref_cavity, 0, seg1)) {
        if (REF_EMPTY != found_seg) { /* found seg twice */
          ref_cavity_state(ref_cavity) = REF_CAVITY_INCONSISTENT;
          return REF_SUCCESS;
        }
        found_seg = seg1;
      }
    }
    RUS(REF_EMPTY, found_seg, "seg missing");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_tri_tet(REF_CAVITY ref_cavity, REF_INT tri) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet0, tet1;
  REF_INT face_nodes[3], cell_face, face_node;
  REF_BOOL all_local, is_tri, already_have_it;

  RSS(ref_cell_nodes(ref_grid_tri(ref_cavity_grid(ref_cavity)), tri, nodes),
      "rm");
  nodes[3] = nodes[0];
  RSS(ref_cell_with_face(ref_grid_tet(ref_cavity_grid(ref_cavity)), nodes,
                         &tet0, &tet1),
      "tet with tri");
  RUS(REF_EMPTY, tet0, "tet0 missing");
  RES(REF_EMPTY, tet1, "tet1 should not be found on boundary");
  RSS(ref_cell_all_local(ref_grid_tet(ref_cavity_grid(ref_cavity)), ref_node,
                         tet0, &all_local),
      "local cell");
  if (!all_local) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), tet0,
                        &already_have_it),
      "have tet0?");
  if (already_have_it) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), tet0), "save tet");
  each_ref_cell_cell_face(ref_grid_tet(ref_cavity_grid(ref_cavity)),
                          cell_face) {
    each_ref_cavity_face_node(ref_cavity, face_node) {
      face_nodes[face_node] =
          ref_cell_f2n(ref_grid_tet(ref_cavity_grid(ref_cavity)), face_node,
                       cell_face, tet0);
    }
    is_tri = REF_TRUE;
    is_tri = is_tri && (nodes[0] == face_nodes[0] ||
                        nodes[0] == face_nodes[1] || nodes[0] == face_nodes[2]);
    is_tri = is_tri && (nodes[1] == face_nodes[0] ||
                        nodes[1] == face_nodes[1] || nodes[1] == face_nodes[2]);
    is_tri = is_tri && (nodes[2] == face_nodes[0] ||
                        nodes[2] == face_nodes[1] || nodes[2] == face_nodes[2]);

    if (!is_tri) {
      RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet side");
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_tri(REF_CAVITY ref_cavity, REF_INT tri) {
  REF_CELL ref_cell = ref_grid_tri(ref_cavity_grid(ref_cavity));
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT cell_edge, node;
  REF_INT seg_nodes[3];
  REF_INT all_local, already_have_it;

  RAS(ref_cell_valid(ref_cell, tri), "invalid tri");

  RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), tri, &already_have_it),
      "have tri?");
  if (already_have_it) return REF_SUCCESS;

  RSS(ref_cell_all_local(ref_cell, ref_node, tri, &all_local), "local cell");
  if (!all_local) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), tri), "save tri");

  each_ref_cell_cell_edge(ref_cell, cell_edge) {
    each_ref_cavity_seg_node(ref_cavity, node) {
      seg_nodes[node] = ref_cell_e2n(ref_cell, node, cell_edge, tri);
    }
    seg_nodes[2] = ref_cell_c2n(ref_cell, ref_cell_node_per(ref_cell), tri);
    RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_tet(REF_CAVITY ref_cavity, REF_INT tet) {
  REF_CELL ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT cell_face, node;
  REF_INT face_nodes[3];
  REF_INT already_have_it;

  RAS(ref_cell_valid(ref_cell, tet), "invalid tet");

  RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), tet, &already_have_it),
      "have tet?");
  if (already_have_it) return REF_SUCCESS;

  RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), tet), "save tet");

  each_ref_cell_cell_face(ref_cell, cell_face) {
    each_ref_cavity_face_node(ref_cavity, node) {
      face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, tet);
      if (!ref_node_owned(ref_node, face_nodes[node])) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
    }
    RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet side");
    if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_check_faces(REF_CAVITY ref_cavity) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face_node, face, node, tet0, tet1, tri_cell;
  REF_INT face_nodes[4];
  REF_INT cell_face;

  node = ref_cavity_node(ref_cavity);
  each_ref_cavity_valid_face(ref_cavity, face) {
    nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
    nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
    nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
    nodes[3] = node;
    if (node == nodes[0] || node == nodes[1] || node == nodes[2])
      continue; /* attached face */
    each_ref_cell_cell_face(ref_cell, cell_face) {
      for (face_node = 0; face_node < 4; face_node++) {
        face_nodes[face_node] =
            nodes[ref_cell_f2n_gen(ref_cell, face_node, cell_face)];
      }
      RSS(ref_cell_with_face(ref_grid_tet(ref_grid), face_nodes, &tet0, &tet1),
          "found too many tets with face_nodes");
      RAS(REF_EMPTY != tet0, "no tet for cavity face");
      if (REF_EMPTY == tet1) {
        RSB(ref_cell_with(ref_grid_tri(ref_grid), face_nodes, &tri_cell),
            "no tri for missing tet1", {
              printf("%d face %d cell face %d %d %d nodes\n", face, cell_face,
                     face_nodes[0], face_nodes[1], face_nodes[2]);
              ref_cavity_inspect(ref_cavity);
            });
      }
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace(REF_CAVITY ref_cavity) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_INT face, seg;
  REF_INT cell, new_cell;
  REF_INT i, item;
  REF_DBL volume;
  REF_LIST ref_list;

  if (ref_cavity_debug(ref_cavity))
    RSS(ref_cavity_tec(ref_cavity, "ref_cavity_replace.tec"),
        "tec for replace fail");

  REIS(REF_CAVITY_VISIBLE, ref_cavity_state(ref_cavity),
       "attempt to replace cavity that is not visible");

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "replace face manifold");
  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "replace seg manifold");

  REIS(REF_CAVITY_VISIBLE, ref_cavity_state(ref_cavity),
       "attempt to replace cavity that is inconsistent");

  node = ref_cavity_node(ref_cavity);
  ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  each_ref_cavity_valid_face(ref_cavity, face) {
    nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
    nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
    nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
    nodes[3] = node;
    if (node == nodes[0] || node == nodes[1] || node == nodes[2])
      continue; /* attached face */
    for (i = 0; i < 4; i++) {
      RAS(ref_node_valid(ref_node, nodes[i]), "cavity tet nodes not valid");
      RAS(ref_node_owned(ref_node, nodes[i]), "cavity tet nodes not local");
    }
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add");
    RSS(ref_node_tet_vol(ref_node, nodes, &volume), "norm");
    if (volume <= ref_node_min_volume(ref_node))
      printf("%d %d %d %d %e\n", nodes[0], nodes[1], nodes[2], nodes[3],
             volume);
  }

  node = ref_cavity_seg_node(ref_cavity);
  ref_cell = ref_grid_tri(ref_cavity_grid(ref_cavity));
  each_ref_cavity_valid_seg(ref_cavity, seg) {
    nodes[0] = ref_cavity_s2n(ref_cavity, 0, seg);
    nodes[1] = ref_cavity_s2n(ref_cavity, 1, seg);
    nodes[2] = node;
    nodes[3] = ref_cavity_s2n(ref_cavity, 2, seg);
    if (node == nodes[0] || node == nodes[1]) continue; /* attached seg */
    for (i = 0; i < 3; i++) {
      RAS(ref_node_valid(ref_node, nodes[i]), "cavity tri nodes not valid");
      RAS(ref_node_owned(ref_node, nodes[i]), "cavity tri nodes not local");
    }
    RSS(ref_cell_add(ref_cell, nodes, &cell), "add");
    /* check validity, area? */
  }

  RSS(ref_list_create(&ref_list), "create removed node list");

  ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  each_ref_list_item(ref_cavity_tet_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tet_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "rm");
    for (i = 0; i < 4; i++) {
      RSS(ref_list_push(ref_list, nodes[i]), "tet list");
      RAS(ref_node_valid(ref_node, nodes[i]), "cavity rm tet nodes not valid");
      RAS(ref_node_owned(ref_node, nodes[i]), "cavity rm tet nodes not local");
    }
    RSS(ref_cell_remove(ref_cell, cell), "rm tet");
  }

  ref_cell = ref_grid_tri(ref_cavity_grid(ref_cavity));
  each_ref_list_item(ref_cavity_tri_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tri_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "rm");
    for (i = 0; i < 3; i++) {
      RSS(ref_list_push(ref_list, nodes[i]), "tri list");
      RAS(ref_node_valid(ref_node, nodes[i]), "cavity rm tri nodes not valid");
      RAS(ref_node_owned(ref_node, nodes[i]), "cavity rm tri nodes not local");
    }
    RSS(ref_cell_remove(ref_cell, cell), "rm tri");
  }

  each_ref_list_item(ref_list, item) {
    node = ref_list_value(ref_list, item);
    if (ref_node_valid(ref_node, node)) {
      if (ref_cell_node_empty(ref_grid_tri(ref_grid), node)) {
        RSS(ref_geom_remove_all(ref_geom, node), "remove geom");
      }
      if (ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
          ref_cell_node_empty(ref_grid_tet(ref_grid), node)) {
        RSS(ref_node_remove(ref_node, node), "remove node");
      }
    }
  }

  if (ref_cavity->split_node0 != REF_EMPTY &&
      ref_cavity->split_node1 != REF_EMPTY) {
    nodes[0] = ref_cavity->split_node0;
    nodes[1] = ref_cavity->split_node1;
    ref_cell = ref_grid_edg(ref_cavity_grid(ref_cavity));
    RXS(ref_cell_with(ref_cell, nodes, &cell), REF_NOT_FOUND, "find edg");
    if (REF_EMPTY != cell) {
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      RSS(ref_cell_remove(ref_cell, cell), "remove");
      nodes[0] = ref_cavity->split_node0;
      nodes[1] = ref_cavity_seg_node(ref_cavity);
      RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
      nodes[0] = ref_cavity_seg_node(ref_cavity);
      nodes[1] = ref_cavity->split_node1;
      RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
    }
  }

  if (ref_cavity->collapse_node0 != REF_EMPTY &&
      ref_cavity->collapse_node1 != REF_EMPTY) {
    nodes[0] = ref_cavity->collapse_node0;
    nodes[1] = ref_cavity->collapse_node1;
    ref_cell = ref_grid_edg(ref_cavity_grid(ref_cavity));
    RXS(ref_cell_with(ref_cell, nodes, &cell), REF_NOT_FOUND, "find edg");
    if (REF_EMPTY != cell) {
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      RSS(ref_cell_remove(ref_cell, cell), "remove");
      RSS(ref_cell_replace_node(ref_cell, ref_cavity->collapse_node1,
                                ref_cavity->collapse_node0),
          "replace node");
    }
  }

  /* self check to catch invalid node early */
  node = ref_cavity_node(ref_cavity);
  ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "");
    for (i = 0; i < 4; i++) {
      RAS(ref_node_valid(ref_node, nodes[i]),
          "cavity replaced tet nodes not valid");
    }
  }

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid))
    ref_cell = ref_grid_tri(ref_grid);

  each_ref_list_item(ref_list, item) {
    node = ref_list_value(ref_list, item);
    if (ref_node_valid(ref_node, node)) {
      RAS(!ref_cell_node_empty(ref_cell, node),
          "element missing for valid node");
    }
  }

  RSS(ref_list_free(ref_list), "list free");

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
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell_face, face_node, cell;
  REF_BOOL has_node;
  REF_BOOL already_have_it, all_local;
  REF_INT face_nodes[3], seg_nodes[3];
  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "init form empty");
  if (!ref_node_owned(ref_node, node)) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    RAS(!already_have_it, "added tet twice?");
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local), "local cell");
    if (!all_local) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }
    each_ref_cell_cell_face(ref_cell, cell_face) {
      each_ref_cavity_face_node(ref_cavity, face_node) {
        face_nodes[face_node] =
            ref_cell_f2n(ref_cell, face_node, cell_face, cell);
      }
      has_node = (node == face_nodes[0] || node == face_nodes[1] ||
                  node == face_nodes[2]);
      if (!has_node) {
        RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet side");
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    RAS(!already_have_it, "added tri twice?");
    RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), cell), "save tri");
    RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local), "local cell");
    if (!all_local) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }

    seg_nodes[2] = ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell);
    seg_nodes[0] = ref_cell_c2n(ref_cell, 1, cell);
    seg_nodes[1] = ref_cell_c2n(ref_cell, 2, cell);
    has_node = (node == seg_nodes[0] || node == seg_nodes[1]);
    if (!has_node) {
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    }
    seg_nodes[0] = ref_cell_c2n(ref_cell, 2, cell);
    seg_nodes[1] = ref_cell_c2n(ref_cell, 0, cell);
    has_node = (node == seg_nodes[0] || node == seg_nodes[1]);
    if (!has_node) {
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    }
    seg_nodes[0] = ref_cell_c2n(ref_cell, 0, cell);
    seg_nodes[1] = ref_cell_c2n(ref_cell, 1, cell);
    has_node = (node == seg_nodes[0] || node == seg_nodes[1]);
    if (!has_node) {
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    }
  }

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "ball face manifold");
  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "ball seg manifold");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_edge_swap(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                     REF_INT node0, REF_INT node1,
                                     REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell_node, cell;
  REF_INT node2, node3, face_node, cell_face;
  REF_BOOL has_triangle, has_node0, has_node1;
  REF_BOOL already_have_it, all_local;
  REF_INT face_nodes[3], seg_nodes[3];
  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "init form empty");
  if (!ref_node_owned(ref_node, node0) || !ref_node_owned(ref_node, node1) ||
      !ref_node_owned(ref_node, node)) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    RAS(!already_have_it, "added tet twice?");
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local), "local cell");
    if (!all_local) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }
    each_ref_cell_cell_face(ref_cell, cell_face) {
      each_ref_cavity_face_node(ref_cavity, face_node) {
        face_nodes[face_node] =
            ref_cell_f2n(ref_cell, face_node, cell_face, cell);
      }
      has_node0 = (node0 == face_nodes[0] || node0 == face_nodes[1] ||
                   node0 == face_nodes[2]);
      has_node1 = (node1 == face_nodes[0] || node1 == face_nodes[1] ||
                   node1 == face_nodes[2]);
      if (!(has_node0 && has_node1)) {
        RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet side");
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_has_side(ref_cell, node0, node1, &has_triangle),
      "triangle side");
  if (has_triangle) {
    RSS(ref_swap_node23(ref_grid, node0, node1, &node2, &node3),
        "nodes 2 and 3");

    ref_cavity_surf_node(ref_cavity) = node2;

    seg_nodes[2] = REF_EMPTY;
    each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
      RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), cell), "save tri");
      RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local),
          "local cell");
      if (!all_local) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
      seg_nodes[2] = ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell);
    }
    REIS(2, ref_list_n(ref_cavity_tri_list(ref_cavity)), "expect two tri");
    RUS(REF_EMPTY, seg_nodes[2], "faceid not set");
    seg_nodes[0] = node0;
    seg_nodes[1] = node3;
    RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    seg_nodes[0] = node3;
    seg_nodes[1] = node1;
    RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    seg_nodes[0] = node1;
    seg_nodes[1] = node2;
    RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
    seg_nodes[0] = node2;
    seg_nodes[1] = node0;
    RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
  }

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "swap face manifold");
  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "swap seg manifold");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_edge_split(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                      REF_INT node0, REF_INT node1,
                                      REF_INT new_node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell_node, cell;
  REF_INT face_node, cell_face;
  REF_BOOL has_triangle, has_node0, has_node1;
  REF_BOOL already_have_it, all_local;
  REF_INT face_nodes[3], seg_nodes[3];
  REF_INT cell_edge, n0, n1;

  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, new_node), "init form empty");
  if (!ref_node_owned(ref_node, node0) || !ref_node_owned(ref_node, node1) ||
      !ref_node_owned(ref_node, new_node)) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cavity->split_node0 = node0;
  ref_cavity->split_node1 = node1;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    RAS(!already_have_it, "added tet twice?");
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local), "local cell");
    if (!all_local) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
      return REF_SUCCESS;
    }
    each_ref_cell_cell_face(ref_cell, cell_face) {
      each_ref_cavity_face_node(ref_cavity, face_node) {
        face_nodes[face_node] =
            ref_cell_f2n(ref_cell, face_node, cell_face, cell);
      }
      has_node0 = (node0 == face_nodes[0] || node0 == face_nodes[1] ||
                   node0 == face_nodes[2]);
      has_node1 = (node1 == face_nodes[0] || node1 == face_nodes[1] ||
                   node1 == face_nodes[2]);
      if (!(has_node0 && has_node1)) {
        RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet side");
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_has_side(ref_cell, node0, node1, &has_triangle),
      "triangle side");
  if (has_triangle) {
    each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
      RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), cell), "save tri");
      RSS(ref_cell_all_local(ref_cell, ref_node, cell, &all_local),
          "local cell");
      if (!all_local) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
      each_ref_cell_cell_edge(ref_cell, cell_edge) {
        n0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
        n1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
        if (MIN(node0, node1) == MIN(n0, n1) &&
            MAX(node0, node1) == MAX(n0, n1)) {
          /* skip the edge to be split */
        } else {
          seg_nodes[0] = n0;
          seg_nodes[1] = n1;
          seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
          RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
        }
      }
    }
    RAS(1 == ref_list_n(ref_cavity_tri_list(ref_cavity)) ||
            2 == ref_list_n(ref_cavity_tri_list(ref_cavity)),
        "expect one or two tri");
    if (1 == ref_list_n(ref_cavity_tri_list(ref_cavity))) {
      each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node,
                                 cell) {
        each_ref_cell_cell_edge(ref_cell, cell_edge) {
          n0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
          n1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
          if (MIN(node0, node1) == MIN(n0, n1) &&
              MAX(node0, node1) == MAX(n0, n1)) {
            /* close topo with explicit edge split */
            seg_nodes[0] = n0;
            seg_nodes[1] = new_node;
            seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
            RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
            seg_nodes[0] = new_node;
            seg_nodes[1] = n1;
            seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
            RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri side");
          }
        }
      }
    }
  }

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "split face manifold");
  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "split seg manifold");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_form_edge_collapse(REF_CAVITY ref_cavity,
                                         REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell;

  REF_INT cell_face;
  REF_INT node, face_nodes[3], seg_nodes[3];
  REF_BOOL will_be_collapsed, already_have_it;
  REF_BOOL has_node0, has_node1;

  RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node0), "init form empty");
  if (!ref_node_owned(ref_node, node0) || !ref_node_owned(ref_node, node1)) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cavity->collapse_node0 = node0;
  ref_cavity->collapse_node1 = node1;

  ref_cell = ref_grid_tet(ref_grid);

  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    if (already_have_it) continue;
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    has_node0 = REF_FALSE;
    has_node1 = REF_FALSE;
    each_ref_cell_cell_node(ref_cell, node) {
      has_node0 = has_node0 || (node0 == ref_cell_c2n(ref_cell, node, cell));
      has_node1 = has_node1 || (node1 == ref_cell_c2n(ref_cell, node, cell));
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
    }
    will_be_collapsed = has_node0 && has_node1;
    if (will_be_collapsed) continue;
    each_ref_cell_cell_face(ref_cell, cell_face) {
      has_node0 = REF_FALSE;
      each_ref_cavity_face_node(ref_cavity, node) {
        face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
        has_node0 = has_node0 || (node0 == face_nodes[node]);
      }
      if (!has_node0)
        RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet face");
    }
  }

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    if (already_have_it) continue;
    RSS(ref_list_push(ref_cavity_tet_list(ref_cavity), cell), "save tet");
    has_node0 = REF_FALSE;
    has_node1 = REF_FALSE;
    each_ref_cell_cell_node(ref_cell, node) {
      has_node0 = has_node0 || (node0 == ref_cell_c2n(ref_cell, node, cell));
      has_node1 = has_node1 || (node1 == ref_cell_c2n(ref_cell, node, cell));
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
    }
    will_be_collapsed = has_node0 && has_node1;
    if (will_be_collapsed) continue;
    each_ref_cell_cell_face(ref_cell, cell_face) {
      has_node1 = REF_FALSE;
      each_ref_cavity_face_node(ref_cavity, node) {
        face_nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
        has_node1 = has_node1 || (node1 == face_nodes[node]);
      }
      if (!has_node1)
        RSS(ref_cavity_insert_face(ref_cavity, face_nodes), "tet face");
    }
  }

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    if (already_have_it) continue;
    RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), cell), "save tri");
    has_node0 = REF_FALSE;
    has_node1 = REF_FALSE;
    each_ref_cell_cell_node(ref_cell, node) {
      has_node0 = has_node0 || (node0 == ref_cell_c2n(ref_cell, node, cell));
      has_node1 = has_node1 || (node1 == ref_cell_c2n(ref_cell, node, cell));
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
    }
    will_be_collapsed = has_node0 && has_node1;
    if (will_be_collapsed) continue;
    if (node0 != ref_cell_c2n(ref_cell, 0, cell) &&
        node0 != ref_cell_c2n(ref_cell, 1, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 0, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 1, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 0 side 01");
    }
    if (node0 != ref_cell_c2n(ref_cell, 1, cell) &&
        node0 != ref_cell_c2n(ref_cell, 2, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 1, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 2, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 0 side 12");
    }
    if (node0 != ref_cell_c2n(ref_cell, 2, cell) &&
        node0 != ref_cell_c2n(ref_cell, 0, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 2, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 0, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 0 side 20");
    }
  }

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), cell,
                          &already_have_it),
        "have tet?");
    if (already_have_it) continue;
    RSS(ref_list_push(ref_cavity_tri_list(ref_cavity), cell), "save tri");
    has_node0 = REF_FALSE;
    has_node1 = REF_FALSE;
    each_ref_cell_cell_node(ref_cell, node) {
      has_node0 = has_node0 || (node0 == ref_cell_c2n(ref_cell, node, cell));
      has_node1 = has_node1 || (node1 == ref_cell_c2n(ref_cell, node, cell));
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        ref_cavity_state(ref_cavity) = REF_CAVITY_PARTITION_CONSTRAINED;
        return REF_SUCCESS;
      }
    }
    will_be_collapsed = has_node0 && has_node1;
    if (will_be_collapsed) continue;
    if (node1 != ref_cell_c2n(ref_cell, 0, cell) &&
        node1 != ref_cell_c2n(ref_cell, 1, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 0, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 1, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 1 side 01");
    }
    if (node1 != ref_cell_c2n(ref_cell, 1, cell) &&
        node1 != ref_cell_c2n(ref_cell, 2, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 1, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 2, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 1 side 12");
    }
    if (node1 != ref_cell_c2n(ref_cell, 2, cell) &&
        node1 != ref_cell_c2n(ref_cell, 0, cell)) {
      seg_nodes[0] = ref_cell_c2n(ref_cell, 2, cell);
      seg_nodes[1] = ref_cell_c2n(ref_cell, 0, cell);
      seg_nodes[2] = ref_cell_c2n(ref_cell, 3, cell);
      RSS(ref_cavity_insert_seg(ref_cavity, seg_nodes), "tri 1 side 20");
    }
  }

  RSB(ref_cavity_verify_face_manifold(ref_cavity), "collapse face manifold", {
    printf("at %f %f %f\n",
           ref_node_xyz(ref_node, 0, ref_cavity_node(ref_cavity)),
           ref_node_xyz(ref_node, 1, ref_cavity_node(ref_cavity)),
           ref_node_xyz(ref_node, 2, ref_cavity_node(ref_cavity)));
  });
  RSB(ref_cavity_verify_seg_manifold(ref_cavity), "collapse seg manifold", {
    printf("at %f %f %f\n",
           ref_node_xyz(ref_node, 0, ref_cavity_node(ref_cavity)),
           ref_node_xyz(ref_node, 1, ref_cavity_node(ref_cavity)),
           ref_node_xyz(ref_node, 2, ref_cavity_node(ref_cavity)));
  });

  return REF_SUCCESS;
}

static REF_STATUS ref_cavity_manifold(REF_CAVITY ref_cavity,
                                      REF_BOOL *manifold) {
  REF_INT node;
  REF_CELL ref_cell;
  REF_INT seg, face;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL contains;

  *manifold = REF_FALSE;

  node = ref_cavity_seg_node(ref_cavity);
  ref_cell = ref_grid_tri(ref_cavity_grid(ref_cavity));
  each_ref_cavity_valid_seg(ref_cavity, seg) {
    /* skip a seg attached to node */
    if (node == ref_cavity_s2n(ref_cavity, 0, seg) ||
        node == ref_cavity_s2n(ref_cavity, 1, seg))
      continue;

    nodes[0] = ref_cavity_s2n(ref_cavity, 0, seg);
    nodes[1] = ref_cavity_s2n(ref_cavity, 1, seg);
    nodes[2] = node;
    nodes[3] = ref_cavity_s2n(ref_cavity, 2, seg);

    RXS(ref_cell_with(ref_cell, nodes, &cell), REF_NOT_FOUND,
        "with manifold seach failed");
    if (REF_EMPTY != cell) {
      RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), cell, &contains),
          "contains a plan to remove");
      if (!contains) {
        *manifold = REF_FALSE;
        return REF_SUCCESS;
      }
    }
  }

  node = ref_cavity_node(ref_cavity);
  ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  each_ref_cavity_valid_face(ref_cavity, face) {
    /* skip a face attached to node */
    if (node == ref_cavity_f2n(ref_cavity, 0, face) ||
        node == ref_cavity_f2n(ref_cavity, 1, face) ||
        node == ref_cavity_f2n(ref_cavity, 2, face))
      continue;

    nodes[0] = ref_cavity_f2n(ref_cavity, 0, face);
    nodes[1] = ref_cavity_f2n(ref_cavity, 1, face);
    nodes[2] = ref_cavity_f2n(ref_cavity, 2, face);
    nodes[3] = node;

    RXS(ref_cell_with(ref_cell, nodes, &cell), REF_NOT_FOUND,
        "with manifold seach failed");
    if (REF_EMPTY != cell) {
      RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), cell, &contains),
          "contains a plan to remove");
      if (!contains) {
        *manifold = REF_FALSE;
        return REF_SUCCESS;
      }
    }
  }

  *manifold = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_combined(REF_CAVITY ref_cavity) {
  RSS(ref_cavity_enlarge_conforming(ref_cavity), "enlarge boundary");
  if (REF_CAVITY_VISIBLE != ref_cavity_state(ref_cavity)) return REF_SUCCESS;
  ref_cavity_state(ref_cavity) = REF_CAVITY_UNKNOWN;
  RSS(ref_cavity_enlarge_visible(ref_cavity), "enlarge volume");
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_conforming(REF_CAVITY ref_cavity, REF_INT seg,
                                 REF_BOOL *conforming) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normdev;
  REF_DBL sign_uv_area, uv_area;

  node = ref_cavity_seg_node(ref_cavity);

  *conforming = REF_FALSE;

  nodes[0] = ref_cavity_s2n(ref_cavity, 0, seg);
  nodes[1] = ref_cavity_s2n(ref_cavity, 1, seg);
  nodes[2] = node;
  nodes[3] = ref_cavity_s2n(ref_cavity, 2, seg);

  RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "old");

  if (normdev <= ref_cavity_min_normdev(ref_cavity)) return REF_SUCCESS;

  if (ref_geom_meshlinked(ref_geom)) {
    *conforming = REF_TRUE;
    return REF_SUCCESS;
  }

  RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
  RSS(ref_geom_uv_area_sign(ref_grid, nodes[3], &sign_uv_area), "sign");
  uv_area *= sign_uv_area;

  if (uv_area <= ref_node_min_uv_area(ref_node)) return REF_SUCCESS;

  *conforming = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_conforming(REF_CAVITY ref_cavity) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT node;
  REF_INT seg;
  REF_BOOL conforming, manifold;
  REF_BOOL keep_growing;

  node = ref_cavity_seg_node(ref_cavity);

  RAS(ref_node_owned(ref_node, node), "cavity part must own node");

  if (ref_cavity_debug(ref_cavity))
    printf(" conforming start %d tris %d segs\n",
           ref_list_n(ref_cavity_tri_list(ref_cavity)),
           ref_cavity_nseg(ref_cavity));

  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "initial seg manifold");

  keep_growing = REF_TRUE;
  while (keep_growing) {
    keep_growing = REF_FALSE;
    each_ref_cavity_valid_seg(ref_cavity, seg) {
      /* skip a seg attached to node */
      if (node == ref_cavity_s2n(ref_cavity, 0, seg) ||
          node == ref_cavity_s2n(ref_cavity, 1, seg))
        continue;

      RSS(ref_cavity_conforming(ref_cavity, seg, &conforming), "free");
      if (!conforming) {
        RSS(ref_cavity_enlarge_seg(ref_cavity, seg), "enlarge seg");
        if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) {
          if (ref_cavity_debug(ref_cavity)) {
            RSS(ref_cavity_tec(ref_cavity, "enlarge.tec"),
                "tec for enlarge_seg fail");
          }
          return REF_SUCCESS;
        }
        keep_growing = REF_TRUE;
      }
    }
  }

  if (ref_cavity_debug(ref_cavity))
    printf(" conforming final %d tris %d segs\n",
           ref_list_n(ref_cavity_tri_list(ref_cavity)),
           ref_cavity_nseg(ref_cavity));

  if (ref_cavity_debug(ref_cavity)) RSS(ref_cavity_topo(ref_cavity), "topo");

  RSS(ref_cavity_manifold(ref_cavity, &manifold), "manifold");
  if (!manifold) {
    if (ref_cavity_debug(ref_cavity)) printf(" conforming not manifold\n");
    ref_cavity_state(ref_cavity) = REF_CAVITY_MANIFOLD_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cavity_state(ref_cavity) = REF_CAVITY_VISIBLE;

  RSS(ref_cavity_verify_seg_manifold(ref_cavity), "final seg manifold");

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
  REF_BOOL visible;
  REF_BOOL keep_growing;
  REF_BOOL manifold;

  RAS(ref_node_owned(ref_node, node), "cavity part must own node");

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge start %d tets %d faces\n",
           ref_list_n(ref_cavity_tet_list(ref_cavity)),
           ref_cavity_nface(ref_cavity));

  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "initial manifold check");

  keep_growing = REF_TRUE;
  while (keep_growing) {
    keep_growing = REF_FALSE;
    each_ref_cavity_valid_face(ref_cavity, face) {
      /* skip a face attached to node */
      if (node == ref_cavity_f2n(ref_cavity, 0, face) ||
          node == ref_cavity_f2n(ref_cavity, 1, face) ||
          node == ref_cavity_f2n(ref_cavity, 2, face))
        continue;

      RSS(ref_cavity_visible(ref_cavity, face, &visible), "free");
      if (!visible) {
        RSS(ref_cavity_enlarge_face(ref_cavity, face), "enlarge face");
        if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) {
          if (ref_cavity_debug(ref_cavity)) {
            RSS(ref_cavity_tec(ref_cavity, "enlarge.tec"),
                "tec for enlarge_face fail");
          }
          return REF_SUCCESS;
        }
        keep_growing = REF_TRUE;
      }
    }
  }

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge final %d tets %d faces\n",
           ref_list_n(ref_cavity_tet_list(ref_cavity)),
           ref_cavity_nface(ref_cavity));

  if (ref_cavity_debug(ref_cavity)) RSS(ref_cavity_topo(ref_cavity), "topo");

  RSS(ref_cavity_manifold(ref_cavity, &manifold), "manifold");
  if (!manifold) {
    if (ref_cavity_debug(ref_cavity)) printf(" conforming not manifold\n");
    ref_cavity_state(ref_cavity) = REF_CAVITY_MANIFOLD_CONSTRAINED;
    return REF_SUCCESS;
  }

  ref_cavity_state(ref_cavity) = REF_CAVITY_VISIBLE;

  RSS(ref_cavity_verify_face_manifold(ref_cavity), "final manifold check");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_check_visible(REF_CAVITY ref_cavity) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_INT face;
  REF_BOOL visible;

  RAS(ref_node_owned(ref_node, node), "cavity part must own node");

  if (REF_CAVITY_UNKNOWN != ref_cavity_state(ref_cavity)) return REF_SUCCESS;

  each_ref_cavity_valid_face(ref_cavity, face) {
    /* skip a face attached to node */
    if (node == ref_cavity_f2n(ref_cavity, 0, face) ||
        node == ref_cavity_f2n(ref_cavity, 1, face) ||
        node == ref_cavity_f2n(ref_cavity, 2, face))
      continue;

    RSS(ref_cavity_visible(ref_cavity, face, &visible), "free");
    if (!visible) {
      ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
      return REF_SUCCESS;
    }
  }

  ref_cavity_state(ref_cavity) = REF_CAVITY_VISIBLE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_seg(REF_CAVITY ref_cavity, REF_INT seg) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_BOOL have_cell0, have_cell1;
  REF_INT ntri, tri_list[2];

  RAS(ref_cavity_valid_seg(ref_cavity, seg), "invalid seg");

  RSS(ref_cell_list_with2(tri, ref_cavity_s2n(ref_cavity, 0, seg),
                          ref_cavity_s2n(ref_cavity, 1, seg), 2, &ntri,
                          tri_list),
      "tri with2");
  if (2 != ntri) { /* for nonmanifold wirebody */
    ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), tri_list[0],
                        &have_cell0),
      "cell0");
  RSS(ref_list_contains(ref_cavity_tri_list(ref_cavity), tri_list[1],
                        &have_cell1),
      "cell1");
  if (have_cell0 == have_cell1) THROW("cavity same seg-tri state");
  if (have_cell0) RSS(ref_cavity_add_tri(ref_cavity, tri_list[1]), "add c1");
  if (have_cell1) RSS(ref_cavity_add_tri(ref_cavity, tri_list[0]), "add c0");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_face(REF_CAVITY ref_cavity, REF_INT face) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT face_node, face_nodes[4], node;
  REF_BOOL have_cell0, have_cell1;
  REF_INT tet0, tet1;

  RAS(ref_cavity_valid_face(ref_cavity, face), "invalid face");

  /* make sure all face nodes are owned */
  each_ref_cavity_face_node(ref_cavity, face_node) {
    node = ref_cavity_f2n(ref_cavity, face_node, face);
    if (!ref_node_owned(ref_node, node)) {
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
        printf("%d face_nodes %d %d %d %d\n%f %f %f\n", face, face_nodes[0],
               face_nodes[1], face_nodes[2], face_nodes[3],
               ref_node_xyz(ref_node, 0, face_nodes[0]),
               ref_node_xyz(ref_node, 1, face_nodes[0]),
               ref_node_xyz(ref_node, 2, face_nodes[0]));
        ref_cavity_tec(ref_cavity, "ref_cavity_too_many_tet.tec");
      });
  if (REF_EMPTY == tet0) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
    return REF_SUCCESS;
  }
  if (REF_EMPTY == tet1) {
    ref_cavity_state(ref_cavity) = REF_CAVITY_BOUNDARY_CONSTRAINED;
    return REF_SUCCESS;
  }

  RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), tet0, &have_cell0),
      "cell0");
  RSS(ref_list_contains(ref_cavity_tet_list(ref_cavity), tet1, &have_cell1),
      "cell1");
  if (have_cell0) RSS(ref_cavity_add_tet(ref_cavity, tet1), "add c1");
  if (have_cell1) RSS(ref_cavity_add_tet(ref_cavity, tet0), "add c0");

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_tec(REF_CAVITY ref_cavity, const char *filename) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_DICT node_dict, face_dict;
  REF_CELL ref_cell;
  REF_INT face, face_node, seg, seg_node;
  REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local;
  REF_DBL xyz_phys[3];
  const char *zonetype;
  FILE *f;

  zonetype = "fetetrahedron";

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine cavity\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_list_item(ref_cavity_tet_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tet_list(ref_cavity), item);
    RSS(ref_dict_store(face_dict, cell, 0), "store");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
    }
  }

  if (0 < ref_dict_n(node_dict) && 0 < ref_dict_n(face_dict)) {
    fprintf(f,
            "zone t=\"old-tet\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
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
  }
  RSS(ref_dict_free(face_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  RSS(ref_dict_store(node_dict, node, 0), "store");
  each_ref_cavity_valid_face(ref_cavity, face) {
    RSS(ref_dict_store(face_dict, face, 0), "store");
    each_ref_cavity_face_node(ref_cavity, face_node) {
      RSS(ref_dict_store(node_dict, ref_cavity_f2n(ref_cavity, face_node, face),
                         0),
          "store");
    }
  }

  if (0 < ref_dict_n(node_dict) && 0 < ref_dict_n(face_dict)) {
    fprintf(f,
            "zone t=\"new-tet\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
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
  }

  RSS(ref_dict_free(face_dict), "free face");
  RSS(ref_dict_free(node_dict), "free nodes");

  zonetype = "fetriangle";

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_list_item(ref_cavity_tri_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tri_list(ref_cavity), item);
    RSS(ref_dict_store(face_dict, cell, 0), "store");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
    }
  }

  if (0 < ref_dict_n(node_dict) && 0 < ref_dict_n(face_dict)) {
    fprintf(f,
            "zone t=\"old-tri\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
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
  }
  RSS(ref_dict_free(face_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&face_dict), "create faces");

  RSS(ref_dict_store(node_dict, node, 0), "store");
  each_ref_cavity_valid_seg(ref_cavity, face) {
    RSS(ref_dict_store(face_dict, face, 0), "store");
    each_ref_cavity_seg_node(ref_cavity, seg_node) {
      RSS(ref_dict_store(node_dict, ref_cavity_s2n(ref_cavity, seg_node, face),
                         0),
          "store");
    }
  }

  if (0 < ref_dict_n(node_dict) && 0 < ref_dict_n(face_dict)) {
    fprintf(f,
            "zone t=\"new-tri\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            ref_dict_n(node_dict), ref_dict_n(face_dict), "point", zonetype);
    for (item = 0; item < ref_dict_n(node_dict); item++) {
      local = ref_dict_key(node_dict, item);
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local);
      fprintf(f, " %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
    }

    for (item = 0; item < ref_dict_n(face_dict); item++) {
      seg = ref_dict_key(face_dict, item);
      RSS(ref_dict_location(node_dict, node, &local), "center node");
      fprintf(f, " %d", local + 1);
      each_ref_cavity_seg_node(ref_cavity, seg_node) {
        RSS(ref_dict_location(
                node_dict, ref_cavity_s2n(ref_cavity, seg_node, seg), &local),
            "ret");
        fprintf(f, " %d", local + 1);
      }
      fprintf(f, "\n");
    }
  }

  RSS(ref_dict_free(face_dict), "free face");
  RSS(ref_dict_free(node_dict), "free nodes");
  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_validate(REF_CAVITY ref_cavity) {
  REF_NODE ref_node = ref_grid_node(ref_cavity_grid(ref_cavity));
  REF_INT node;
  REF_INT face, face_node;
  REF_INT seg, seg_node;

  RAS(ref_node_valid(ref_node, ref_cavity_node(ref_cavity)),
      "cavity node not valid");
  if (REF_EMPTY != ref_cavity_surf_node(ref_cavity))
    RAS(ref_node_valid(ref_node, ref_cavity_surf_node(ref_cavity)),
        "cavity surf node not valid");

  each_ref_cavity_valid_face(ref_cavity, face) {
    ; /* semi to force format */
    each_ref_cavity_face_node(ref_cavity, face_node) {
      node = ref_cavity_f2n(ref_cavity, face_node, face);
      RAS(ref_node_valid(ref_node, node), "cavity face node not valid");
    }
  }

  each_ref_cavity_valid_seg(ref_cavity, seg) {
    ; /* semi to force format */
    each_ref_cavity_seg_node(ref_cavity, seg_node) {
      node = ref_cavity_s2n(ref_cavity, seg_node, seg);
      RAB(ref_node_valid(ref_node, node), "cavity segment node not valid",
          { printf("for node %d\n", node); });
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_ratio(REF_CAVITY ref_cavity, REF_BOOL *allowed) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node = ref_cavity_node(ref_cavity);
  REF_DBL ratio;
  REF_INT face, face_node;
  REF_BOOL skip;

  *allowed = REF_TRUE;

  each_ref_cavity_valid_face(ref_cavity, face) {
    skip = REF_FALSE;
    /* skip a collapsed triangle that in on the boundary of cavity */
    each_ref_cavity_face_node(ref_cavity, face_node) {
      if (node == ref_cavity_f2n(ref_cavity, face_node, face)) {
        skip = REF_TRUE;
      }
    }
    if (skip) continue;
    each_ref_cavity_face_node(ref_cavity, face_node) {
      RSS(ref_node_ratio(ref_node, node,
                         ref_cavity_f2n(ref_cavity, face_node, face), &ratio),
          "ratio");
      if (ratio < ref_grid_adapt(ref_grid, post_min_ratio) ||
          ratio > ref_grid_adapt(ref_grid, post_max_ratio)) {
        *allowed = REF_FALSE;
        return REF_SUCCESS;
      }
    }
  }

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
  each_ref_list_item(ref_cavity_tet_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tet_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "new qual");
    n++;
    total_quality += quality;
    min_quality = MIN(min_quality, quality);
  }
  if (ref_cavity_debug(ref_cavity) && n > 0)
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
    each_ref_cavity_face_node(ref_cavity, face_node) {
      if (node == ref_cavity_f2n(ref_cavity, face_node, face)) {
        skip = REF_TRUE;
      }
    }
    if (skip) continue;
    each_ref_cavity_face_node(ref_cavity, face_node) {
      nodes[face_node] = ref_cavity_f2n(ref_cavity, face_node, face);
    }
    nodes[3] = node;
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "new qual");
    n++;
    total_quality += quality;
    min_quality = MIN(min_quality, quality);
  }
  if (ref_cavity_debug(ref_cavity) && n > 0)
    printf("+ min %12.8f avg %12.8f n %d\n", min_quality,
           total_quality / ((REF_DBL)n), n);
  *min_add = min_quality;
  n_add = n;

  if (ref_cavity_debug(ref_cavity))
    printf(" min %12.8f <- %12.8f diff %12.8f n %d <- %d\n", *min_add, *min_del,
           *min_add - *min_del, n_add, n_del);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_normdev(REF_CAVITY ref_cavity, REF_BOOL *improved) {
  REF_GRID ref_grid = ref_cavity_grid(ref_cavity);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT node;
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normdev, min_normdev, old_normdev;
  REF_DBL sign_uv_area, uv_area, min_uv_area;
  REF_INT seg, seg_node;
  REF_BOOL skip;

  node = ref_cavity_seg_node(ref_cavity);

  *improved = REF_TRUE;

  min_normdev = 2.0;
  min_uv_area = REF_DBL_MAX;
  each_ref_list_item(ref_cavity_tri_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tri_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev),
        "old tri normdev");
    min_normdev = MIN(min_normdev, normdev);
    if (!ref_geom_meshlinked(ref_geom)) {
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
      RSS(ref_geom_uv_area_sign(ref_grid, nodes[3], &sign_uv_area), "sign");
      uv_area *= sign_uv_area;
      min_uv_area = MIN(min_uv_area, uv_area);
    }
  }
  old_normdev = min_normdev;

  if (ref_cavity_debug(ref_cavity))
    printf("- min %12.8f %12.8f\n", min_normdev, min_uv_area);

  min_normdev = 2.0;
  min_uv_area = REF_DBL_MAX;
  each_ref_cavity_valid_seg(ref_cavity, seg) {
    skip = REF_FALSE;
    /* skip a collapsed triangle that in on the boundary of cavity */
    each_ref_cavity_seg_node(ref_cavity, seg_node) {
      if (node == ref_cavity_s2n(ref_cavity, seg_node, seg)) {
        skip = REF_TRUE;
      }
    }
    if (skip) continue;
    each_ref_cavity_seg_node(ref_cavity, seg_node) {
      nodes[seg_node] = ref_cavity_s2n(ref_cavity, seg_node, seg);
    }
    nodes[2] = node;
    nodes[3] = ref_cavity_s2n(ref_cavity, 2, seg);
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev),
        "new tri normdev");
    min_normdev = MIN(min_normdev, normdev);
    if (!ref_geom_meshlinked(ref_geom)) {
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
      RSS(ref_geom_uv_area_sign(ref_grid, nodes[3], &sign_uv_area), "sign");
      uv_area *= sign_uv_area;
      min_uv_area = MIN(min_uv_area, uv_area);
    }
  }
  if (ref_cavity_debug(ref_cavity))
    printf("+ min %12.8f %12.8f\n", min_normdev, min_uv_area);

  *improved = (min_uv_area > ref_node_min_uv_area(ref_node) &&
               min_normdev > old_normdev);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_topo(REF_CAVITY ref_cavity) {
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT item, cell, face, face_node, seg;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_cell = ref_grid_tet(ref_cavity_grid(ref_cavity));
  each_ref_list_item(ref_cavity_tet_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tet_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    printf("old tet %d %d %d %d (%d)\n", nodes[0], nodes[1], nodes[2], nodes[3],
           cell);
  }

  node = ref_cavity_node(ref_cavity);
  each_ref_cavity_valid_face(ref_cavity, face) {
    if (node == ref_cavity_f2n(ref_cavity, 0, face)) continue;
    if (node == ref_cavity_f2n(ref_cavity, 1, face)) continue;
    if (node == ref_cavity_f2n(ref_cavity, 2, face)) continue;
    printf("new tet ");
    for (face_node = 0; face_node < 3; face_node++)
      printf(" %d ", ref_cavity_f2n(ref_cavity, face_node, face));
    printf(" %d ", node);
    printf(" [%d]\n", face);
  }

  ref_cell = ref_grid_tri(ref_cavity_grid(ref_cavity));
  each_ref_list_item(ref_cavity_tri_list(ref_cavity), item) {
    cell = ref_list_value(ref_cavity_tri_list(ref_cavity), item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell");
    printf("old tri %d %d %d faceid %d (%d)\n", nodes[0], nodes[1], nodes[2],
           nodes[3], cell);
  }

  node = ref_cavity_seg_node(ref_cavity);
  each_ref_cavity_valid_seg(ref_cavity, seg) {
    if (node == ref_cavity_s2n(ref_cavity, 0, seg)) continue;
    if (node == ref_cavity_s2n(ref_cavity, 1, seg)) continue;
    printf("new tri ");
    printf(" %d ", ref_cavity_s2n(ref_cavity, 0, seg));
    printf(" %d ", ref_cavity_s2n(ref_cavity, 1, seg));
    printf(" %d ", node);
    printf(" faceid %d ", ref_cavity_s2n(ref_cavity, 2, seg));
    printf(" [%d]\n", seg);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_cavity_edge_swap_boundary(REF_GRID ref_grid,
                                                REF_INT node0, REF_INT node1,
                                                REF_BOOL *allowed) {
  REF_BOOL has_triangle, same_face, manifold, conforming, improved, in_limits;

  *allowed = REF_FALSE;

  RSS(ref_cell_has_side(ref_grid_tri(ref_grid), node0, node1, &has_triangle),
      "triangle side");
  if (!has_triangle) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_same_faceid(ref_grid, node0, node1, &same_face),
      "not allowed if a side of a edge or diff faceid");
  if (!same_face) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_manifold(ref_grid, node0, node1, &manifold),
      "not allowed if creates a triangle or edge that exists");
  if (!manifold) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_conforming(ref_grid, node0, node1, &conforming),
      "normals and uv area must conform to geom");
  if (!conforming) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_quality(ref_grid, node0, node1, &improved),
      "require tri quality improvement");
  if (!improved) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  RSS(ref_swap_ratio(ref_grid, node0, node1, &in_limits),
      "require tri ratio in limits");
  if (!in_limits) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  *allowed = REF_TRUE;
  return REF_SUCCESS;
}

static REF_STATUS ref_cavity_swap_tet_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, min_del, min_add, best;
  REF_INT best_other;
  REF_CAVITY ref_cavity;
  REF_BOOL allowed;
  REF_INT degree;
  REF_INT n0, n1, n2;
  REF_INT other;
  REF_INT others[12][3] = {
      {0, 1, 2}, {0, 1, 3}, {0, 2, 1}, {0, 2, 3}, {0, 3, 1}, {0, 3, 2},
      {1, 2, 0}, {1, 2, 3}, {1, 3, 0}, {1, 3, 2}, {2, 3, 0}, {2, 3, 1},
  };

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    if (quality < ref_grid_adapt(ref_grid, swap_min_quality)) {
      best_other = REF_EMPTY;
      best = -2.0;
      for (other = 0; other < 12; other++) {
        n0 = others[other][0];
        n1 = others[other][1];
        n2 = others[other][2];
        RSS(ref_cell_local_gem(ref_cell, ref_node, nodes[n0], nodes[n1],
                               &allowed),
            "local gem");
        if (!allowed) continue;
        RSS(ref_cavity_edge_swap_boundary(ref_grid, nodes[n0], nodes[n1],
                                          &allowed),
            "surface geom and topo");
        if (!allowed) continue;
        RSS(ref_cell_degree_with2(ref_cell, nodes[n0], nodes[n1], &degree),
            "edge degree");
        if (degree > ref_grid_adapt(ref_grid, swap_max_degree)) continue;
        RSS(ref_cavity_create(&ref_cavity), "create");
        if (REF_SUCCESS != ref_cavity_form_edge_swap(ref_cavity, ref_grid,
                                                     nodes[n0], nodes[n1],
                                                     nodes[n2])) {
          REF_WHERE("form edge swap"); /* note but skip cavity failures */
          RSS(ref_cavity_free(ref_cavity), "free");
          continue;
        }
        if (REF_CAVITY_INCONSISTENT == ref_cavity_state(ref_cavity)) {
          /* skip cavity failures */
          RSS(ref_cavity_free(ref_cavity), "free");
          continue;
        }
        if (REF_SUCCESS != ref_cavity_check_visible(ref_cavity)) {
          REF_WHERE("check visible"); /* note but skip cavity failures */
          RSS(ref_cavity_free(ref_cavity), "free");
          continue;
        }
        if (REF_CAVITY_VISIBLE == ref_cavity_state(ref_cavity)) {
          RSS(ref_cavity_ratio(ref_cavity, &allowed), "post ratio limits");
          if (!allowed) {
            RSS(ref_cavity_free(ref_cavity), "free");
            continue;
          }
          RSS(ref_cavity_change(ref_cavity, &min_del, &min_add), "change");
          if (min_add - min_del > 0.0001) {
            if (best < min_add) {
              best = min_add;
              best_other = other;
            }
          }
        }
        RSS(ref_cavity_free(ref_cavity), "free");
      }
      if (REF_EMPTY != best_other) {
        RSS(ref_cavity_create(&ref_cavity), "create");
        n0 = others[best_other][0];
        n1 = others[best_other][1];
        n2 = others[best_other][2];
        RSS(ref_cavity_form_edge_swap(ref_cavity, ref_grid, nodes[n0],
                                      nodes[n1], nodes[n2]),
            "cavity gem");
        RSS(ref_cavity_check_visible(ref_cavity), "enlarge viz");
        if (ref_cavity_debug(ref_cavity)) {
          RSS(ref_cavity_change(ref_cavity, &min_del, &min_add), "change");
          printf("cavity accepted %f -> %f\n", min_del, min_add);
        }
        RSS(ref_cavity_replace(ref_cavity), "replace");
        RSS(ref_cavity_free(ref_cavity), "free");
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_cavity_surf_geom_edge_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL edg = ref_grid_edg(ref_grid);
  REF_INT node0, node1, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell;
  REF_INT edge_tri[2];
  REF_INT i, tri_cell;
  REF_DBL normdev;
  REF_CAVITY ref_cavity;
  REF_BOOL improved;

  if (!ref_grid_surf(ref_grid)) return REF_SUCCESS;

  each_ref_cell_valid_cell_with_nodes(edg, cell, nodes) {
    node0 = nodes[0];
    node1 = nodes[1];

    /* expands around node0, should try node1 too? */

    /* skip unless node0 is owned */
    if (!ref_node_owned(ref_node, node0)) {
      continue;
    }

    RSB(ref_cell_list_with2(tri, node0, node1, 2, &ncell, edge_tri), "tris", {
      REF_DBL xyz_phys[3];
      REF_INT local;
      local = node0;
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local);
      printf(" %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
      local = node1;
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local);
      printf(" %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
      ref_export_tec_surf(ref_grid, "ref_cavity_geom_edge_fail.tec");
    });
    for (i = 0; i < ncell; i++) {
      tri_cell = edge_tri[i];
      RSS(ref_cell_nodes(tri, tri_cell, nodes), "cell nodes");
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "nd");
      if (normdev < 0.5) {
        RSS(ref_cavity_create(&ref_cavity), "create");
        RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node0), "insert ball");
        RSS(ref_cavity_add_tri(ref_cavity, tri_cell), "insert tri");
        RSS(ref_cavity_enlarge_conforming(ref_cavity), "enlarge tri");
        if (REF_CAVITY_VISIBLE == ref_cavity_state(ref_cavity)) {
          RSS(ref_cavity_normdev(ref_cavity, &improved), "normdev tri");
          if (improved) {
            RSS(ref_cavity_replace(ref_cavity), "replace tri");
          }
        }
        RSS(ref_cavity_free(ref_cavity), "free");
      }
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_cavity_surf_geom_face_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normdev;
  REF_CAVITY ref_cavity;
  REF_BOOL improved, geom_edge;
  if (!ref_grid_surf(ref_grid)) return REF_SUCCESS;
  each_ref_cell_valid_cell_with_nodes(tri, cell, nodes) {
    if (!ref_node_owned(ref_node, nodes[0])) {
      continue;
    }
    RSS(ref_geom_is_a(ref_grid_geom(ref_grid), nodes[0], REF_GEOM_EDGE,
                      &geom_edge),
        "n");
    if (geom_edge) {
      continue;
    }
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "nd");
    if (normdev < 0.1) {
      RSS(ref_cavity_create(&ref_cavity), "create");
      RSS(ref_cavity_form_ball(ref_cavity, ref_grid, nodes[0]), "insert ball");
      RSS(ref_cavity_enlarge_conforming(ref_cavity), "enlarge tri");
      if (REF_CAVITY_VISIBLE == ref_cavity_state(ref_cavity)) {
        RSS(ref_cavity_normdev(ref_cavity, &improved), "normdev tri");
        if (improved) {
          RSS(ref_cavity_replace(ref_cavity), "replace tri");
        }
      }
      RSS(ref_cavity_free(ref_cavity), "free");
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_pass(REF_GRID ref_grid) {
  RSS(ref_cavity_swap_tet_pass(ref_grid), "cavity swap pass");
  RSS(ref_cavity_surf_geom_edge_pass(ref_grid), "cavity geom edge");
  RSS(ref_cavity_surf_geom_face_pass(ref_grid), "cavity geom edge");
  return REF_SUCCESS;
}
