
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

#ifndef REF_CAVITY_H
#define REF_CAVITY_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CAVITY_STRUCT REF_CAVITY_STRUCT;
typedef REF_CAVITY_STRUCT *REF_CAVITY;
typedef enum REF_CAVITY_STATES {
  /* 0 */ REF_CAVITY_UNKNOWN,
  /* 1 */ REF_CAVITY_VISIBLE,
  /* 2 */ REF_CAVITY_BOUNDARY_CONSTRAINED,
  /* 3 */ REF_CAVITY_PARTITION_CONSTRAINED,
  /* 4 */ REF_CAVITY_MANIFOLD_CONSTRAINED,
  /* 5 */ REF_CAVITY_STATUS_LAST
} REF_CAVITY_STATE;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_list.h"

BEGIN_C_DECLORATION

struct REF_CAVITY_STRUCT {
  REF_CAVITY_STATE state;
  REF_GRID ref_grid;
  REF_INT node;
  REF_INT nseg;
  REF_INT maxseg;
  REF_INT blankseg;
  REF_INT *s2n;
  REF_INT nface;
  REF_INT maxface;
  REF_INT blankface;
  REF_INT *f2n;
  REF_LIST tri_list;
  REF_LIST tet_list;
  REF_INT faceid;
  REF_BOOL debug;
};

REF_STATUS ref_cavity_create(REF_CAVITY *ref_cavity);
REF_STATUS ref_cavity_free(REF_CAVITY ref_cavity);
REF_STATUS ref_cavity_inspect(REF_CAVITY ref_cavity);

#define ref_cavity_state(ref_cavity) ((ref_cavity)->state)
#define ref_cavity_node(ref_cavity) ((ref_cavity)->node)
#define ref_cavity_grid(ref_cavity) ((ref_cavity)->ref_grid)

#define ref_cavity_nseg(ref_cavity) ((ref_cavity)->nseg)
#define ref_cavity_maxseg(ref_cavity) ((ref_cavity)->maxseg)
#define ref_cavity_blankseg(ref_cavity) ((ref_cavity)->blankseg)
#define ref_cavity_s2n(ref_cavity, node, cavity) \
  ((ref_cavity)->s2n[(node) + 3 * (cavity)])

#define ref_cavity_nface(ref_cavity) ((ref_cavity)->nface)
#define ref_cavity_maxface(ref_cavity) ((ref_cavity)->maxface)
#define ref_cavity_blankface(ref_cavity) ((ref_cavity)->blankface)
#define ref_cavity_f2n(ref_cavity, node, cavity) \
  ((ref_cavity)->f2n[(node) + 3 * (cavity)])

#define ref_cavity_tri_list(ref_cavity) ((ref_cavity)->tri_list)
#define ref_cavity_tet_list(ref_cavity) ((ref_cavity)->tet_list)

#define ref_cavity_debug(ref_cavity) ((ref_cavity)->debug)

#define ref_cavity_valid_seg(ref_cavity, seg)             \
  ((seg) >= 0 && (seg) < ref_cavity_maxseg(ref_cavity) && \
   REF_EMPTY != ref_cavity_s2n(ref_cavity, 0, seg))

#define each_ref_cavity_valid_seg(ref_cavity, seg)                \
  for ((seg) = 0; (seg) < ref_cavity_maxseg(ref_cavity); (seg)++) \
    if (ref_cavity_valid_seg(ref_cavity, seg))

#define each_ref_cavity_seg_node(ref_cavity, seg_node) \
  for ((seg_node) = 0; (seg_node) < 2; (seg_node)++)

#define ref_cavity_valid_face(ref_cavity, face)              \
  ((face) >= 0 && (face) < ref_cavity_maxface(ref_cavity) && \
   REF_EMPTY != ref_cavity_f2n(ref_cavity, 0, face))

#define each_ref_cavity_valid_face(ref_cavity, face)                  \
  for ((face) = 0; (face) < ref_cavity_maxface(ref_cavity); (face)++) \
    if (ref_cavity_valid_face(ref_cavity, face))

#define each_ref_cavity_face_node(ref_cavity, face_node) \
  for ((face_node) = 0; (face_node) < 3; (face_node)++)

REF_STATUS ref_cavity_insert_seg(REF_CAVITY ref_cavity, REF_INT *nodes);
REF_STATUS ref_cavity_delete_seg(REF_CAVITY ref_cavity, REF_INT seg);
REF_STATUS ref_cavity_find_seg(REF_CAVITY ref_cavity, REF_INT *nodes,
                               REF_INT *found_seg, REF_BOOL *reversed);
REF_STATUS ref_cavity_insert_face(REF_CAVITY ref_cavity, REF_INT *nodes);
REF_STATUS ref_cavity_find_face(REF_CAVITY ref_cavity, REF_INT *nodes,
                                REF_INT *found_face, REF_BOOL *reversed);

REF_STATUS ref_cavity_add_tri(REF_CAVITY ref_cavity, REF_INT tri);
REF_STATUS ref_cavity_replace_tri(REF_CAVITY ref_cavity);

REF_STATUS ref_cavity_add_tet(REF_CAVITY ref_cavity, REF_INT tet);
REF_STATUS ref_cavity_rm_tet(REF_CAVITY ref_cavity, REF_INT tet);
REF_STATUS ref_cavity_replace_tet(REF_CAVITY ref_cavity);

REF_STATUS ref_cavity_form_empty(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                 REF_INT node);
REF_STATUS ref_cavity_form_ball(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                REF_INT node);
REF_STATUS ref_cavity_form_gem(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                               REF_INT node0, REF_INT node1, REF_INT node);
REF_STATUS ref_cavity_form_edge_split(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                      REF_INT node0, REF_INT node1,
                                      REF_INT new_node);
REF_STATUS ref_cavity_form_surf_ball(REF_CAVITY ref_cavity, REF_GRID ref_grid,
                                     REF_INT node);
REF_STATUS ref_cavity_form_surf_edge_split(REF_CAVITY ref_cavity,
                                           REF_GRID ref_grid, REF_INT node0,
                                           REF_INT node1, REF_INT new_node);

REF_STATUS ref_cavity_manifold(REF_CAVITY ref_cavity, REF_BOOL *manifold);
REF_STATUS ref_cavity_conforming(REF_CAVITY ref_cavity, REF_INT seg,
                                 REF_BOOL *conforming);
REF_STATUS ref_cavity_enlarge_conforming(REF_CAVITY ref_cavity);

REF_STATUS ref_cavity_visible(REF_CAVITY ref_cavity, REF_INT face,
                              REF_BOOL *visible);
REF_STATUS ref_cavity_enlarge_visible(REF_CAVITY ref_cavity);

REF_STATUS ref_cavity_enlarge_seg(REF_CAVITY ref_cavity, REF_INT seg);

REF_STATUS ref_cavity_enlarge_face(REF_CAVITY ref_cavity, REF_INT face);

REF_STATUS ref_cavity_tec(REF_CAVITY ref_cavity, const char *filename);

REF_STATUS ref_cavity_local(REF_CAVITY ref_cavity, REF_BOOL *local);
REF_STATUS ref_cavity_change(REF_CAVITY ref_cavity, REF_DBL *min_del,
                             REF_DBL *min_add);
REF_STATUS ref_cavity_normdev(REF_CAVITY ref_cavity, REF_BOOL *improved);
REF_STATUS ref_cavity_topo(REF_CAVITY ref_cavity);

REF_STATUS ref_cavity_pass(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_CAVITY_H */
