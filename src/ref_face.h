
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

#ifndef REF_FACE_H
#define REF_FACE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_FACE_STRUCT REF_FACE_STRUCT;
typedef REF_FACE_STRUCT *REF_FACE;
END_C_DECLORATION

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

struct REF_FACE_STRUCT {
  REF_INT n, max;
  REF_INT *f2n;
  REF_ADJ adj;
};

REF_STATUS ref_face_create(REF_FACE *ref_face, REF_GRID ref_grid);
REF_STATUS ref_face_free(REF_FACE ref_face);

#define ref_face_n(ref_face) ((ref_face)->n)
#define ref_face_max(ref_face) ((ref_face)->max)

#define ref_face_f2n(ref_face, node, face) \
  ((ref_face)->f2n[(node) + 4 * (face)])

#define ref_face_adj(ref_face) ((ref_face)->adj)

#define each_ref_face(ref_face, face) \
  for ((face) = 0; (face) < ref_face_n(ref_face); (face)++)

REF_STATUS ref_face_inspect(REF_FACE ref_face);

REF_STATUS ref_face_with(REF_FACE ref_face, REF_INT *nodes, REF_INT *face);
REF_STATUS ref_face_spanning(REF_FACE ref_face, REF_INT node0, REF_INT node1,
                             REF_INT *face);

REF_STATUS ref_face_add_uniquely(REF_FACE ref_face, REF_INT *nodes);

REF_STATUS ref_face_normal(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz2,
                           REF_DBL *xyz3, REF_DBL *normal);

REF_STATUS ref_face_open_node(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz2,
                              REF_DBL *xyz3, REF_INT *open_node, REF_DBL *dot);

REF_STATUS ref_face_part(REF_FACE ref_face, REF_NODE ref_node, REF_INT face,
                         REF_INT *part);

END_C_DECLORATION

#endif /* REF_FACE_H */
