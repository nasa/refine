
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

#ifndef REF_SHARD_H
#define REF_SHARD_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_SHARD_STRUCT REF_SHARD_STRUCT;
typedef REF_SHARD_STRUCT *REF_SHARD;
END_C_DECLORATION

#include "ref_face.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_SHARD_STRUCT {
  REF_GRID grid;
  REF_FACE face;
  REF_INT *mark;
};

#define ref_shard_grid(ref_shard) ((ref_shard)->grid)
#define ref_shard_face(ref_shard) ((ref_shard)->face)

REF_STATUS ref_shard_create(REF_SHARD *ref_shard, REF_GRID ref_grid);
REF_STATUS ref_shard_free(REF_SHARD ref_shard);

#define ref_shard_mark(ref_shard, face) ((ref_shard)->mark[face])

REF_STATUS ref_shard_mark_to_split(REF_SHARD ref_shard, REF_INT node0,
                                   REF_INT node1);
REF_STATUS ref_shard_marked(REF_SHARD ref_shard, REF_INT node0, REF_INT node1,
                            REF_BOOL *marked);

REF_STATUS ref_shard_mark_n(REF_SHARD ref_shard, REF_INT *face_marks,
                            REF_INT *hex_marks);

REF_STATUS ref_shard_mark_cell_edge_split(REF_SHARD ref_shard, REF_INT cell,
                                          REF_INT cell_edge);

REF_STATUS ref_shard_mark_relax(REF_SHARD ref_shard);
REF_STATUS ref_shard_split(REF_SHARD ref_shard);

REF_STATUS ref_shard_prism_into_tet(REF_GRID ref_grid, REF_INT keeping_n_layers,
                                    REF_INT of_faceid);

REF_STATUS ref_shard_extract_tri(REF_GRID ref_grid, REF_CELL *ref_cell);

END_C_DECLORATION

#endif /* REF_SHARD_H */
