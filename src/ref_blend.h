
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

#ifndef REF_BLEND_H
#define REF_BLEND_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_BLEND_STRUCT REF_BLEND_STRUCT;
typedef REF_BLEND_STRUCT *REF_BLEND;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_search.h"

BEGIN_C_DECLORATION

struct REF_BLEND_STRUCT {
  REF_GRID grid;
  REF_DBL *displacement;
  REF_BOOL *strong_bc;
  REF_SEARCH *edge_search;
  REF_SEARCH *face_search;
};

#define ref_blend_grid(ref_blend) ((ref_blend)->grid)

REF_STATUS ref_blend_create(REF_BLEND *ref_blend, REF_GRID freeable_ref_grid);
REF_STATUS ref_blend_free(REF_BLEND ref_blend);

REF_STATUS ref_blend_attach(REF_GRID ref_grid);
REF_STATUS ref_blend_import(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_blend_enclosing(REF_BLEND ref_blend, REF_INT type, REF_INT id,
                               REF_DBL *param, REF_INT *cell, REF_DBL *bary);

REF_STATUS ref_blend_eval_at(REF_BLEND ref_blend, REF_INT type, REF_INT id,
                             REF_DBL *params, REF_DBL *xyz, REF_DBL *dxyz_dtuv);
REF_STATUS ref_blend_inverse_eval(REF_BLEND ref_blend, REF_INT type, REF_INT id,
                                  REF_DBL *xyz, REF_DBL *param);

REF_STATUS ref_blend_tec(REF_BLEND ref_blend, const char *filename);

REF_STATUS ref_blend_max_distance(REF_BLEND ref_blend, REF_DBL *distance);

END_C_DECLORATION

#endif /* REF_BLEND_H */
