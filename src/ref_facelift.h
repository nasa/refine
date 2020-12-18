
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

#ifndef REF_FACELIFT_H
#define REF_FACELIFT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_FACELIFT_STRUCT REF_FACELIFT_STRUCT;
typedef REF_FACELIFT_STRUCT *REF_FACELIFT;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_search.h"

BEGIN_C_DECLORATION

struct REF_FACELIFT_STRUCT {
  REF_GRID grid;
  REF_DBL *displacement;
  REF_BOOL *strong_bc;
  REF_SEARCH *edge_search;
  REF_SEARCH *face_search;
};

#define ref_facelift_grid(ref_facelift) ((ref_facelift)->grid)

#define ref_facelift_direct(ref_facelift) (NULL == (ref_facelift)->displacement)

#define ref_facelift_displacement(ref_facelift, ixyz, geom) \
  ((ref_facelift)->displacement[(ixyz) + 3 * (geom)])
#define ref_facelift_distance(ref_facelift, geom)                  \
  (sqrt(pow(ref_facelift_displacement(ref_facelift, 0, geom), 2) + \
        pow(ref_facelift_displacement(ref_facelift, 1, geom), 2) + \
        pow(ref_facelift_displacement(ref_facelift, 2, geom), 2)))

REF_STATUS ref_facelift_create(REF_FACELIFT *ref_facelift,
                               REF_GRID freeable_ref_grid, REF_BOOL direct);
REF_STATUS ref_facelift_free(REF_FACELIFT ref_facelift);

REF_STATUS ref_facelift_attach(REF_GRID ref_grid);
REF_STATUS ref_facelift_import(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_facelift_enclosing(REF_FACELIFT ref_facelift, REF_INT type,
                                  REF_INT id, REF_DBL *param, REF_INT *cell,
                                  REF_DBL *bary);

REF_STATUS ref_facelift_eval_at(REF_FACELIFT ref_facelift, REF_INT type,
                                REF_INT id, REF_DBL *params, REF_DBL *xyz,
                                REF_DBL *dxyz_dtuv);
REF_STATUS ref_facelift_inverse_eval(REF_FACELIFT ref_facelift, REF_INT type,
                                     REF_INT id, REF_DBL *xyz, REF_DBL *param);

REF_STATUS ref_facelift_tec(REF_FACELIFT ref_facelift, const char *filename);

REF_STATUS ref_facelift_max_distance(REF_FACELIFT ref_facelift,
                                     REF_DBL *distance);
REF_STATUS ref_facelift_multiscale(REF_GRID ref_grid, REF_DBL complexity);

END_C_DECLORATION

#endif /* REF_FACELIFT_H */
