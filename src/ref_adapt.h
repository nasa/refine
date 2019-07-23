
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

#ifndef REF_ADAPT_H
#define REF_ADAPT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ADAPT_STRUCT REF_ADAPT_STRUCT;
typedef REF_ADAPT_STRUCT *REF_ADAPT;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_ADAPT_STRUCT {
  REF_INT split_per_pass;
  REF_BOOL split_ratio_growth;
  REF_DBL split_ratio;
  REF_DBL split_quality_absolute;
  REF_DBL split_quality_relative;

  REF_INT collapse_per_pass;
  REF_DBL collapse_ratio;
  REF_DBL collapse_quality_absolute;

  REF_INT smooth_per_pass;
  REF_DBL smooth_min_quality;

  REF_INT swap_max_degree;
  REF_DBL swap_min_quality;

  REF_DBL post_min_normdev;
  REF_DBL post_min_ratio;
  REF_DBL post_max_ratio;

  REF_BOOL instrument;
  REF_BOOL watch_param;
};

REF_STATUS ref_adapt_create(REF_ADAPT *ref_adapt);
REF_STATUS ref_adapt_deep_copy(REF_ADAPT *ref_adapt_ptr, REF_ADAPT original);
REF_STATUS ref_adapt_free(REF_ADAPT ref_adapt);

REF_STATUS ref_adapt_pass(REF_GRID ref_grid, REF_BOOL *all_done);

REF_STATUS ref_adapt_surf_to_geom(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_ADAPT_H */
