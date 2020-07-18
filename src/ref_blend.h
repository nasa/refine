
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

BEGIN_C_DECLORATION

struct REF_BLEND_STRUCT {
  REF_GRID grid;
  REF_DBL *displacement;
};

REF_STATUS ref_blend_create(REF_BLEND *ref_blend, REF_GRID ref_grid);
REF_STATUS ref_blend_free(REF_BLEND ref_blend);

END_C_DECLORATION

#endif /* REF_BLEND_H */
