
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

#include "ref_blend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_egads.h"
#include "ref_malloc.h"

#define ref_blend_grid(ref_blend) ((ref_blend)->grid)
#define ref_blend_displacement(ref_blend, ixyz, geom) \
  ((ref_blend)->displacement[(ixyz) + 3 * (geom)])

REF_STATUS ref_blend_create(REF_BLEND *ref_blend_ptr, REF_GRID ref_grid) {
  REF_BLEND ref_blend;
  REF_INT n;

  ref_malloc(*ref_blend_ptr, 1, REF_BLEND_STRUCT);

  ref_blend = *ref_blend_ptr;

  RSS(ref_grid_deep_copy(&ref_blend_grid(ref_blend), ref_grid), "deep copy");
  n = ref_geom_max(ref_grid_geom(ref_blend_grid(ref_blend)));
  ref_malloc_init(ref_blend->displacement, 3 * n, REF_DBL, 0.0);

  return REF_SUCCESS;
}

REF_STATUS ref_blend_free(REF_BLEND ref_blend) {
  if (NULL == (void *)ref_blend) return REF_NULL;

  ref_free(ref_blend->displacement);
  ref_grid_free(ref_blend_grid(ref_blend));
  ref_free(ref_blend);

  return REF_SUCCESS;
}

REF_STATUS ref_blend_initialize(REF_BLEND ref_blend) {
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DBL edge_xyz[3], face_xyz[3];
  REF_INT i, edge_geom, face_geom, node, item;

  /* also displace edge and face to geom nodes */

  each_ref_geom_edge(ref_geom, edge_geom) {
    RSS(ref_egads_eval(ref_geom, edge_geom, edge_xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    each_ref_geom_having_node(ref_geom, node, item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, node)) {
        RSS(ref_egads_eval(ref_geom, face_geom, face_xyz, NULL), "eval face");
        for (i = 0; i < 3; i++) {
          ref_blend_displacement(ref_blend, i, face_geom) =
              edge_xyz[i] - face_xyz[i];
        }
      }
    }
  }
  return REF_SUCCESS;
}
