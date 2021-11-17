
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

#ifndef REF_CLOUD_H
#define REF_CLOUD_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CLOUD_STRUCT REF_CLOUD_STRUCT;
typedef REF_CLOUD_STRUCT *REF_CLOUD;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_CLOUD_STRUCT {
  REF_INT n, max, naux;
  REF_GLOB *global;
  REF_DBL *aux;
};

REF_STATUS ref_cloud_create(REF_CLOUD *ref_cloud, REF_INT naux);
REF_STATUS ref_cloud_free(REF_CLOUD ref_cloud);

REF_STATUS ref_cloud_deep_copy(REF_CLOUD *ref_cloud, REF_CLOUD original);

#define ref_cloud_n(ref_cloud) ((ref_cloud)->n)
#define ref_cloud_max(ref_cloud) ((ref_cloud)->max)
#define ref_cloud_naux(ref_cloud) ((ref_cloud)->naux)
#define ref_cloud_global(ref_cloud, item) ((ref_cloud)->global[(item)])
#define ref_cloud_aux(ref_cloud, aux_index, item) \
  ((ref_cloud)->aux[(aux_index) + ref_cloud_naux(ref_cloud) * (item)])

#define ref_cloud_valid(ref_cloud, item) \
  (0 <= (item) && (item) < ref_cloud_n(ref_cloud))

#define ref_cloud_safe_global(ref_cloud, item) \
  (ref_cloud_valid(ref_cloud, item) ? (ref_cloud)->global[(item)] : REF_EMPTY)

#define each_ref_cloud_item(ref_cloud, item) \
  for ((item) = 0; (item) < ref_cloud_n(ref_cloud); (item)++)

#define each_ref_cloud_aux(ref_cloud, aux_index) \
  for ((aux_index) = 0; (aux_index) < ref_cloud_naux(ref_cloud); (aux_index)++)

#define each_ref_cloud_global(ref_cloud, item, global)                \
  for ((item) = 0, (global) = ref_cloud_safe_global(ref_cloud, item); \
       (item) < ref_cloud_n(ref_cloud);                               \
       (item)++, (global) = ref_cloud_safe_global(ref_cloud, item))

REF_STATUS ref_cloud_store(REF_CLOUD ref_cloud, REF_GLOB global, REF_DBL *aux);
REF_STATUS ref_cloud_push(REF_CLOUD ref_cloud, REF_GLOB global, REF_DBL *aux);
REF_STATUS ref_cloud_item(REF_CLOUD ref_cloud, REF_GLOB global, REF_INT *item);
REF_BOOL ref_cloud_has_global(REF_CLOUD ref_cloud, REF_GLOB global);

END_C_DECLORATION

#endif /* REF_CLOUD_H */
