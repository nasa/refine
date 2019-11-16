
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

#include "ref_cloud.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"

REF_STATUS ref_cloud_create(REF_CLOUD *ref_cloud_ptr, REF_INT naux) {
  REF_CLOUD ref_cloud;

  ref_malloc(*ref_cloud_ptr, 1, REF_CLOUD_STRUCT);

  ref_cloud = (*ref_cloud_ptr);

  ref_cloud_n(ref_cloud) = 0;
  ref_cloud_max(ref_cloud) = 20;
  ref_cloud_naux(ref_cloud) = naux;

  ref_malloc(ref_cloud->global, ref_cloud_max(ref_cloud), REF_GLOB);
  ref_malloc(ref_cloud->aux,
             ref_cloud_naux(ref_cloud) * ref_cloud_max(ref_cloud), REF_DBL);

  return REF_SUCCESS;
}

REF_STATUS ref_cloud_free(REF_CLOUD ref_cloud) {
  if (NULL == (void *)ref_cloud) return REF_NULL;
  ref_free(ref_cloud->aux);
  ref_free(ref_cloud->global);
  ref_free(ref_cloud);
  return REF_SUCCESS;
}

REF_STATUS ref_cloud_deep_copy(REF_CLOUD *ref_cloud_ptr, REF_CLOUD original) {
  REF_CLOUD ref_cloud;
  REF_INT item, i;
  REF_GLOB global;

  ref_malloc(*ref_cloud_ptr, 1, REF_CLOUD_STRUCT);

  ref_cloud = (*ref_cloud_ptr);

  ref_cloud_n(ref_cloud) = ref_cloud_n(original);
  ref_cloud_max(ref_cloud) = ref_cloud_max(original);
  ref_cloud_naux(ref_cloud) = ref_cloud_naux(original);

  ref_malloc(ref_cloud->global, ref_cloud_max(ref_cloud), REF_GLOB);
  ref_malloc(ref_cloud->aux,
             ref_cloud_naux(ref_cloud) * ref_cloud_max(ref_cloud), REF_DBL);

  each_ref_cloud_global(original, item, global) {
    ref_cloud_global(ref_cloud, item) = global;
    each_ref_cloud_aux(ref_cloud, i) {
      ref_cloud_aux(ref_cloud, i, item) = ref_cloud_aux(original, i, item);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cloud_store(REF_CLOUD ref_cloud, REF_GLOB global, REF_DBL *aux) {
  REF_INT item, i, insert_point;

  if (ref_cloud_max(ref_cloud) == ref_cloud_n(ref_cloud)) {
    ref_cloud_max(ref_cloud) += 100;

    ref_realloc(ref_cloud->global, ref_cloud_max(ref_cloud), REF_GLOB);
    ref_realloc(ref_cloud->aux,
                ref_cloud_naux(ref_cloud) * ref_cloud_max(ref_cloud), REF_DBL);
  }

  insert_point = 0;
  for (item = ref_cloud_n(ref_cloud) - 1; item >= 0; item--) {
    if (ref_cloud_global(ref_cloud, item) == global) { /* replace aux */
      each_ref_cloud_aux(ref_cloud, i) {
        ref_cloud_aux(ref_cloud, i, item) = aux[i];
      }
      return REF_SUCCESS;
    }
    if (ref_cloud_global(ref_cloud, item) < global) {
      insert_point = item + 1;
      break;
    }
  }
  /* shift to open up insert_point */
  for (item = ref_cloud_n(ref_cloud); item > insert_point; item--) {
    ref_cloud_global(ref_cloud, item) = ref_cloud_global(ref_cloud, item - 1);
  }
  for (item = ref_cloud_n(ref_cloud); item > insert_point; item--) {
    each_ref_cloud_aux(ref_cloud, i) {
      ref_cloud_aux(ref_cloud, i, item) = ref_cloud_aux(ref_cloud, i, item - 1);
    }
  }
  /* fill insert_point */
  ref_cloud_n(ref_cloud)++;
  item = insert_point;
  ref_cloud_global(ref_cloud, item) = global;
  each_ref_cloud_aux(ref_cloud, i) {
    ref_cloud_aux(ref_cloud, i, item) = aux[i];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cloud_item(REF_CLOUD ref_cloud, REF_GLOB global, REF_INT *item) {
  REF_INT i;

  *item = REF_EMPTY;

  each_ref_cloud_item(ref_cloud, i) {
    if (global == ref_cloud_global(ref_cloud, i)) {
      *item = i;
      return REF_SUCCESS;
    }
  }

  return REF_NOT_FOUND;
}

REF_BOOL ref_cloud_has_global(REF_CLOUD ref_cloud, REF_GLOB global) {
  REF_INT i;

  each_ref_cloud_item(ref_cloud, i) {
    if (global == ref_cloud_global(ref_cloud, i)) {
      return REF_TRUE;
    }
  }

  return REF_FALSE;
}
