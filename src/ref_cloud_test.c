
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  REF_CLOUD ref_cloud;

  {
    REIS(REF_NULL, ref_cloud_free(NULL), "dont free NULL");
    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    REIS(0, ref_cloud_n(ref_cloud), "init zero");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* missing fails safely */
    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    REIS(REF_EMPTY, ref_cloud_safe_global(ref_cloud, -1), "missing");
    REIS(REF_EMPTY, ref_cloud_safe_global(ref_cloud, 0), "missing");
    REIS(REF_EMPTY, ref_cloud_safe_global(ref_cloud, 5), "missing");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store one */
    REF_GLOB global = 2;
    REF_DBL aux[4] = {0.0, 1.0, 2.0, 3.0};
    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    REIS(global, ref_cloud_safe_global(ref_cloud, 0), "missing");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store two, recall, and has (in global order) */
    REF_GLOB global;
    REF_INT item;
    REF_DBL aux[4] = {0.0, 1.0, 2.0, 3.0};
    RSS(ref_cloud_create(&ref_cloud, 4), "create");

    global = 2;
    RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    REIS(2, ref_cloud_safe_global(ref_cloud, 0), "missing");

    global = 1;
    RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    REIS(1, ref_cloud_safe_global(ref_cloud, 0), "missing");
    REIS(2, ref_cloud_safe_global(ref_cloud, 1), "missing");

    global = 1;
    RSS(ref_cloud_item(ref_cloud, global, &item), "find");
    REIS(0, item, "missing");
    RAS(ref_cloud_has_global(ref_cloud, global), "has");

    global = 2;
    RSS(ref_cloud_item(ref_cloud, global, &item), "find");
    REIS(1, item, "missing");
    RAS(ref_cloud_has_global(ref_cloud, global), "has");

    global = 5;
    RAS(!ref_cloud_has_global(ref_cloud, global), "has");

    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store lots */
    REF_GLOB global;
    REF_INT max;
    REF_DBL aux[4] = {0.0, 1.0, 2.0, 3.0};
    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    max = ref_cloud_max(ref_cloud);
    for (global = 0; global <= max; global++) {
      RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    }
    RAS(ref_cloud_max(ref_cloud) > max, "more?");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store key only once with latest value */
    REF_GLOB global = 7;
    REF_DBL aux0[4] = {0.0, 1.0, 2.0, 3.0};
    REF_DBL aux1[4] = {10.0, 11.0, 12.0, 13.0};
    RSS(ref_cloud_create(&ref_cloud, 4), "create");

    RSS(ref_cloud_store(ref_cloud, global, aux0), "store");
    REIS(global, ref_cloud_safe_global(ref_cloud, 0), "set");
    RWDS(aux0[0], ref_cloud_aux(ref_cloud, 0, 0), -1.0, "not set");

    RSS(ref_cloud_store(ref_cloud, global, aux1), "store");
    REIS(global, ref_cloud_safe_global(ref_cloud, 0), "set");
    RWDS(aux1[0], ref_cloud_aux(ref_cloud, 0, 0), -1.0, "not set");

    REIS(1, ref_cloud_n(ref_cloud), "one keys");

    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store two, deep copy */
    REF_CLOUD deep_copy;
    REF_GLOB global;
    REF_DBL aux[4] = {0.0, 1.0, 2.0, 3.0};

    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    global = 5;
    RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    global = 8;
    RSS(ref_cloud_store(ref_cloud, global, aux), "store");

    RSS(ref_cloud_deep_copy(&deep_copy, ref_cloud), "copy");

    global = 8;
    REIS(global, ref_cloud_safe_global(deep_copy, 1), "set");
    RWDS(aux[0], ref_cloud_aux(deep_copy, 0, 1), -1.0, "not set");

    RSS(ref_cloud_free(deep_copy), "free");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  { /* store lots, deep copy */
    REF_CLOUD deep_copy;
    REF_GLOB global;
    REF_INT max;
    REF_DBL aux[4] = {0.0, 1.0, 2.0, 3.0};
    RSS(ref_cloud_create(&ref_cloud, 4), "create");
    max = ref_cloud_max(ref_cloud);
    for (global = 0; global <= max; global++) {
      RSS(ref_cloud_store(ref_cloud, global, aux), "store");
    }
    RSS(ref_cloud_deep_copy(&deep_copy, ref_cloud), "copy");
    RSS(ref_cloud_free(deep_copy), "free");
    RSS(ref_cloud_free(ref_cloud), "free");
  }

  return 0;
}
