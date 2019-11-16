
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

#include "ref_dict.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  REF_DICT ref_dict;

  {
    REIS(REF_NULL, ref_dict_free(NULL), "dont free NULL");
    RSS(ref_dict_create(&ref_dict), "create");
    REIS(0, ref_dict_n(ref_dict), "init zero");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* missing fails */
    REF_INT key, value;
    key = 2;
    value = 5;
    RSS(ref_dict_create(&ref_dict), "create");
    REIS(REF_NOT_FOUND, ref_dict_value(ref_dict, key, &value), "missing");
    REIS(5, value, "get value");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store one */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");
    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    RSS(ref_dict_value(ref_dict, key, &value), "retrieve");
    REIS(5, value, "get value");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store two */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");
    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    key = 1;
    value = 3;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    key = 1;
    RSS(ref_dict_value(ref_dict, key, &value), "retrieve");
    REIS(3, value, "get value");
    key = 2;
    RSS(ref_dict_value(ref_dict, key, &value), "retrieve");
    REIS(5, value, "get value");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* remove */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");
    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    RSS(ref_dict_remove(ref_dict, key), "remove");
    REIS(0, ref_dict_n(ref_dict), "back to zero");
    REIS(REF_NOT_FOUND, ref_dict_value(ref_dict, key, &value),
         "should not retrieve");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store lots */
    REF_INT key, value, max;
    RSS(ref_dict_create(&ref_dict), "create");
    max = ref_dict_max(ref_dict);
    for (key = 0; key <= max; key++) {
      value = 10 * key;
      RSS(ref_dict_store(ref_dict, key, value), "store");
    }
    RAS(ref_dict_max(ref_dict) > max, "more?");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store key only once with latest value */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");
    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    key = 1;
    value = 3;
    RSS(ref_dict_store(ref_dict, key, value), "store");

    key = 2;
    value = 7;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    key = 2;
    RSS(ref_dict_value(ref_dict, key, &value), "retrieve");
    REIS(7, value, "get value");

    REIS(2, ref_dict_n(ref_dict), "two keys");

    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* has key has value*/
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");

    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");

    REIS(REF_TRUE, ref_dict_has_key(ref_dict, 2), "not found");
    REIS(REF_FALSE, ref_dict_has_key(ref_dict, 1), "found?");

    REIS(REF_TRUE, ref_dict_has_value(ref_dict, 5), "not found");
    REIS(REF_FALSE, ref_dict_has_value(ref_dict, 7), "found?");

    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store two, deep copy */
    REF_DICT deep_copy;
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict), "create");
    key = 2;
    value = 5;
    RSS(ref_dict_store(ref_dict, key, value), "store");
    key = 1;
    value = 3;
    RSS(ref_dict_store(ref_dict, key, value), "store");

    RSS(ref_dict_deep_copy(&deep_copy, ref_dict), "copy");

    key = 1;
    RSS(ref_dict_value(deep_copy, key, &value), "retrieve");
    REIS(3, value, "get value");
    key = 2;
    RSS(ref_dict_value(deep_copy, key, &value), "retrieve");
    REIS(5, value, "get value");
    RSS(ref_dict_free(deep_copy), "free");
    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* store lots, deep copy */
    REF_DICT deep_copy;
    REF_INT key, value, max;
    RSS(ref_dict_create(&ref_dict), "create");
    max = ref_dict_max(ref_dict);
    for (key = 0; key <= max; key++) {
      value = 10 * key;
      RSS(ref_dict_store(ref_dict, key, value), "store");
    }
    RSS(ref_dict_deep_copy(&deep_copy, ref_dict), "copy");
    RSS(ref_dict_free(deep_copy), "free");
    RSS(ref_dict_free(ref_dict), "free");
  }

  return 0;
}
