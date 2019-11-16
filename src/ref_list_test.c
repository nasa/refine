
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

#include "ref_list.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_mpi.h"

int main(int argc, char *argv[]) {
  REF_LIST ref_list;
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  {
    REIS(REF_NULL, ref_list_free(NULL), "dont free NULL");
    RSS(ref_list_create(&ref_list), "create");
    REIS(0, ref_list_n(ref_list), "init zero");
    RSS(ref_list_free(ref_list), "free");
  }

  { /* deep copy empty */
    REF_LIST original;
    RSS(ref_list_create(&original), "create");
    RSS(ref_list_deep_copy(&ref_list, original), "deep copy");
    RSS(ref_list_free(original), "free");
    RSS(ref_list_free(ref_list), "free");
  }

  { /* store one */
    REF_INT item;
    RSS(ref_list_create(&ref_list), "create");
    item = 27;
    RSS(ref_list_push(ref_list, item), "add");
    REIS(1, ref_list_n(ref_list), "has one");
    RSS(ref_list_free(ref_list), "free");
  }

  { /* remove */
    REF_INT item, last;
    RSS(ref_list_create(&ref_list), "create");

    REIS(REF_FAILURE, ref_list_pop(ref_list, &last), "rm");

    item = 27;
    RSS(ref_list_push(ref_list, item), "add");
    RSS(ref_list_pop(ref_list, &last), "rm");
    REIS(0, ref_list_n(ref_list), "has none");

    RSS(ref_list_free(ref_list), "free");
  }

  { /* store lots */
    REF_INT item, max;
    RSS(ref_list_create(&ref_list), "create");
    max = ref_list_max(ref_list);
    for (item = 0; item <= max; item++) {
      RSS(ref_list_push(ref_list, item), "store");
    }
    RAS(ref_list_max(ref_list) > max, "more?");
    RSS(ref_list_free(ref_list), "free");
  }

  { /* erase */
    RSS(ref_list_create(&ref_list), "create");
    RSS(ref_list_push(ref_list, 20), "store");
    RSS(ref_list_push(ref_list, 10), "store");

    RSS(ref_list_erase(ref_list), "rm -rf");

    REIS(0, ref_list_n(ref_list), "has none");

    RSS(ref_list_free(ref_list), "free");
  }

  { /* contains */
    REF_INT item;
    REF_BOOL contains;
    RSS(ref_list_create(&ref_list), "create");
    item = 27;
    RSS(ref_list_push(ref_list, item), "add");
    RSS(ref_list_contains(ref_list, item, &contains), "have");
    REIS(REF_TRUE, contains, "does have");
    RSS(ref_list_contains(ref_list, 5, &contains), "have");
    REIS(REF_FALSE, contains, "does have");
    RSS(ref_list_free(ref_list), "free");
  }

  { /* delete first */
    REF_INT item;
    RSS(ref_list_create(&ref_list), "create");

    item = 5;
    REIS(REF_NOT_FOUND, ref_list_delete(ref_list, item), "rm");

    item = 21;
    RSS(ref_list_push(ref_list, item), "add");
    item = 22;
    RSS(ref_list_push(ref_list, item), "add");
    item = 23;
    RSS(ref_list_push(ref_list, item), "add");
    REIS(3, ref_list_n(ref_list), "has 3");

    item = 21;
    RSS(ref_list_delete(ref_list, item), "have");
    REIS(2, ref_list_n(ref_list), "has 2");

    item = 23;
    RSS(ref_list_delete(ref_list, item), "have");
    REIS(1, ref_list_n(ref_list), "has 1");

    item = 22;
    RSS(ref_list_delete(ref_list, item), "have");
    REIS(0, ref_list_n(ref_list), "has 0");

    RSS(ref_list_free(ref_list), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
