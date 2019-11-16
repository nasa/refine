
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

#include "ref_adj.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
  {
    REF_ADJ ref_adj;
    REIS(REF_NULL, ref_adj_free(NULL), "dont free NULL");
    RSS(ref_adj_create(&ref_adj), "create");
    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* deep copy empty */
    REF_ADJ original;
    REF_ADJ ref_adj;
    RSS(ref_adj_create(&original), "create");
    RSS(ref_adj_deep_copy(&ref_adj, original), "create");
    RSS(ref_adj_free(original), "free");
    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* add and count */
    REF_ADJ ref_adj;
    REF_INT item;
    RSS(ref_adj_create(&ref_adj), "create");

    RAS(!ref_adj_valid(ref_adj_first(ref_adj, 0)), "empty");

    RSS(ref_adj_add(ref_adj, 0, 12), "add");

    item = ref_adj_first(ref_adj, 0);
    RES(12, ref_adj_safe_ref(ref_adj, item), "added ref");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* remove*/
    REF_ADJ ref_adj;
    REF_INT item;
    RSS(ref_adj_create(&ref_adj), "create");

    RSS(ref_adj_add(ref_adj, 0, 12), "add");
    REIS(REF_INVALID, ref_adj_remove(ref_adj, 0, 13), "remove missing");
    RSS(ref_adj_remove(ref_adj, 0, 12), "remove added");

    item = ref_adj_first(ref_adj, 0);
    RES(REF_EMPTY, ref_adj_safe_ref(ref_adj, item), "added ref");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* iterate */
    REF_ADJ ref_adj;
    REF_INT item, ref, degree;
    RSS(ref_adj_create(&ref_adj), "create");

    degree = 0;
    each_ref_adj_node_item_with_ref(ref_adj, 0, item, ref) degree++;
    RES(0, degree, "empty degree");

    RSS(ref_adj_add(ref_adj, 0, 14), "add");

    degree = 0;
    each_ref_adj_node_item_with_ref(ref_adj, 0, item, ref) {
      degree++;
      RES(14, ref, "check ref");
    }
    RES(1, degree, "node degree");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* empty */
    REF_ADJ ref_adj;
    RSS(ref_adj_create(&ref_adj), "create");

    RAS(ref_adj_empty(ref_adj, 0), "starts empty");
    RSS(ref_adj_add(ref_adj, 0, 14), "add");
    RAS(!ref_adj_empty(ref_adj, 0), "not empty anymore");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* negative node */
    REF_ADJ ref_adj;
    RSS(ref_adj_create(&ref_adj), "create");

    RES(REF_EMPTY, ref_adj_first(ref_adj, -1), "negative first");
    REIS(REF_INVALID, ref_adj_add(ref_adj, -1, 21), "negative add");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* reallocate nodes */
    REF_ADJ ref_adj;
    REF_INT node, item;
    RSS(ref_adj_create(&ref_adj), "create");
    node = ref_adj_nnode(ref_adj);

    RSS(ref_adj_add(ref_adj, node, 15), "add and realloc nodes");

    RAS(ref_adj_nnode(ref_adj) > node, "nodes bigger");

    item = ref_adj_first(ref_adj, node);
    RES(15, ref_adj_safe_ref(ref_adj, item), "realloc has ref");

    item = ref_adj_first(ref_adj, node + 1);
    RES(REF_EMPTY, ref_adj_safe_ref(ref_adj, item), "realloc init empty");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* reallocate adj */
    REF_ADJ ref_adj;
    REF_INT nitem, item;
    RSS(ref_adj_create(&ref_adj), "create");
    nitem = ref_adj_nitem(ref_adj);
    for (item = 0; item < nitem + 1; item++)
      RSS(ref_adj_add(ref_adj, 0, item), "add requiring item realloc");

    RAS(ref_adj_nitem(ref_adj) > nitem, "item bigger");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* add uniquely */
    REF_ADJ ref_adj;
    REF_INT item;
    RSS(ref_adj_create(&ref_adj), "create");

    RAS(!ref_adj_valid(ref_adj_first(ref_adj, 0)), "empty");

    RSS(ref_adj_add_uniquely(ref_adj, 0, 12), "add");
    RSS(ref_adj_add_uniquely(ref_adj, 0, 12), "add");

    item = ref_adj_first(ref_adj, 0);
    REIS(12, ref_adj_safe_ref(ref_adj, item), "added ref");
    REIS(REF_FALSE, ref_adj_valid(ref_adj_item_next(ref_adj, item)), "no next");

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* degree */
    REF_ADJ ref_adj;
    REF_INT degree;
    RSS(ref_adj_create(&ref_adj), "create");

    RSS(ref_adj_degree(ref_adj, 0, &degree), "deg");
    REIS(0, degree, "zeroth degree");

    RSS(ref_adj_add(ref_adj, 0, 12), "add");

    RSS(ref_adj_degree(ref_adj, 0, &degree), "deg");
    REIS(1, degree, "first degree")

    RSS(ref_adj_add(ref_adj, 0, 17), "add");

    RSS(ref_adj_degree(ref_adj, 0, &degree), "deg");
    REIS(2, degree, "second degree")

    RSS(ref_adj_free(ref_adj), "free");
  }

  { /* min degree node */
    REF_ADJ ref_adj;
    REF_INT min_degree, min_degree_node;
    RSS(ref_adj_create(&ref_adj), "create");

    RSS(ref_adj_min_degree_node(ref_adj, &min_degree, &min_degree_node), "deg");
    REIS(REF_EMPTY, min_degree, "empty min degree");
    REIS(REF_EMPTY, min_degree_node, "empty min degree node");

    RSS(ref_adj_add(ref_adj, 0, 12), "add");

    RSS(ref_adj_min_degree_node(ref_adj, &min_degree, &min_degree_node), "deg");
    REIS(1, min_degree, "one min degree");
    REIS(0, min_degree_node, "one min degree node");

    RSS(ref_adj_add(ref_adj, 0, 14), "add");
    RSS(ref_adj_add(ref_adj, 0, 16), "add");
    RSS(ref_adj_add(ref_adj, 5, 27), "add");
    RSS(ref_adj_add(ref_adj, 5, 28), "add");

    RSS(ref_adj_min_degree_node(ref_adj, &min_degree, &min_degree_node), "deg");
    REIS(2, min_degree, "two min degree");
    REIS(5, min_degree_node, "two min degree node");

    RSS(ref_adj_add(ref_adj, 10, 38), "add");
    RSS(ref_adj_add(ref_adj, 10, 32), "add");
    RSS(ref_adj_add(ref_adj, 10, 35), "add");

    RSS(ref_adj_min_degree_node(ref_adj, &min_degree, &min_degree_node), "deg");
    REIS(2, min_degree, "two min degree");
    REIS(5, min_degree_node, "two min degree node");

    RSS(ref_adj_free(ref_adj), "free");
  }

  return 0;
}
