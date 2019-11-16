
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

#include "ref_node.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  REIS(REF_NULL, ref_node_free(NULL), "dont free NULL");

  { /* init */
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    REIS(0, ref_node_n(ref_node), "init zero nodes");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* deep copy empty */
    REF_NODE original, copy;
    RSS(ref_node_create(&original, ref_mpi), "create");
    RSS(ref_node_deep_copy(&copy, original), "deep copy");

    RSS(ref_node_free(original), "free");
    RSS(ref_node_free(copy), "free");
  }

  {
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RES(REF_EMPTY, ref_node_global(ref_node, 0),
        "global empty for missing node");

    /* first add in order */
    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "first add");
    RES(0, node, "first node is zero");
    RES(1, ref_node_n(ref_node), "count incremented");
    RES(global, ref_node_global(ref_node, 0), "global match for first node");

    global = 20;
    RSS(ref_node_add(ref_node, global, &node), "second add");
    RES(1, node, "second node is one");
    RES(2, ref_node_n(ref_node), "count incremented");
    RES(global, ref_node_global(ref_node, 1), "global match for second node");

    /* removed node invalid */
    REIS(REF_INVALID, ref_node_remove(ref_node, -1), "remove invalid node");
    REIS(REF_INVALID, ref_node_remove(ref_node, 2), "remove invalid node");

    RSS(ref_node_remove(ref_node, 0), "remove first node");
    RES(REF_EMPTY, ref_node_global(ref_node, 0),
        "global empty for removed node");
    RES(1, ref_node_n(ref_node), "count decremented");

    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "replace");
    RES(0, node, "reuse removed node");
    RES(global, ref_node_global(ref_node, node),
        "global match for replaced node");
    RES(2, ref_node_n(ref_node), "count incremented");

    RES(20, ref_node_global(ref_node, 1), "global match for second node");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* remove max node */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    max = ref_node_max(ref_node);
    for (global = 0; global < max; global += 1)
      RSS(ref_node_add(ref_node, global, &node), "realloc");

    RSS(ref_node_remove(ref_node, max - 1), "remove last node");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* remove max node without global */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    max = ref_node_max(ref_node);
    for (global = 0; global < max; global += 1)
      RSS(ref_node_add(ref_node, global, &node), "realloc");

    RSS(ref_node_remove_without_global(ref_node, max - 1), "remove last node");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* add bunch testing realloc */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    max = ref_node_max(ref_node);
    for (global = 10; global < 10 * (max + 2); global += 10)
      RSS(ref_node_add(ref_node, global, &node), "realloc");

    RAS(max < ref_node_max(ref_node), "grow max");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* lookup local from global */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "realloc");

    node = 0;
    REIS(REF_NOT_FOUND, ref_node_local(ref_node, -1, &node), "invalid global");
    RES(REF_EMPTY, node, "expect node empty for invalid global");
    REIS(REF_NOT_FOUND, ref_node_local(ref_node, 5, &node), "invalid global");
    REIS(REF_NOT_FOUND, ref_node_local(ref_node, 200, &node), "invalid global");

    RSS(ref_node_local(ref_node, 10, &node), "return global");
    REIS(0, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* lookup local from global after remove */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "add");
    global = 20;
    RSS(ref_node_add(ref_node, global, &node), "add");

    RSS(ref_node_remove(ref_node, 0), "remove");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(1, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* compact nodes */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT *o2n, *n2o;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RSS(ref_node_add(ref_node, 1, &node), "add");
    RSS(ref_node_add(ref_node, 3, &node), "add");
    RSS(ref_node_add(ref_node, 2, &node), "add");
    RSS(ref_node_remove(ref_node, 1), "remove");

    RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

    REIS(0, o2n[0], "o2n");
    REIS(REF_EMPTY, o2n[1], "o2n");
    REIS(1, o2n[2], "o2n");

    REIS(0, n2o[0], "n2o");
    REIS(2, n2o[1], "n2o");

    ref_free(n2o);
    ref_free(o2n);

    RSS(ref_node_free(ref_node), "free");
  }

  { /* compact local nodes first */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT *o2n, *n2o;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RSS(ref_node_add(ref_node, 1, &node), "add");
    ref_node_part(ref_node, node) = ref_mpi_rank(ref_mpi) + 1;
    RSS(ref_node_add(ref_node, 3, &node), "add");
    RSS(ref_node_add(ref_node, 2, &node), "add");
    RSS(ref_node_remove(ref_node, 1), "remove");

    RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

    REIS(1, o2n[0], "o2n");
    REIS(REF_EMPTY, o2n[1], "o2n");
    REIS(0, o2n[2], "o2n");

    REIS(2, n2o[0], "n2o");
    REIS(0, n2o[1], "n2o");

    ref_free(n2o);
    ref_free(o2n);

    RSS(ref_node_free(ref_node), "free");
  }

  { /* valid */
    REF_INT node;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RAS(!ref_node_valid(ref_node, 0), "empty invalid");
    RSS(ref_node_add(ref_node, 0, &node), "add 0 global");
    RAS(ref_node_valid(ref_node, 0), "zero is valid global");
    RES(0, ref_node_global(ref_node, 0), "zero global");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* unique */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "first");
    global = 20;
    RSS(ref_node_add(ref_node, global, &node), "second");

    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "first");
    REIS(0, node, "return first");
    global = 20;
    RSS(ref_node_add(ref_node, global, &node), "second");
    REIS(1, node, "return second");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* sorted_global rebuild */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 20;
    RSS(ref_node_add(ref_node, global, &node), "realloc");

    global = 10;
    RSS(ref_node_add(ref_node, global, &node), "realloc");

    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "realloc");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(0, node, "wrong local");
    RSS(ref_node_local(ref_node, 10, &node), "return global");
    REIS(1, node, "wrong local");
    RSS(ref_node_local(ref_node, 30, &node), "return global");
    REIS(2, node, "wrong local");

    RSS(ref_node_rebuild_sorted_global(ref_node), "rebuild");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(0, node, "wrong local");
    RSS(ref_node_local(ref_node, 10, &node), "return global");
    REIS(1, node, "wrong local");
    RSS(ref_node_local(ref_node, 30, &node), "return global");
    REIS(2, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* add many to empty */
    REF_INT n = 2, node;
    REF_GLOB global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node, n, global), "many");

    REIS(2, ref_node_n(ref_node), "init zero nodes");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(0, node, "wrong local");
    RSS(ref_node_local(ref_node, 10, &node), "return global");
    REIS(1, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* add many to existing */
    REF_INT n = 2, node;
    REF_GLOB global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RSS(ref_node_add(ref_node, 10, &node), "many");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node, n, global), "many");

    REIS(2, ref_node_n(ref_node), "init zero nodes");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(1, node, "wrong local");
    RSS(ref_node_local(ref_node, 10, &node), "return global");
    REIS(0, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* add many duplicates */
    REF_INT n = 2, node;
    REF_GLOB global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global[0] = 20;
    global[1] = 20;

    RSS(ref_node_add_many(ref_node, n, global), "many");

    REIS(1, ref_node_n(ref_node), "init zero nodes");

    RSS(ref_node_local(ref_node, 20, &node), "return global");
    REIS(0, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* reuse removed global */
    REF_NODE ref_node;
    REF_INT node;
    REF_GLOB global, next;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 3542;
    RSS(ref_node_add(ref_node, global, &node), "add orig");

    RSS(ref_node_remove(ref_node, node), "remove node");

    RSS(ref_node_next_global(ref_node, &next), "next gloabal");
    REIS(global, next, "not reused");

    RSS(ref_node_free(ref_node), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* shift globals */
    REF_INT local, node;
    REF_GLOB global;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 10;
    RSS(ref_node_add(ref_node, global, &local), "add");
    global = 20;
    RSS(ref_node_add(ref_node, global, &local), "add");

    RSS(ref_node_initialize_n_global(ref_node, 30), "init n glob");

    RSS(ref_node_next_global(ref_node, &global), "next");
    REIS(30, global, "expected n global");
    RSS(ref_node_add(ref_node, global, &local), "add");

    RSS(ref_node_shift_new_globals(ref_node), "shift");

    RSS(ref_node_local(ref_node, 30, &node), "return global");
    REIS(2, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* eliminate unused globals */
    REF_INT local, global, node;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 10;
    RSS(ref_node_add(ref_node, global, &local), "add");
    global = 20;
    RSS(ref_node_add(ref_node, global, &local), "add");
    global = 30;
    RSS(ref_node_add(ref_node, global, &local), "add");

    RSS(ref_node_remove(ref_node, 1), "rm");

    RSS(ref_node_initialize_n_global(ref_node, 30), "init n glob");

    RSS(ref_node_eliminate_unused_globals(ref_node), "unused");

    RSS(ref_node_local(ref_node, 29, &node), "return global");
    REIS(2, node, "wrong local");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* ghost int */
    REF_NODE ref_node;
    REF_INT local, ghost, global;
    REF_INT data[2];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &local), "add");
    ref_node_part(ref_node, local) = global;
    data[local] = ref_mpi_rank(ref_mpi);

    global = ref_mpi_rank(ref_mpi) + 1;
    if (global >= ref_mpi_n(ref_mpi)) global = 0;
    if (ref_mpi_para(ref_mpi)) {
      RSS(ref_node_add(ref_node, global, &ghost), "add");
      ref_node_part(ref_node, ghost) = global;
      data[ghost] = REF_EMPTY;
    }

    RSS(ref_node_ghost_int(ref_node, data, 1), "update ghosts");

    global = ref_mpi_rank(ref_mpi);
    REIS(global, data[local], "local changed");
    if (ref_mpi_para(ref_mpi)) {
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global = 0;
      REIS(global, data[ghost], "local changed");
    }
    RSS(ref_node_free(ref_node), "free");
  }

  { /* ghost dbl */
    REF_NODE ref_node;
    REF_INT local, ghost, global;
    REF_INT ldim = 2;
    REF_DBL data[4];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &local), "add");
    ref_node_part(ref_node, local) = global;
    data[0 + ldim * local] = (REF_DBL)ref_mpi_rank(ref_mpi);
    data[1 + ldim * local] = 10.0 * (REF_DBL)ref_mpi_rank(ref_mpi);

    global = ref_mpi_rank(ref_mpi) + 1;
    if (global >= ref_mpi_n(ref_mpi)) global = 0;
    if (ref_mpi_para(ref_mpi)) {
      RSS(ref_node_add(ref_node, global, &ghost), "add");
      ref_node_part(ref_node, ghost) = global;
      data[0 + ldim * ghost] = -1.0;
      data[1 + ldim * ghost] = -1.0;
    }

    RSS(ref_node_ghost_dbl(ref_node, data, ldim), "update ghosts");

    global = ref_mpi_rank(ref_mpi);
    RWDS((REF_DBL)global, data[0 + ldim * local], -1.0, "local changed");
    RWDS(10.0 * (REF_DBL)global, data[1 + ldim * local], -1.0, "local changed");
    if (ref_mpi_para(ref_mpi)) {
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global = 0;
      RWDS((REF_DBL)global, data[0 + ldim * ghost], -1.0, "local changed");
      RWDS(10.0 * (REF_DBL)global, data[1 + ldim * ghost], -1.0,
           "local changed");
    }
    RSS(ref_node_free(ref_node), "free");
  }

  { /* localize ghost int */
    REF_NODE ref_node;
    REF_INT local, ghost, global;
    REF_INT data[2];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &local), "add");
    ref_node_part(ref_node, local) = global;
    data[local] = 1;

    if (ref_mpi_para(ref_mpi)) {
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global = 0;
      RSS(ref_node_add(ref_node, global, &ghost), "add");
      ref_node_part(ref_node, ghost) = global;
      data[ghost] = 1 + global;
    }

    RSS(ref_node_localize_ghost_int(ref_node, data), "update ghosts");

    if (ref_mpi_para(ref_mpi)) {
      global = ref_mpi_rank(ref_mpi);
      REIS(1 + 1 + global, data[local],
           "sum of original (1) and ghost (1+global)");
      REIS(0, data[ghost], "ghost not set to zero changed");
    } else {
      REIS(1, data[local], "local changed");
    }
    RSS(ref_node_free(ref_node), "free");
  }

  { /* add initializes metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    RSS(ref_node_metric_get(ref_node, node, m), "get");
    RWDS(1.0, m[0], -1.0, "m[0]");
    RWDS(0.0, m[1], -1.0, "m[1]");
    RWDS(0.0, m[2], -1.0, "m[2]");
    RWDS(1.0, m[3], -1.0, "m[3]");
    RWDS(0.0, m[4], -1.0, "m[4]");
    RWDS(1.0, m[5], -1.0, "m[5]");

    RSS(ref_node_metric_get_log(ref_node, node, m), "get");
    RWDS(0.0, m[0], -1.0, "m[0]");
    RWDS(0.0, m[1], -1.0, "m[1]");
    RWDS(0.0, m[2], -1.0, "m[2]");
    RWDS(0.0, m[3], -1.0, "m[3]");
    RWDS(0.0, m[4], -1.0, "m[4]");
    RWDS(0.0, m[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* form get metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    RSS(ref_node_metric_form(ref_node, node, 10, 1, 2, 20, 3, 30), "from");
    RSS(ref_node_metric_get(ref_node, node, m), "get");
    RWDS(10.0, m[0], -1.0, "m[0]");
    RWDS(1.0, m[1], -1.0, "m[1]");
    RWDS(2.0, m[2], -1.0, "m[2]");
    RWDS(20.0, m[3], -1.0, "m[3]");
    RWDS(3.0, m[4], -1.0, "m[4]");
    RWDS(30.0, m[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* set get metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m0[6], m1[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    m0[0] = 10;
    m0[1] = 1;
    m0[2] = 2;
    m0[3] = 20;
    m0[4] = 3;
    m0[5] = 30;
    RSS(ref_node_metric_set(ref_node, node, m0), "set");
    RSS(ref_node_metric_get(ref_node, node, m1), "get");
    RWDS(m0[0], m1[0], -1.0, "m[0]");
    RWDS(m0[1], m1[1], -1.0, "m[1]");
    RWDS(m0[2], m1[2], -1.0, "m[2]");
    RWDS(m0[3], m1[3], -1.0, "m[3]");
    RWDS(m0[4], m1[4], -1.0, "m[4]");
    RWDS(m0[5], m1[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* set log get log metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m0[6], m1[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    m0[0] = 2.3; /* exp(10) */
    m0[1] = 1;
    m0[2] = 2;
    m0[3] = 3.0; /* exp(20) */
    m0[4] = 3;
    m0[5] = 3.4; /* exp(30) */
    RSS(ref_node_metric_set_log(ref_node, node, m0), "set");
    RSS(ref_node_metric_get_log(ref_node, node, m1), "get");
    RWDS(m0[0], m1[0], -1.0, "m[0]");
    RWDS(m0[1], m1[1], -1.0, "m[1]");
    RWDS(m0[2], m1[2], -1.0, "m[2]");
    RWDS(m0[3], m1[3], -1.0, "m[3]");
    RWDS(m0[4], m1[4], -1.0, "m[4]");
    RWDS(m0[5], m1[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* set log get metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m0[6], m1[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    m0[0] = 1;
    m0[1] = 0;
    m0[2] = 0;
    m0[3] = 2;
    m0[4] = 0;
    m0[5] = 3;
    RSS(ref_node_metric_set_log(ref_node, node, m0), "set");
    RSS(ref_node_metric_get(ref_node, node, m1), "get");
    RWDS(exp(1), m1[0], -1.0, "m[0]");
    RWDS(0, m1[1], -1.0, "m[1]");
    RWDS(0, m1[2], -1.0, "m[2]");
    RWDS(exp(2), m1[3], -1.0, "m[3]");
    RWDS(0, m1[4], -1.0, "m[4]");
    RWDS(exp(3), m1[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* set get log metric */
    REF_NODE ref_node;
    REF_INT global, node;
    REF_DBL m0[6], m1[6];
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 30;
    RSS(ref_node_add(ref_node, global, &node), "add");
    m0[0] = exp(1);
    m0[1] = 0;
    m0[2] = 0;
    m0[3] = exp(2);
    m0[4] = 0;
    m0[5] = exp(3);
    RSS(ref_node_metric_set(ref_node, node, m0), "set");
    RSS(ref_node_metric_get_log(ref_node, node, m1), "get");
    RWDS(1, m1[0], -1.0, "m[0]");
    RWDS(0, m1[1], -1.0, "m[1]");
    RWDS(0, m1[2], -1.0, "m[2]");
    RWDS(2, m1[3], -1.0, "m[3]");
    RWDS(0, m1[4], -1.0, "m[4]");
    RWDS(3, m1[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* distance in zero metric */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    global = 1;
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 1.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    RSS(ref_node_metric_form(ref_node, node0, 0, 0, 0, 0, 0, 0), "node0 met");
    RSS(ref_node_metric_form(ref_node, node1, 0, 0, 0, 0, 0, 0), "node1 met");
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(0.0, ratio, -1.0, "ratio expected");
  }

  { /* distance in metric */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio, h;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    global = 1;
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 0.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(0.0, ratio, -1.0, "ratio expected");
    RSS(ref_node_ratio_node0(ref_node, node0, node1, &ratio), "ratio0");
    RWDS(0.0, ratio, -1.0, "ratio expected");
    RSS(ref_node_ratio_node0(ref_node, node1, node0, &ratio), "ratio1");
    RWDS(0.0, ratio, -1.0, "ratio expected");

    ref_node_xyz(ref_node, 0, node1) = 1.0;
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(1.0, ratio, -1.0, "ratio expected");
    RSS(ref_node_ratio_node0(ref_node, node0, node1, &ratio), "ratio0");
    RWDS(1.0, ratio, -1.0, "ratio expected");
    RSS(ref_node_ratio_node0(ref_node, node1, node0, &ratio), "ratio1");
    RWDS(1.0, ratio, -1.0, "ratio expected");

    h = 0.5;
    RSS(ref_node_metric_form(ref_node, node0, 1.0 / (h * h), 0, 0, 1, 0, 1),
        "node0 met");
    RSS(ref_node_metric_form(ref_node, node1, 1.0 / (h * h), 0, 0, 1, 0, 1),
        "node1 met");
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(2.0, ratio, -1.0, "ratio expected");

    h = 0.1;
    RSS(ref_node_metric_form(ref_node, node0, 1.0 / (h * h), 0, 0, 1, 0, 1),
        "node0 met");
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(4.970679, ratio, 0.00001, "ratio expected");

    RSS(ref_node_ratio_node0(ref_node, node0, node1, &ratio), "ratio0");
    RWDS(10.0, ratio, -1.0, "ratio expected");
    RSS(ref_node_ratio_node0(ref_node, node1, node0, &ratio), "ratio1");
    RWDS(2.0, ratio, -1.0, "ratio expected");

    RSS(ref_node_free(ref_node), "free");
  }

#define FD_NODE0(xfuncx)                                         \
  {                                                              \
    REF_DBL f, d[3];                                             \
    REF_DBL fd[3], x0, step = 1.0e-7, tol = 1.0e-6;              \
    REF_INT dir;                                                 \
    RSS(xfuncx(ref_node, node0, node1, &f, d), "fd0");           \
    for (dir = 0; dir < 3; dir++) {                              \
      x0 = ref_node_xyz(ref_node, dir, node0);                   \
      ref_node_xyz(ref_node, dir, node0) = x0 + step;            \
      RSS(xfuncx(ref_node, node0, node1, &(fd[dir]), d), "fd+"); \
      fd[dir] = (fd[dir] - f) / step;                            \
      ref_node_xyz(ref_node, dir, node0) = x0;                   \
    }                                                            \
    RSS(xfuncx(ref_node, node0, node1, &f, d), "exact");         \
    RWDS(fd[0], d[0], tol, "dx expected");                       \
    RWDS(fd[1], d[1], tol, "dy expected");                       \
    RWDS(fd[2], d[2], tol, "dz expected");                       \
  }

  { /* derivative of node0 distance in metric */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio;
    REF_DBL f_ratio, d_ratio[3];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    global = 1;
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 0.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    /* same node */
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RSS(ref_node_dratio_dnode0(ref_node, node0, node1, &f_ratio, d_ratio),
        "ratio deriv");
    RWDS(ratio, f_ratio, -1.0, "ratio expected");
    RWDS(0.0, d_ratio[0], -1.0, "dx expected");
    RWDS(0.0, d_ratio[1], -1.0, "dy expected");
    RWDS(0.0, d_ratio[2], -1.0, "dz expected");

    /* length one in x */
    ref_node_xyz(ref_node, 0, node1) = 1.0;

    FD_NODE0(ref_node_dratio_dnode0);
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RSS(ref_node_dratio_dnode0(ref_node, node0, node1, &f_ratio, d_ratio),
        "ratio deriv");
    RWDS(ratio, f_ratio, -1.0, "ratio expected");

    /* length one in xyz */
    ref_node_xyz(ref_node, 0, node1) = 1.0;
    ref_node_xyz(ref_node, 1, node1) = 1.0;
    ref_node_xyz(ref_node, 2, node1) = 1.0;

    FD_NODE0(ref_node_dratio_dnode0);
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RSS(ref_node_dratio_dnode0(ref_node, node0, node1, &f_ratio, d_ratio),
        "ratio deriv");
    RWDS(ratio, f_ratio, -1.0, "ratio expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* derivative of node0 distance in metric  gen */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio, f_ratio, d_ratio[3];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;
    RSS(ref_node_metric_form(ref_node, node0, 1.0, 1.3, 0.4, 1.8, 0.5, 0.5),
        "node0 set");

    global = 1;
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 0.6;
    ref_node_xyz(ref_node, 1, node1) = 0.7;
    ref_node_xyz(ref_node, 2, node1) = 0.8;

    FD_NODE0(ref_node_dratio_dnode0);

    RSS(ref_node_dratio_dnode0(ref_node, node0, node1, &f_ratio, d_ratio),
        "ratio deriv");
    RSS(ref_node_ratio(ref_node, node0, node1, &ratio), "ratio");
    RWDS(ratio, f_ratio, -1.0, "ratio expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* tet volume */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL *xyzs[4];
    REF_DBL vol;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    for (global = 0; global < 4; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0;

    xyzs[0] = ref_node_xyz_ptr(ref_node, nodes[0]);
    xyzs[1] = ref_node_xyz_ptr(ref_node, nodes[1]);
    xyzs[2] = ref_node_xyz_ptr(ref_node, nodes[2]);
    xyzs[3] = ref_node_xyz_ptr(ref_node, nodes[3]);

    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS(1.0 / 6.0, vol, -1.0, "vol expected");
    RSS(ref_node_xyz_vol(xyzs, &vol), "vol");
    RWDS(1.0 / 6.0, vol, -1.0, "vol expected");

    /* inverted tet is negative volume */
    ref_node_xyz(ref_node, 2, nodes[3]) = -1.0;
    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS(-1.0 / 6.0, vol, -1.0, "vol expected");
    RSS(ref_node_xyz_vol(xyzs, &vol), "vol");
    RWDS(-1.0 / 6.0, vol, -1.0, "vol expected");

    RSS(ref_node_free(ref_node), "free");
  }

#define FD_NODES0(xfuncx)                                 \
  {                                                       \
    REF_DBL f, d[3];                                      \
    REF_DBL fd[3], x0, step = 1.0e-8, tol = 5.0e-6;       \
    REF_INT dir;                                          \
    RSS(xfuncx(ref_node, nodes, &f, d), "fd0");           \
    for (dir = 0; dir < 3; dir++) {                       \
      x0 = ref_node_xyz(ref_node, dir, nodes[0]);         \
      ref_node_xyz(ref_node, dir, nodes[0]) = x0 + step;  \
      RSS(xfuncx(ref_node, nodes, &(fd[dir]), d), "fd+"); \
      fd[dir] = (fd[dir] - f) / step;                     \
      ref_node_xyz(ref_node, dir, nodes[0]) = x0;         \
    }                                                     \
    RSS(xfuncx(ref_node, nodes, &f, d), "exact");         \
    RWDS(fd[0], d[0], tol, "dx expected");                \
    RWDS(fd[1], d[1], tol, "dy expected");                \
    RWDS(fd[2], d[2], tol, "dz expected");                \
  }

  { /* tet volume quality deriv */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL f_vol, d_vol[3], vol;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    ref_node_xyz(ref_node, 0, nodes[0]) = 0.1;
    ref_node_xyz(ref_node, 1, nodes[0]) = 0.2;
    ref_node_xyz(ref_node, 2, nodes[0]) = 0.3;

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[1]) = 0.5;
    ref_node_xyz(ref_node, 2, nodes[1]) = 0.7;

    ref_node_xyz(ref_node, 0, nodes[2]) = 0.3;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.4;
    ref_node_xyz(ref_node, 2, nodes[2]) = 0.6;

    ref_node_xyz(ref_node, 0, nodes[3]) = 0.2;
    ref_node_xyz(ref_node, 1, nodes[3]) = 0.7;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.9;

    FD_NODES0(ref_node_tet_dvol_dnode0);

    RSS(ref_node_tet_dvol_dnode0(ref_node, nodes, &f_vol, d_vol),
        "ratio deriv");
    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS(vol, f_vol, -1.0, "vol expected");

    for (global = 0; global < 4; global++) {
      RSS(ref_node_metric_form(ref_node, nodes[global], 20, 0, 0, 20, 0, 20),
          "set 20 iso metric");
    }

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;

    FD_NODES0(ref_node_tet_dquality_dnode0);

    RSS(ref_node_tet_dquality_dnode0(ref_node, nodes, &f_quality, d_quality),
        "deriv");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    RWDS(quality, f_quality, -1.0, "same expected");

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;

    FD_NODES0(ref_node_tet_dquality_dnode0);

    RSS(ref_node_tet_dquality_dnode0(ref_node, nodes, &f_quality, d_quality),
        "deriv");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    RWDS(quality, f_quality, -1.0, "same expected");

    /* test negative tet */
    ref_node_xyz(ref_node, 2, nodes[3]) = -1.9;

    FD_NODES0(ref_node_tet_dquality_dnode0);

    RSS(ref_node_tet_dquality_dnode0(ref_node, nodes, &f_quality, d_quality),
        "deriv");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    RWDS(quality, f_quality, -1.0, "vol expected");

    RSS(ref_node_free(ref_node), "free");
  }

  /* FIXME, break this test up into pieces */

  { /* right tri normal, area, centroid */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL norm[3], centroid[3];
    REF_DBL area;
    REF_BOOL valid;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 2.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 2.0;

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS(2.0, area, -1.0, "expected area");
    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(-2.0, area, -1.0, "expected area");

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS(0.5, area, -1.0, "expected area");
    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(-0.5, area, -1.0, "expected area");
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS(!valid, "expected invalid");

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS(0.0, norm[0], -1.0, "nx");
    RWDS(-1.0, norm[1], -1.0, "ny");
    RWDS(0.0, norm[2], -1.0, "nz");

    RSS(ref_node_tri_centroid(ref_node, nodes, centroid), "norm");
    RWDS(1.0 / 3.0, centroid[0], -1.0, "cx");
    RWDS(0, centroid[1], -1.0, "cy");
    RWDS(1.0 / 3.0, centroid[2], -1.0, "cz");

    global = nodes[2];
    nodes[2] = nodes[1];
    nodes[1] = global;

    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(0.5, area, -1.0, "expected area");
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS(valid, "expected valid");
    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS(0.0, norm[0], -1.0, "nx");
    RWDS(1.0, norm[1], -1.0, "ny");
    RWDS(0.0, norm[2], -1.0, "nz");

    ref_node_xyz(ref_node, 0, nodes[1]) = 0.5;
    ref_node_xyz(ref_node, 2, nodes[1]) = 0.0;

    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(0.0, area, -1.0, "expected area");
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS(!valid, "expected zero area is invalid");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* right tri node angle */
    REF_NODE ref_node;
    REF_INT nodes[3], global, node;
    REF_DBL angle;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 2.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 2.0;

    node = 5;
    REIS(REF_NOT_FOUND, ref_node_tri_node_angle(ref_node, nodes, node, &angle),
         "node should not be in nodes");

    node = 0;
    RSS(ref_node_tri_node_angle(ref_node, nodes, node, &angle), "angle");
    RWDS(0.50 * ref_math_pi, angle, -1.0, "expected 90 degree angle");

    node = 1;
    RSS(ref_node_tri_node_angle(ref_node, nodes, node, &angle), "angle");
    RWDS(0.25 * ref_math_pi, angle, -1.0, "expected 45 degree angle");

    node = 2;
    RSS(ref_node_tri_node_angle(ref_node, nodes, node, &angle), "angle");
    RWDS(0.25 * ref_math_pi, angle, -1.0, "expected 45 degree angle");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* right tri qual epic */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    ref_node->tri_quality = REF_NODE_EPIC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }
    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(0.5 * sqrt(3.0), qual, -1.0, "qual expected");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
      RSS(ref_node_metric_form(ref_node, nodes[global], 100, 0, 0, 100, 0, 100),
          "set 100 iso metric");
    }
    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(0.5 * sqrt(3.0), qual, -1.0, "qual expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* right tri qual jac */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    ref_node->tri_quality = REF_NODE_JAC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }
    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 0, nodes[2]) = 0.5;
    ref_node_xyz(ref_node, 1, nodes[2]) = 0.5 * sqrt(3.0);

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(1.0, qual, -1.0, "qual expected");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }
    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(0.5 * sqrt(3.0), qual, -1.0, "qual expected");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
      RSS(ref_node_metric_form(ref_node, nodes[global], 100, 0, 0, 100, 0, 100),
          "set 100 iso metric");
    }
    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(0.5 * sqrt(3.0), qual, -1.0, "qual expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* tri area quality deriv */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL f_area, d_area[3], area;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      RSS(ref_node_metric_form(ref_node, nodes[global], 10, 0, 0, 10, 0, 10),
          "set 10 iso metric");
    }

    ref_node_xyz(ref_node, 0, nodes[0]) = 0.1;
    ref_node_xyz(ref_node, 1, nodes[0]) = 0.2;
    ref_node_xyz(ref_node, 2, nodes[0]) = 0.3;

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.1;
    ref_node_xyz(ref_node, 1, nodes[1]) = 0.5;
    ref_node_xyz(ref_node, 2, nodes[1]) = 0.7;

    ref_node_xyz(ref_node, 0, nodes[2]) = 0.4;
    ref_node_xyz(ref_node, 1, nodes[2]) = 2.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 0.6;

    FD_NODES0(ref_node_tri_darea_dnode0);
    RSS(ref_node_tri_darea_dnode0(ref_node, nodes, &f_area, d_area),
        "area deriv");
    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS(area, f_area, -1.0, "expected area");

    /* quality */

    FD_NODES0(ref_node_tri_dquality_dnode0);
    RSS(ref_node_tri_dquality_dnode0(ref_node, nodes, &f_quality, d_quality),
        "qual deriv");
    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual deriv");
    RWDS(quality, f_quality, -1.0, "expected quality");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* tri quality jac deriv */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    ref_node->tri_quality = REF_NODE_JAC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      RSS(ref_node_metric_form(ref_node, nodes[global], 14, 25, 40, 45, 71,
                               115),
          "set gen metric");
    }

    ref_node_xyz(ref_node, 0, nodes[0]) = 0.1;
    ref_node_xyz(ref_node, 1, nodes[0]) = 0.2;
    ref_node_xyz(ref_node, 2, nodes[0]) = 0.3;

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.1;
    ref_node_xyz(ref_node, 1, nodes[1]) = 0.5;
    ref_node_xyz(ref_node, 2, nodes[1]) = 0.7;

    ref_node_xyz(ref_node, 0, nodes[2]) = 0.4;
    ref_node_xyz(ref_node, 1, nodes[2]) = 2.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 0.6;

    /* quality */

    FD_NODES0(ref_node_tri_dquality_dnode0);
    RSS(ref_node_tri_dquality_dnode0(ref_node, nodes, &f_quality, d_quality),
        "qual deriv");
    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual deriv");
    RWDS(quality, f_quality, -1.0, "expected quality");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* equilateral tri normal, area, qual */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL norm[3];
    REF_DBL area, qual;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 0, nodes[2]) = 0.5;
    ref_node_xyz(ref_node, 1, nodes[2]) = 0.5 * sqrt(3.0);

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS(0.25 * sqrt(3.0), area, -1.0, "expected area");

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS(1.0, qual, -1.0, "qual expected");

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS(0.0, norm[0], -1.0, "nx");
    RWDS(0.0, norm[1], -1.0, "ny");
    RWDS(0.5 * sqrt(3.0), norm[2], -1.0, "nz");

    global = nodes[2];
    nodes[2] = nodes[1];
    nodes[1] = global;

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS(0.0, norm[0], -1.0, "nx");
    RWDS(0.0, norm[1], -1.0, "ny");
    RWDS(-0.5 * sqrt(3.0), norm[2], -1.0, "nz");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* right tet quality */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    for (global = 0; global < 4; global++) {
      ref_node_xyz(ref_node, 0, global) = 0.0;
      ref_node_xyz(ref_node, 1, global) = 0.0;
      ref_node_xyz(ref_node, 2, global) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0;

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.839947, qual, 0.00001, "epic qual expected");

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.839947, qual, 0.00001, "jac qual expected");

    for (global = 0; global < 4; global++) {
      RSS(ref_node_metric_form(ref_node, nodes[global], 100, 0, 0, 100, 0, 100),
          "set 100 iso metric");
    }

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.839947, qual, 0.00001, "qual expected not metric dep");

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.839947, qual, 0.00001, "jac qual expected");

    /* inverted tet is negative volume */
    ref_node_xyz(ref_node, 2, nodes[3]) = -1.0;
    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(-1.0 / 6.0, qual, -1.0, "qual expected");
    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    /* mapped volume */
    RWDS(-1.0 / 6.0, qual, -1.0, "qual expected");

    RSS(ref_node_free(ref_node), "free");
  }

  if (REF_FALSE) { /* right tet jac quality sign */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    for (global = 0; global < 4; global++) {
      ref_node_xyz(ref_node, 0, global) = 0.0;
      ref_node_xyz(ref_node, 1, global) = 0.0;
      ref_node_xyz(ref_node, 2, global) = 0.0;
      RSS(ref_node_metric_form(ref_node, nodes[global], 5.569680186165702e+00,
                               -1.669546513836191e-01, -5.647583642671333e-02,
                               5.645045675922971e+00, 1.230436787883563e+00,
                               2.881886210688973e+00),
          "set gen metric");
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0;

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.81577123, qual, 0.00001, "epic qual expected");

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(0.839947, qual, 0.00001, "jac qual expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* Regular Tetrahedron vol, quality, ratio */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL *xyzs[4];
    REF_DBL qual, vol, ratio;
    REF_DBL a;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    a = 1.0;

    ref_node_xyz(ref_node, 0, nodes[0]) = 1.0 / 3.0 * sqrt(3.0) * a;
    ref_node_xyz(ref_node, 1, nodes[0]) = 0.0;
    ref_node_xyz(ref_node, 2, nodes[0]) = 0.0;

    ref_node_xyz(ref_node, 0, nodes[1]) = -1.0 / 6.0 * sqrt(3.0) * a;
    ref_node_xyz(ref_node, 1, nodes[1]) = 0.5 * a;
    ref_node_xyz(ref_node, 2, nodes[1]) = 0.0;

    ref_node_xyz(ref_node, 0, nodes[2]) = -1.0 / 6.0 * sqrt(3.0) * a;
    ref_node_xyz(ref_node, 1, nodes[2]) = -0.5 * a;
    ref_node_xyz(ref_node, 2, nodes[2]) = 0.0;

    ref_node_xyz(ref_node, 0, nodes[3]) = 0.0;
    ref_node_xyz(ref_node, 1, nodes[3]) = 0.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0 / 3.0 * sqrt(6.0) * a;

    xyzs[0] = ref_node_xyz_ptr(ref_node, nodes[0]);
    xyzs[1] = ref_node_xyz_ptr(ref_node, nodes[1]);
    xyzs[2] = ref_node_xyz_ptr(ref_node, nodes[2]);
    xyzs[3] = ref_node_xyz_ptr(ref_node, nodes[3]);

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(1.0, qual, -1.0, "qual expected");
    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS(1.0, qual, -1.0, "qual expected");

    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS(1.0 / 12.0 * sqrt(2.0), vol, -1.0, "vol expected");
    RSS(ref_node_xyz_vol(xyzs, &vol), "vol");
    RWDS(1.0 / 12.0 * sqrt(2.0), vol, -1.0, "vol expected");

    RSS(ref_node_ratio(ref_node, nodes[2], nodes[3], &ratio), "ratio");
    RWDS(1.0, ratio, -1.0, "ratio expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* interpolate */
    REF_NODE ref_node;
    REF_INT node0, node1, new_node;
    REF_GLOB global;
    REF_DBL m[6];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 1.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &new_node), "add");
    RSS(ref_node_interpolate_edge(ref_node, node0, node1, 0.5, new_node),
        "interp");

    RWDS(0.5, ref_node_xyz(ref_node, 0, new_node), -1.0, "x");
    RWDS(0.0, ref_node_xyz(ref_node, 1, new_node), -1.0, "y");
    RWDS(0.0, ref_node_xyz(ref_node, 2, new_node), -1.0, "z");

    RSS(ref_node_metric_get(ref_node, new_node, m), "get");
    RWDS(1.0, m[0], -1.0, "m[0]");
    RWDS(0.0, m[1], -1.0, "m[1]");
    RWDS(0.0, m[2], -1.0, "m[2]");
    RWDS(1.0, m[3], -1.0, "m[3]");
    RWDS(0.0, m[4], -1.0, "m[4]");
    RWDS(1.0, m[5], -1.0, "m[5]");

    RSS(ref_node_metric_form(ref_node, node1, 1.0 / (0.1 * 0.1), 0, 0, 1, 0, 1),
        "set 0.1 metric");

    RSS(ref_node_interpolate_edge(ref_node, node0, node1, 0.5, new_node),
        "interp");

    RSS(ref_node_metric_get(ref_node, new_node, m), "get");
    RWDS(1.0 / 0.1, m[0], -1.0, "m[0]");
    RWDS(0.0, m[1], -1.0, "m[1]");
    RWDS(0.0, m[2], -1.0, "m[2]");
    RWDS(1.0, m[3], -1.0, "m[3]");
    RWDS(0.0, m[4], -1.0, "m[4]");
    RWDS(1.0, m[5], -1.0, "m[5]");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* interpolate aux */
    REF_NODE ref_node;
    REF_INT node0, node1, new_node;
    REF_GLOB global;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    ref_node_naux(ref_node) = 2;
    RSS(ref_node_resize_aux(ref_node), "resize aux");

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &node0), "add");
    ref_node_xyz(ref_node, 0, node0) = 0.0;
    ref_node_xyz(ref_node, 1, node0) = 0.0;
    ref_node_xyz(ref_node, 2, node0) = 0.0;

    ref_node_aux(ref_node, 0, node0) = 1.0;
    ref_node_aux(ref_node, 1, node0) = 20.0;

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &node1), "add");
    ref_node_xyz(ref_node, 0, node1) = 1.0;
    ref_node_xyz(ref_node, 1, node1) = 0.0;
    ref_node_xyz(ref_node, 2, node1) = 0.0;

    ref_node_aux(ref_node, 0, node1) = 2.0;
    ref_node_aux(ref_node, 1, node1) = 40.0;

    RSS(ref_node_next_global(ref_node, &global), "next_global");
    RSS(ref_node_add(ref_node, global, &new_node), "add");
    RSS(ref_node_interpolate_edge(ref_node, node0, node1, 0.5, new_node),
        "interp");

    RWDS(1.5, ref_node_aux(ref_node, 0, new_node), -1.0, "a");
    RWDS(30.0, ref_node_aux(ref_node, 1, new_node), -1.0, "a");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* twod tri bary */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL xyz[3];
    REF_DBL bary[3];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[2]) = 1.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(1.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(1.0, bary[2], -1.0, "b2");

    xyz[0] = 0.5;
    xyz[1] = 0.0;
    xyz[2] = 0.5;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.5, bary[1], -1.0, "b1");
    RWDS(0.5, bary[2], -1.0, "b2");

    xyz[0] = 1.0 / 3.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0 / 3.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0 / 3.0, bary[0], -1.0, "b0");
    RWDS(1.0 / 3.0, bary[1], -1.0, "b1");
    RWDS(1.0 / 3.0, bary[2], -1.0, "b2");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* threed tri bary */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL xyz[3];
    REF_DBL bary[3];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(1.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(1.0, bary[2], -1.0, "b2");

    xyz[0] = 0.5;
    xyz[1] = 0.5;
    xyz[2] = 2.0; /* out of plane */

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.5, bary[1], -1.0, "b1");
    RWDS(0.5, bary[2], -1.0, "b2");

    xyz[0] = 1.0 / 3.0;
    xyz[1] = 1.0 / 3.0;
    xyz[2] = -5.0; /* out of plane */

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0 / 3.0, bary[0], -1.0, "b0");
    RWDS(1.0 / 3.0, bary[1], -1.0, "b1");
    RWDS(1.0 / 3.0, bary[2], -1.0, "b2");

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    xyz[2] = -5.0; /* out of plane */

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(-1.0, bary[0], -1.0, "b0");
    RWDS(1.0, bary[1], -1.0, "b1");
    RWDS(1.0, bary[2], -1.0, "b2");

    /* check scale */
    ref_node_xyz(ref_node, 0, nodes[1]) = 2.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 2.0;

    xyz[0] = 2.0;
    xyz[1] = 2.0;
    xyz[2] = 7.0; /* out of plane */

    RSS(ref_node_bary3d(ref_node, nodes, xyz, bary), "bary");
    RWDS(-1.0, bary[0], -1.0, "b0");
    RWDS(1.0, bary[1], -1.0, "b1");
    RWDS(1.0, bary[2], -1.0, "b2");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* threed tet bary */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL xyz[3];
    REF_DBL bary[4];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    for (global = 0; global < 4; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(1.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");
    RWDS(0.0, bary[3], -1.0, "b3");

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(1.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");
    RWDS(0.0, bary[3], -1.0, "b3");

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(1.0, bary[2], -1.0, "b2");
    RWDS(0.0, bary[3], -1.0, "b3");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(0.0, bary[1], -1.0, "b1");
    RWDS(0.0, bary[2], -1.0, "b2");
    RWDS(1.0, bary[3], -1.0, "b3");

    xyz[0] = 1.0 / 3.0;
    xyz[1] = 1.0 / 3.0;
    xyz[2] = 1.0 / 3.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.0, bary[0], -1.0, "b0");
    RWDS(1.0 / 3.0, bary[1], -1.0, "b1");
    RWDS(1.0 / 3.0, bary[2], -1.0, "b2");
    RWDS(1.0 / 3.0, bary[3], -1.0, "b3");

    xyz[0] = 0.25;
    xyz[1] = 0.25;
    xyz[2] = 0.25;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary), "bary");
    RWDS(0.25, bary[0], -1.0, "b0");
    RWDS(0.25, bary[1], -1.0, "b1");
    RWDS(0.25, bary[2], -1.0, "b2");
    RWDS(0.25, bary[3], -1.0, "b3");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* interpolate only bary */
    REF_DBL orig_bary[4] = {1.0, 1.0, -0.5, -0.5};
    REF_DBL clip_bary[4];

    RSS(ref_node_clip_bary4(orig_bary, clip_bary), "clip");

    RWDS(0.5, clip_bary[0], -1.0, "b 0");
    RWDS(0.5, clip_bary[1], -1.0, "b 1");
    RWDS(0.0, clip_bary[2], -1.0, "b 2");
    RWDS(0.0, clip_bary[3], -1.0, "b 3");
  }

  { /* signed tri projection */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL xyz[3];
    REF_DBL projection;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");

    for (global = 0; global < 3; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;

    xyz[0] = 0.2;
    xyz[1] = 0.3;
    xyz[2] = 0.0;

    RSS(ref_node_tri_projection(ref_node, nodes, xyz, &projection), "bary");
    RWDS(0.0, projection, -1.0, "b0");

    xyz[0] = 0.4;
    xyz[1] = 0.5;
    xyz[2] = 1.0;

    RSS(ref_node_tri_projection(ref_node, nodes, xyz, &projection), "bary");
    RWDS(1.0, projection, -1.0, "b0");

    xyz[0] = 0.6;
    xyz[1] = 0.7;
    xyz[2] = -1.0;

    RSS(ref_node_tri_projection(ref_node, nodes, xyz, &projection), "bary");
    RWDS(-1.0, projection, -1.0, "b0");

    /* check scale */
    ref_node_xyz(ref_node, 0, nodes[1]) = 2.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 2.0;

    xyz[0] = 1.0;
    xyz[1] = 3.0;
    xyz[2] = 7.0; /* out of plane */

    RSS(ref_node_tri_projection(ref_node, nodes, xyz, &projection), "bary");
    RWDS(7.0, projection, -1.0, "b0");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* dist to edge */
    REF_NODE ref_node;
    REF_INT nodes[2], global;
    REF_DBL xyz[3];
    REF_DBL dist;

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");

    for (global = 0; global < 2; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;

    xyz[0] = 0.5;
    xyz[1] = 1.0;
    xyz[2] = 0.0;

    RSS(ref_node_dist_to_edge(ref_node, nodes, xyz, &dist), "bary");
    RWDS(1.0, dist, -1.0, "b0");

    /* 3:4:5 traiangle */
    xyz[0] = 2.5;
    xyz[1] = 1.5;
    xyz[2] = 2.0;

    RSS(ref_node_dist_to_edge(ref_node, nodes, xyz, &dist), "bary");
    RWDS(2.5, dist, -1.0, "b0");

    /* 3:4:5 traiangle, shift edge */
    ref_node_xyz(ref_node, 0, nodes[0]) = 5.0;
    ref_node_xyz(ref_node, 0, nodes[1]) = 10.0;
    xyz[0] = 2.5;
    xyz[1] = 1.5;
    xyz[2] = 2.0;

    RSS(ref_node_dist_to_edge(ref_node, nodes, xyz, &dist), "bary");
    RWDS(2.5, dist, -1.0, "b0");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* tet gradient, right tet */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL grad[3], scalar[4];
    REF_DBL *xyzs[4];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    global = 0;
    RSS(ref_node_add(ref_node, global, &(nodes[0])), "add");
    global = 1;
    RSS(ref_node_add(ref_node, global, &(nodes[1])), "add");
    global = 2;
    RSS(ref_node_add(ref_node, global, &(nodes[2])), "add");
    global = 3;
    RSS(ref_node_add(ref_node, global, &(nodes[3])), "add");

    for (global = 0; global < 4; global++) {
      ref_node_xyz(ref_node, 0, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 1, nodes[global]) = 0.0;
      ref_node_xyz(ref_node, 2, nodes[global]) = 0.0;
    }

    ref_node_xyz(ref_node, 0, nodes[1]) = 1.0;
    ref_node_xyz(ref_node, 1, nodes[2]) = 1.0;
    ref_node_xyz(ref_node, 2, nodes[3]) = 1.0;

    xyzs[0] = ref_node_xyz_ptr(ref_node, nodes[0]);
    xyzs[1] = ref_node_xyz_ptr(ref_node, nodes[1]);
    xyzs[2] = ref_node_xyz_ptr(ref_node, nodes[2]);
    xyzs[3] = ref_node_xyz_ptr(ref_node, nodes[3]);

    /* zero gradient */
    scalar[nodes[0]] = 0.0;
    scalar[nodes[1]] = 0.0;
    scalar[nodes[2]] = 0.0;
    scalar[nodes[3]] = 0.0;

    RSS(ref_node_tet_grad_nodes(ref_node, nodes, scalar, grad), "vol");
    RWDS(0.0, grad[0], -1.0, "gradx expected");
    RWDS(0.0, grad[1], -1.0, "grady expected");
    RWDS(0.0, grad[2], -1.0, "gradz expected");
    RSS(ref_node_xyz_grad(xyzs, scalar, grad), "vol");
    RWDS(0.0, grad[0], -1.0, "gradx expected");
    RWDS(0.0, grad[1], -1.0, "grady expected");
    RWDS(0.0, grad[2], -1.0, "gradz expected");

    /* 1-3-5 gradient */
    scalar[nodes[0]] = 0.0;
    scalar[nodes[1]] = 1.0;
    scalar[nodes[2]] = 3.0;
    scalar[nodes[3]] = 5.0;
    RSS(ref_node_tet_grad_nodes(ref_node, nodes, scalar, grad), "vol");
    RWDS(1.0, grad[0], -1.0, "gradx expected");
    RWDS(3.0, grad[1], -1.0, "grady expected");
    RWDS(5.0, grad[2], -1.0, "gradz expected");
    RSS(ref_node_xyz_grad(xyzs, scalar, grad), "vol");
    RWDS(1.0, grad[0], -1.0, "gradx expected");
    RWDS(3.0, grad[1], -1.0, "grady expected");
    RWDS(5.0, grad[2], -1.0, "gradz expected");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* selection */
    REF_NODE ref_node;
    REF_INT node0, node1, center, global, position;
    REF_DBL *elements, median, value;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    center = REF_EMPTY;
    median = 1.0e5;
    if (ref_mpi_once(ref_mpi)) {
      global = 2 * ref_mpi_n(ref_mpi);
      RSS(ref_node_add(ref_node, global, &center), "median");
    }
    global = 0 + 2 * ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &node0), "node0");
    global = 1 + 2 * ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &node1), "node1");
    ref_malloc(elements, ref_node_max(ref_node), REF_DBL);
    elements[node0] = median - (REF_DBL)(10 * ref_mpi_rank(ref_mpi));
    elements[node1] = median + (REF_DBL)(10 * ref_mpi_rank(ref_mpi));
    if (REF_EMPTY != center) {
      elements[center] = median;
    }
    position = ref_mpi_n(ref_mpi) + 1;
    RSS(ref_node_selection(ref_node, elements, position, &value), "median");
    RWDS(median, value, -1.0, "median expected");
    ref_free(elements);
    RSS(ref_node_free(ref_node), "node free");
  }

  { /* push one unused */
    REF_NODE ref_node;
    REF_INT global;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 27;
    RSS(ref_node_push_unused(ref_node, global), "push unused");
    REIS(1, ref_node_n_unused(ref_node), "has one");
    RSS(ref_node_free(ref_node), "node free");
  }

  { /* pop one unused */
    REF_NODE ref_node;
    REF_GLOB global, last;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    global = 27;
    RSS(ref_node_push_unused(ref_node, global), "push unused");
    RSS(ref_node_pop_unused(ref_node, &last), "push unused");
    REIS(0, ref_node_n_unused(ref_node), "has one");
    REIS(global, last, "same one");
    RSS(ref_node_free(ref_node), "node free");
  }

  { /* store lots */
    REF_NODE ref_node;
    REF_INT global;
    REF_INT max;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    max = ref_node_max_unused(ref_node);
    for (global = 0; global <= max; global++) {
      RSS(ref_node_push_unused(ref_node, global), "store");
    }
    RAS(ref_node_max_unused(ref_node) > max, "more?");
    RSS(ref_node_free(ref_node), "free");
  }

  { /* apply offset */
    REF_NODE ref_node;
    REF_GLOB last;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    RSS(ref_node_push_unused(ref_node, 20), "store");
    RSS(ref_node_push_unused(ref_node, 10), "store");
    RSS(ref_node_shift_unused(ref_node, 15, 27), "offset");

    RSS(ref_node_pop_unused(ref_node, &last), "rm");
    REIS(10, last, "has none");

    RSS(ref_node_pop_unused(ref_node, &last), "rm");
    REIS(47, last, "has none");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* sort */
    REF_NODE ref_node;
    REF_GLOB last;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    RSS(ref_node_push_unused(ref_node, 20), "store");
    RSS(ref_node_push_unused(ref_node, 10), "store");
    RSS(ref_node_sort_unused(ref_node), "sort");

    RSS(ref_node_pop_unused(ref_node, &last), "rm");
    REIS(20, last, "has none");
    RSS(ref_node_pop_unused(ref_node, &last), "rm");
    REIS(10, last, "has none");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* erase */
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    RSS(ref_node_push_unused(ref_node, 20), "store");
    RSS(ref_node_push_unused(ref_node, 10), "store");
    RSS(ref_node_erase_unused(ref_node), "sort");

    REIS(0, ref_node_n_unused(ref_node), "has one");

    RSS(ref_node_free(ref_node), "free");
  }

  { /* allgather */
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node, ref_mpi), "create");
    RSS(ref_node_push_unused(ref_node, ref_mpi_rank(ref_mpi)), "store");
    RSS(ref_node_allgather_unused(ref_node), "sort");

    REIS(ref_mpi_n(ref_mpi), ref_node_n_unused(ref_node), "has one");

    RSS(ref_node_free(ref_node), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
