
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

#include "ref_search.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cavity.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  if (argc == 3 && !ref_mpi_para(ref_mpi)) {
    REF_GRID from, to;
    REF_SEARCH ref_search;

    RSS(ref_mpi_stopwatch_start(ref_mpi), "sw start");
    RSS(ref_part_by_extension(&from, ref_mpi, argv[1]), "import");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "read from grid"), "sw start");
    RSS(ref_part_by_extension(&to, ref_mpi, argv[2]), "import");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "read to grid"), "sw start");

    RSS(ref_search_create(&ref_search, ref_node_n(ref_grid_node(from))),
        "mk search");

    RSS(ref_search_free(ref_search), "search free");
    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* create */
    REF_SEARCH ref_search;
    RSS(ref_search_create(&ref_search, 10), "make search");
    RSS(ref_search_free(ref_search), "search free");
  }

  { /* add parent */
    REF_SEARCH ref_search;
    REF_INT item;
    REF_DBL xyz[3], r;

    RSS(ref_search_create(&ref_search, 10), "make search");

    item = 0;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 0.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    RSS(ref_search_free(ref_search), "search free");
  }

  { /* parent touch */
    REF_SEARCH ref_search;
    REF_LIST ref_list;
    REF_INT item;
    REF_DBL xyz[3], r;

    RSS(ref_search_create(&ref_search, 10), "make search");
    RSS(ref_list_create(&ref_list), "make list");

    item = 10;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    xyz[2] = 1.0;
    r = 0.1;
    RSS(ref_search_touching(ref_search, ref_list, xyz, r), "touches");
    REIS(0, ref_list_n(ref_list), "should skip");

    xyz[0] = 1.0;
    xyz[1] = 1.0;
    xyz[2] = 1.0;
    r = 1.0;
    RSS(ref_search_touching(ref_search, ref_list, xyz, r), "touches");
    REIS(1, ref_list_n(ref_list), "should gather");

    RSS(ref_list_free(ref_list), "list free");
    RSS(ref_search_free(ref_search), "search free");
  }

  { /* parent and child touch */
    REF_SEARCH ref_search;
    REF_LIST ref_list;
    REF_INT item;
    REF_DBL xyz[3], r;

    RSS(ref_search_create(&ref_search, 10), "make search");
    RSS(ref_list_create(&ref_list), "make list");

    item = 10;
    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");
    item = 20;
    xyz[0] = 2.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 4.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.5;
    RSS(ref_search_touching(ref_search, ref_list, xyz, r), "touches");
    REIS(1, ref_list_n(ref_list), "should gather");
    REIS(item, ref_list_value(ref_list, 0), "should item");

    RSS(ref_list_free(ref_list), "list free");
    RSS(ref_search_free(ref_search), "search free");
  }

  { /* nearest candidate */
    REF_SEARCH ref_search;
    REF_LIST ref_list;
    REF_INT item;
    REF_DBL xyz[3], r, trim;

    RSS(ref_search_create(&ref_search, 10), "make search");
    RSS(ref_list_create(&ref_list), "make list");

    item = 28;
    xyz[0] = 2.8;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_search_trim_radius(ref_search, xyz, &trim), "touches");
    RWDS(trim, 3.8, -1, "expected");

    item = 32;
    xyz[0] = 3.2;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 1.0;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_search_trim_radius(ref_search, xyz, &trim), "touches");
    RWDS(trim, 3.8, -1, "expected");

    item = 50;
    xyz[0] = 5.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 0.1;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_search_trim_radius(ref_search, xyz, &trim), "touches");
    RWDS(trim, 3.8, -1, "expected");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_search_nearest_candidates(ref_search, ref_list, xyz), "touches");
    REIS(2, ref_list_n(ref_list), "should gather");
    REIS(28, ref_list_value(ref_list, 0), "should item");
    REIS(32, ref_list_value(ref_list, 1), "should item");

    item = 5;
    xyz[0] = 0.5;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    r = 0.1;
    RSS(ref_search_insert(ref_search, item, xyz, r), "make search");

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_search_trim_radius(ref_search, xyz, &trim), "touches");
    RWDS(trim, 0.6, -1, "expected");

    RSS(ref_list_free(ref_list), "list free");
    RSS(ref_search_free(ref_search), "search free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
