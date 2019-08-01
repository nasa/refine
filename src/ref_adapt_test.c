
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_node.h"

#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_dict.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_import.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_adapt.h"
#include "ref_gather.h"

#include "ref_collapse.h"
#include "ref_edge.h"
#include "ref_smooth.h"
#include "ref_split.h"
#include "ref_twod.h"

#include "ref_face.h"
#include "ref_subdiv.h"
#include "ref_validation.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_metric.h"

#include "ref_histogram.h"

#include "ref_cavity.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  { /* adapt twod */
    REF_GRID ref_grid;
    REF_INT i, passes;
    REF_BOOL all_done = REF_FALSE;

    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "set up grid");

    RSS(ref_migrate_to_balance(ref_grid), "balance");

    {
      REF_DBL *metric;
      ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
      RSS(ref_metric_imply_from(metric, ref_grid), "from");
      RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "to");
      RSS(ref_node_ghost_real(ref_grid_node(ref_grid)), "ghost real");
      ref_free(metric);
    }

    passes = 5;
    for (i = 0; i < passes; i++) {
      RSS(ref_adapt_pass(ref_grid, &all_done), "pass");
      RSS(ref_migrate_to_balance(ref_grid), "balance");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
