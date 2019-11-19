
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

#include "ref_clump.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_swap.h"

int main(int argc, char *argv[]) {
  if (1 < argc) {
    REF_MPI ref_mpi;
    REF_GRID ref_grid;
    RSS(ref_mpi_create(&ref_mpi), "create");
    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
    if (2 < argc) {
      RSS(ref_part_metric(ref_grid_node(ref_grid), argv[2]), "part m");
    }
    ref_mpi_stopwatch_stop(ref_mpi, "metric");
    RSS(ref_clump_long_edges(ref_grid, 1.0), "long edge");
    RSS(ref_export_tec_surf(ref_grid, "clump_surf.tec"), "surf");
    ref_mpi_stopwatch_stop(ref_mpi, "output");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
  }

  return 0;
}
