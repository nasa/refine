
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

#include "ref_comprow.h"

#include <math.h>
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
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_subdiv.h"
#include "ref_twod.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  { /* init */
    REF_GRID ref_grid;
    REF_COMPROW ref_comprow;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    RSS(ref_comprow_create(&ref_comprow, ref_grid), "create");

    RSS(ref_comprow_free(ref_comprow), "comprow");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* tet */
    REF_GRID ref_grid;
    REF_COMPROW ref_comprow;
    REF_INT row, col, entry;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "create");
    RSS(ref_comprow_create(&ref_comprow, ref_grid), "create");

    REIS(16, ref_comprow_nnz(ref_comprow), "nnz");
    REIS(0, ref_comprow->col[3], "diag 0");
    REIS(1, ref_comprow->col[7], "diag 1");
    REIS(2, ref_comprow->col[11], "diag 2");
    REIS(3, ref_comprow->col[15], "diag 3");

    row = 0;
    col = 0;
    RSS(ref_comprow_entry(ref_comprow, row, col, &entry), "ent");
    REIS(3, entry, "entry 0,0");

    row = 0;
    col = 3;
    RSS(ref_comprow_entry(ref_comprow, row, col, &entry), "ent");
    REIS(0, entry, "entry 0,0");

    row = 0;
    col = 2;
    RSS(ref_comprow_entry(ref_comprow, row, col, &entry), "ent");
    REIS(1, entry, "entry 0,0");

    RSS(ref_comprow_free(ref_comprow), "comprow");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* 2 tet */
    REF_GRID ref_grid;
    REF_COMPROW ref_comprow;
    REF_INT row, col, entry;

    RSS(ref_fixture_tet2_grid(&ref_grid, ref_mpi), "create");
    RSS(ref_comprow_create(&ref_comprow, ref_grid), "create");

    REIS(23, ref_comprow_nnz(ref_comprow), "nnz");
    REIS(0, ref_comprow->col[3], "diag 0");
    REIS(1, ref_comprow->col[8], "diag 1");
    REIS(2, ref_comprow->col[13], "diag 2");
    REIS(3, ref_comprow->col[18], "diag 3");
    REIS(4, ref_comprow->col[22], "diag 4");

    row = 0;
    col = 4; /* does not exist */
    REIS(REF_NOT_FOUND, ref_comprow_entry(ref_comprow, row, col, &entry),
         "entry 0,4");
    REIS(REF_EMPTY, entry, "entry 0,4");

    RSS(ref_comprow_free(ref_comprow), "comprow");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
