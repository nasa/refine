
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

#include "ref_export.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_validation.h"

#include "ref_adj.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_node.h"

#include "ref_cell.h"

#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_face.h"
#include "ref_mpi.h"
#include "ref_sort.h"

#include "ref_fixture.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  if (argc > 1) {
    printf("validating\n");

    printf("reading %s\n", argv[1]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "from ugrid");
    printf("complete.\n");

    RSS(ref_grid_inspect(ref_grid), "inspection");

    printf("validate.\n");
    RSS(ref_validation_volume_status(ref_grid), "tet volume grid");
    RSS(ref_validation_all(ref_grid), "invalid grid");

    printf("vtk.\n");
    RSS(ref_export_vtk(ref_grid, "validate.vtk"), "vtk");

    printf("tec.\n");
    RSS(ref_export_tec(ref_grid, "validate.tec"), "tec");

    RSS(ref_grid_free(ref_grid), "free");
    printf("done.\n");
  }

  if (!ref_mpi_para(ref_mpi)) {
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "twod brick");
    RSS(ref_validation_twod_outward_normal(ref_grid),
        "twod tri outward normal");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
