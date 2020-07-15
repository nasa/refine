
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

#include "ref_fixture.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"
#include "ref_validation.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--brick", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL xmin, xmax, ymin, ymax, zmin, zmax;
    REF_INT xdim, ydim, zdim;
    REIS(12, argc,
         "required args: --brick grid.ext "
         "xmin xmax ymin ymax zmin xmax xdim ydim zdim");
    REIS(1, pos,
         "required args: --brick grid.ext "
         "xmin xmax ymin ymax zmin zmax xdim ydim zdim");
    xmin = atof(argv[3]);
    xmax = atof(argv[4]);
    ymin = atof(argv[5]);
    ymax = atof(argv[6]);
    zmin = atof(argv[7]);
    zmax = atof(argv[8]);

    xdim = atoi(argv[9]);
    ydim = atoi(argv[10]);
    zdim = atoi(argv[11]);

    RSS(ref_fixture_tet_brick_args_grid(&ref_grid, ref_mpi, xmin, xmax, ymin,
                                        ymax, zmin, zmax, xdim, ydim, zdim),
        "args");
    RSS(ref_export_by_extension(ref_grid, argv[2]), "tec");

    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  {
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "fix");

    RSS(ref_validation_cell_node(ref_grid), "invalid pri");

    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;

    RSS(ref_fixture_tri2_grid(&ref_grid, ref_mpi), "fix");

    RSS(ref_validation_cell_node(ref_grid), "invalid pri");

    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;

    RSS(ref_fixture_tet2_grid(&ref_grid, ref_mpi), "fix");

    RSS(ref_validation_cell_node(ref_grid), "invalid pri");

    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_stack_grid(&ref_grid, ref_mpi), "fix");

    RSS(ref_validation_cell_node(ref_grid), "invalid stack");

    RSS(ref_migrate_to_balance(ref_grid), "bal");

    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "fix");

    RSS(ref_validation_cell_node(ref_grid), "invalid pri");

    RSS(ref_migrate_to_balance(ref_grid), "bal");

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (2 == argc) {
    REF_GRID ref_grid;
    RSS(ref_fixture_twod_square_circle(&ref_grid, ref_mpi), "fix");
    RSS(ref_export_by_extension(ref_grid, argv[1]), "export");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
