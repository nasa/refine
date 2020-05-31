
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

#include "ref_meshlink.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_args.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_mpi.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--fill", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos == 1 && argc == 4) {
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    printf("association %s\n", argv[2]);
    RSS(ref_meshlink_open(ref_grid, argv[2]), "open");
    printf("block %s\n", argv[3]);
    RSS(ref_meshlink_fill(ref_grid, argv[3]), "open");
    RSS(ref_geom_tec(ref_grid, "ref_meshlink_test.tec"), "geom tec");
    RSS(ref_meshlink_close(ref_grid), "close");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (3 < argc) {
    REF_GRID ref_grid;
    printf("grid source %s\n", argv[1]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "argv import");
    RSS(ref_meshlink_open(ref_grid, argv[2]), "open");
    RSS(ref_meshlink_parse(ref_grid, argv[3]), "open");
    /* RSS(ref_meshlink_cache(ref_grid, argv[3]), "cache"); */
    RSS(ref_geom_tec(ref_grid, "ref_meshlink_test.tec"), "geom tec");
    /* RSS(ref_meshlink_examine(ref_grid, argv[3]), "examine"); */
    RSS(ref_export_by_extension(ref_grid, "ref_meshlink_test.meshb"), "meshb");
    RSS(ref_meshlink_close(ref_grid), "close");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
