
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

#include "ref_axi.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_import.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid;

  if (3 != argc) {
    printf("usage: %s input_grid.extension output_grid.extension\n", argv[0]);
    return 0;
  }

  RSS(ref_mpi_create(&ref_mpi), "create");

  printf("importing %s\n", argv[1]);
  RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "from msh");
  printf("complete.\n");

  RSS(ref_grid_inspect(ref_grid), "inspection");

  RSS(ref_axi_wedge(ref_grid), "axi wedge");

  RSS(ref_grid_inspect(ref_grid), "inspection");

  printf("exporting %s\n", argv[2]);
  RSS(ref_export_by_extension(ref_grid, argv[2]), "to file");
  printf("done.\n");

  RSS(ref_grid_free(ref_grid), "free");
  RSS(ref_mpi_free(ref_mpi), "free");
  return 0;
}
