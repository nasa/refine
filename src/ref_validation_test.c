
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

#include "ref_validation.h"

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
#include "ref_fixture.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

REF_STATUS ref_validation_lb8_ugrid_volume(const char *filename) {
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_FILEPOS conn_offset, ibyte = 4;
  REF_DBL *xyz, *xyzs[4];
  REF_INT cell, nodes[4];
  REF_DBL volume, min_volume, max_volume;
  printf("open %s\n", filename);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(1, fread(&nnode, sizeof(REF_INT), 1, file), "nnode");
  REIS(1, fread(&ntri, sizeof(REF_INT), 1, file), "ntri");
  REIS(1, fread(&nqua, sizeof(REF_INT), 1, file), "nqua");
  REIS(1, fread(&ntet, sizeof(REF_INT), 1, file), "ntet");
  REIS(1, fread(&npyr, sizeof(REF_INT), 1, file), "npyr");
  REIS(1, fread(&npri, sizeof(REF_INT), 1, file), "npri");
  REIS(1, fread(&nhex, sizeof(REF_INT), 1, file), "nhex");

  printf("nnode %d\n", nnode);

  ref_malloc(xyz, 3 * nnode, REF_DBL);
  REIS(3 * nnode, fread(xyz, sizeof(REF_DBL), (size_t)(3 * nnode), file),
       "xyz");

  printf("tet %d\n", ntet);
  conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                (REF_FILEPOS)ntri * 4 * ibyte + (REF_FILEPOS)nqua * 5 * ibyte;
  REIS(0, fseeko(file, conn_offset, SEEK_SET), "seek tet failed");
  max_volume = REF_DBL_MIN;
  min_volume = REF_DBL_MAX;
  for (cell = 0; cell < ntet; cell++) {
    REIS(4, fread(nodes, sizeof(REF_INT), 4, file), "tet");
    xyzs[0] = &(xyz[3 * (nodes[0] - 1)]);
    xyzs[1] = &(xyz[3 * (nodes[1] - 1)]);
    xyzs[2] = &(xyz[3 * (nodes[2] - 1)]);
    xyzs[3] = &(xyz[3 * (nodes[3] - 1)]);
    RSS(ref_node_xyz_vol(xyzs, &volume), "vol");
    min_volume = MIN(min_volume, volume);
    max_volume = MAX(max_volume, volume);
  }
  printf("volume range %e %e\n", min_volume, max_volume);
  ref_free(xyz);
  fclose(file);
  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid;
  REF_INT vol_pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--vol", &vol_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != vol_pos) {
    REIS(1, vol_pos, " ref_validation_test --vol grid.lb8.ugrid");
    REIS(3, argc, " ref_validation_test --vol grid.lb8.ugrid");
    RSS(ref_validation_lb8_ugrid_volume(argv[2]), "ugrid vol");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

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
