
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

#include "ref_layer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--fixture", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) { /* add layers to tet fixture */
    REF_LAYER ref_layer;
    REF_GRID ref_grid;
    REF_INT faceid;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "tet brick");
    RSS(ref_layer_create(&ref_layer, ref_mpi), "create");

    faceid = 6;
    RSS(ref_layer_attach(ref_layer, ref_grid, faceid), "attach");
    faceid = 4;
    RSS(ref_layer_attach(ref_layer, ref_grid, faceid), "attach");
    faceid = 3;
    RSS(ref_layer_attach(ref_layer, ref_grid, faceid), "attach");
    faceid = 5;
    RSS(ref_layer_attach(ref_layer, ref_grid, faceid), "attach");
    RSS(ref_layer_puff(ref_layer, ref_grid), "puff");
    RSS(ref_layer_insert(ref_layer, ref_grid), "insert");
    RSS(ref_layer_recon(ref_layer, ref_grid), "insert");

    RSS(ref_export_by_extension(ref_layer_grid(ref_layer),
                                "ref_layer_test_layer.plt"),
        "tec");
    RSS(ref_export_by_extension(ref_grid, "ref_layer_test_all.plt"), "tec");

    RSS(ref_layer_free(ref_layer), "layer");
    RSS(ref_grid_free(ref_grid), "grid");
    RSS(ref_mpi_stop(), "stop");
    RSS(ref_mpi_free(ref_mpi), "free");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--insert", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) { /* add layers to mesh */
    REF_GRID ref_grid;
    REIS(
        5, argc,
        "expected ref_fixture_test --insert mesh.meshb metric.solb geom.egads");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "import mesh");
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[4]), "load egads");
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]), "part metric");

    RSS(ref_layer_align_quad(ref_grid), "ident");

    RSS(ref_export_by_extension(ref_grid, "ref_layer_test_insert.plt"),
        "export");
    RSS(ref_export_by_extension(ref_grid, "ref_layer_test_insert.meshb"),
        "export");

    RSS(ref_grid_free(ref_grid), "grid");
    RSS(ref_mpi_stop(), "stop");
    RSS(ref_mpi_free(ref_mpi), "free");
    return 0;
  }

  { /* make layers */
    REF_LAYER ref_layer;

    RSS(ref_layer_create(&ref_layer, ref_mpi), "create");

    REIS(0, ref_layer_n(ref_layer), "check total layers");

    RSS(ref_layer_free(ref_layer), "layer");
  }

  RSS(ref_mpi_stop(), "stop");
  RSS(ref_mpi_free(ref_mpi), "free");
  return 0;
}
