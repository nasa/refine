
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

#include "ref_blend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  RXS(ref_args_find(argc, argv, "--viz", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REIS(4, argc, "required args: --viz grid.ext geom.egads");
    REIS(1, pos, "required args: --viz grid.ext geom.egads");
    printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    printf("write tec %s\n", "ref_blend_viz.tec");
    RSS(ref_blend_attach(ref_grid), "attach");
    {
      REF_BLEND ref_blend = ref_geom_blend(ref_grid_geom(ref_grid));
      RSS(ref_blend_tec(ref_blend, "ref_blend_viz.tec"), "blend tec");
      RSS(ref_geom_tec(ref_grid, "ref_blend_geom.tec"), "blend tec");
      RSS(ref_export_tec_surf(ref_blend_grid(ref_blend), "ref_blend_surf.tec"),
          "blend tec");
      RSS(ref_export_by_extension(ref_blend_grid(ref_blend),
                                  "ref_blend_surrogate.meshb"),
          "blend export");
    }
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--metric", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REIS(4, argc, "required args: --metric grid.ext geom.egads");
    REIS(1, pos, "required args: --metric grid.ext geom.egads");
    printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
    printf("write tec %s\n", "ref_blend_viz.tec");
    RSS(ref_blend_attach(ref_grid), "attach");
    {
      REF_BLEND ref_blend = ref_geom_blend(ref_grid_geom(ref_grid));
      RSS(ref_blend_tec(ref_blend, "ref_blend_viz.tec"), "blend tec");
      RSS(ref_geom_tec(ref_grid, "ref_blend_geom.tec"), "blend tec");
      RSS(ref_export_tec_surf(ref_blend_grid(ref_blend), "ref_blend_surf.tec"),
          "blend tec");
      RSS(ref_export_by_extension(ref_blend_grid(ref_blend),
                                  "ref_blend_surrogate.meshb"),
          "blend export");
      RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
      RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_blend_curve"), "al");
      RSS(ref_blend_multiscale(ref_blend), "blend multiscale");
      RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_blend_multiscale"),
          "al");
    }
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--import", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REIS(5, argc, "required args: --import grid.ext geom.egads blend.meshb");
    REIS(1, pos, "required args: --import grid.ext geom.egads blend.meshb");
    printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    printf("import blend %s\n", argv[4]);
    RSS(ref_blend_import(ref_grid, argv[4]), "attach");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "blend load");
    {
      REF_BLEND ref_blend = ref_geom_blend(ref_grid_geom(ref_grid));
      RSS(ref_blend_tec(ref_blend, "ref_blend_import_viz.tec"), "blend tec");
      RSS(ref_geom_tec(ref_grid, "ref_blend_import_geom.tec"), "blend tec");
      RSS(ref_export_tec_surf(ref_blend_grid(ref_blend),
                              "ref_blend_import_surf.tec"),
          "blend tec");
    }
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  {
    REF_BLEND ref_blend;
    REF_GRID freeable_ref_grid;
    RSS(ref_grid_create(&freeable_ref_grid, ref_mpi), "create");
    RSS(ref_blend_create(&ref_blend, freeable_ref_grid), "create");
    RSS(ref_blend_free(ref_blend), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
