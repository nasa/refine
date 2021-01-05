
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

#include "ref_facelift.h"

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
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_swap.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos = REF_EMPTY;
  REF_BOOL full = REF_FALSE;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  RXS(ref_args_find(argc, argv, "--full", &pos), REF_NOT_FOUND, "arg search");
  full = (pos != REF_EMPTY);

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
    printf("write tec %s\n", "ref_facelift_viz.tec");
    RSS(ref_facelift_attach(ref_grid), "attach");
    {
      REF_FACELIFT ref_facelift = ref_geom_facelift(ref_grid_geom(ref_grid));
      RSS(ref_facelift_tec(ref_facelift, "ref_facelift_viz.tec"),
          "facelift tec");
      RSS(ref_geom_tec(ref_grid, "ref_facelift_geom.tec"), "facelift tec");
      RSS(ref_export_tec_surf(ref_facelift_grid(ref_facelift),
                              "ref_facelift_surf.tec"),
          "facelift tec");
      RSS(ref_export_by_extension(ref_facelift_grid(ref_facelift),
                                  "ref_facelift_surrogate.meshb"),
          "facelift export");
    }
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--surrogate", &pos), REF_NOT_FOUND,
      "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REIS(5, argc,
         "required args: --surrogate orig.ext geom.egads surrogate.ext");
    REIS(1, pos,
         "required args: --surrogate orig.ext geom.egads surrogate.ext");
    if (ref_mpi_once(ref_mpi)) printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    if (ref_mpi_once(ref_mpi)) printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom mark");
    if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
    if (ref_mpi_once(ref_mpi)) printf("underlying egads tolerance\n");
    RSS(ref_geom_verify_topo(ref_grid), "geom topo");
    RSS(ref_geom_verify_param(ref_grid), "geom param");
    ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");
    if (ref_mpi_once(ref_mpi)) printf("underlying egads tolerance\n");
    if (ref_mpi_once(ref_mpi)) printf("ref_facelift_orig_geom.tec\n");
    if (ref_mpi_once(ref_mpi))
      RSS(ref_geom_tec(ref_grid, "ref_facelift_orig_geom.tec"), "geom export");
    printf("facelift surrogate %s\n", argv[4]);
    RSS(ref_facelift_surrogate(ref_grid, argv[4]), "surrogate");
    if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
    if (ref_mpi_once(ref_mpi)) printf("surrogate tolerance\n");
    RSS(ref_geom_verify_topo(ref_grid), "geom topo");
    RSS(ref_geom_verify_param(ref_grid), "geom param");
    ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");
    if (ref_mpi_once(ref_mpi)) printf("ref_facelift_surr_geom.tec\n");
    if (ref_mpi_once(ref_mpi))
      RSS(ref_geom_tec(ref_grid, "ref_facelift_surr_geom.tec"), "geom export");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--metric", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL complexity;
    REIS(5, argc, "required args: --metric grid.ext geom.egads complexity");
    REIS(1, pos, "required args: --metric grid.ext geom.egads complexity");
    printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    complexity = atof(argv[4]);
    printf("complexity %f\n", complexity);
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_facelift_curve"), "al");
    RSS(ref_facelift_multiscale(ref_grid, complexity), "facelift multiscale");
    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_facelift_multiscale"),
        "al");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--import", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REIS(5, argc, "required args: --import grid.ext geom.egads facelift.meshb");
    REIS(1, pos, "required args: --import grid.ext geom.egads facelift.meshb");
    printf("import grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "grid import");
    printf("load geom %s\n", argv[3]);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "geom load");
    printf("import facelift %s\n", argv[4]);
    RSS(ref_facelift_import(ref_grid, argv[4]), "attach");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "facelift load");
    {
      REF_FACELIFT ref_facelift = ref_geom_facelift(ref_grid_geom(ref_grid));
      RSS(ref_facelift_tec(ref_facelift, "ref_facelift_import_viz.tec"),
          "facelift tec");
      RSS(ref_geom_tec(ref_grid, "ref_facelift_import_geom.tec"),
          "facelift tec");
      RSS(ref_export_tec_surf(ref_facelift_grid(ref_facelift),
                              "ref_facelift_import_surf.tec"),
          "facelift tec");
    }
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  {
    REF_FACELIFT ref_facelift;
    REF_GRID freeable_ref_grid;
    RSS(ref_grid_create(&freeable_ref_grid, ref_mpi), "create");
    RSS(ref_facelift_create(&ref_facelift, freeable_ref_grid, REF_FALSE),
        "create");
    RSS(ref_facelift_free(ref_facelift), "free");
  }

  if (ref_egads_allows_construction()) { /* steinmetz P2 surrogate */
    REF_GRID ref_grid, surrogate;
    REF_FACELIFT ref_facelift;
    REF_DBL gap;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "steinmetz"),
        "create geom");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
    RSS(ref_geom_enrich2(surrogate), "enrich2");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    RAB(gap < 1.0e-13, "expected watertight", { printf("gap %e\n", gap); });
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* steinmetz P3 surrogate */
    REF_GRID ref_grid, surrogate;
    REF_FACELIFT ref_facelift;
    REF_DBL gap;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "steinmetz"),
        "create geom");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
    RSS(ref_geom_enrich3(surrogate), "enrich2");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    RAB(gap < 1.0e-13, "expected watertight", { printf("gap %e\n", gap); });
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* sphere (w/ degen) P3 surrogate */
    REF_GRID ref_grid, surrogate;
    REF_FACELIFT ref_facelift;
    REF_DBL gap;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "sphere"), "create geom");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
    RSS(ref_geom_enrich3(surrogate), "enrich2");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    RAB(gap < 1.0e-13, "expected watertight", { printf("gap %e\n", gap); });
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* steinmetz P2 swap */
    REF_GRID ref_grid, surrogate;
    REF_FACELIFT ref_facelift;
    REF_DBL gap;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "steinmetz"),
        "create geom");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
    RSS(ref_geom_enrich2(surrogate), "enrich2");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    RSS(ref_swap_tri_pass(ref_grid), "swap pass");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    RAB(gap < 1.0e-13, "expected watertight", { printf("gap %e\n", gap); });
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (full && ref_egads_allows_construction()) { /* steinmetz P2 smooth */
    REF_GRID ref_grid, surrogate;
    REF_FACELIFT ref_facelift;
    REF_DBL gap;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "steinmetz"),
        "create geom");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
    RSS(ref_geom_enrich2(surrogate), "enrich2");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    RSS(ref_smooth_pass(ref_grid), "smooth pass");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    RAB(gap < 1.0e-12, "expected watertight", { printf("gap %e\n", gap); });
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
