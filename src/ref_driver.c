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
#include <unistd.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_node.h"

#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_dict.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_import.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_adapt.h"
#include "ref_gather.h"

#include "ref_collapse.h"
#include "ref_edge.h"
#include "ref_smooth.h"
#include "ref_split.h"
#include "ref_twod.h"

#include "ref_face.h"
#include "ref_subdiv.h"
#include "ref_validation.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_metric.h"

#include "ref_histogram.h"

#include "ref_cavity.h"

#include "ref_gather.h"
#include "ref_mpi.h"
#include "ref_part.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define VERSION "not available"
#endif

static void echo_argv(int argc, char *argv[]) {
  int pos;
  printf("\n");
  for (pos = 0; pos < argc; pos++) printf(" %s", argv[pos]);
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid = NULL;
  REF_GRID background_grid = NULL;
  int opt;
  int passes = 15, pass;
  REF_BOOL tecplot_movie = REF_FALSE;
  REF_BOOL sanitize_metric = REF_FALSE;
  REF_BOOL curvature_metric = REF_TRUE;
  REF_BOOL curvature_constraint = REF_FALSE;
  REF_BOOL debug_verbose = REF_FALSE;
  char output_project[1004];
  char output_filename[1024];
  REF_INT ngeom;
  REF_BOOL all_done;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");
  ref_mpi_stopwatch_start(ref_mpi);

  output_project[0] = '\0';

  if (ref_mpi_once(ref_mpi)) {
    printf("version %s, on or after 1.8.20\n", VERSION);
    echo_argv(argc, argv);
  }

  while ((opt = getopt(argc, argv, "i:m:g:r:o:x:s:ltd")) != -1) {
    switch (opt) {
      case 'i':
        if (ref_mpi_para(ref_mpi)) {
          RSS(ref_part_by_extension(&ref_grid, ref_mpi, optarg), "part");
        } else {
          RSS(ref_import_by_extension(&ref_grid, ref_mpi, optarg), "import");
        }
        break;
      case 'g':
        RNS(ref_grid, "input grid must be loaded before geom");
        RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), optarg), "ld e");
        break;
      case 'r':
        RNS(ref_grid, "input grid must be loaded before geom");
        ref_geom_segments_per_radian_of_curvature(ref_grid_geom(ref_grid)) =
            atof(optarg);
        break;
      case 'm':
        RSS(ref_part_metric(ref_grid_node(ref_grid), optarg), "part m");
        curvature_metric = REF_FALSE;
        break;
      case 'o':
        snprintf(output_project, 1004, "%s", optarg);
        break;
      case 'x':
        break;
      case 's':
        passes = atoi(optarg);
        break;
      case 'l':
        sanitize_metric = REF_TRUE;
        break;
      case 't':
        tecplot_movie = REF_TRUE;
        break;
      case 'd':
        debug_verbose = REF_TRUE;
        ref_mpi->debug = REF_TRUE;
        break;
      case '?':
      default:
        printf("parse error -%c\n", optopt);
        printf("usage: \n %s\n", argv[0]);
        printf("       [-i input_grid.ext]\n");
        printf("       [-g geometry.egads]\n");
        printf("       [-r segments_per_curvature_radian]\n");
        printf(
            "       [-m input_project.{solb,sol,metric}] (curvature metric "
            "when missing)\n");
        printf("       [-s number_of_adaptation_sweeps] default is 15\n");
        printf("       [-o output_project]\n");
        printf("       [-x export_grid.ext]\n");
        printf("       [-l] limit metric change\n");
        printf("       [-t] tecplot movie\n");
        printf("       [-d] debug verbose\n");
        return 1;
    }
  }

  RNS(ref_grid, "input grid required");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "read inputs");

  ref_grid_adapt(ref_grid, watch_param) = REF_TRUE;
  ref_grid_adapt(ref_grid, instrument) = REF_TRUE; /* timing datails */
  ref_grid_adapt(ref_grid, collapse_per_pass) = 5; /* timing datails */

  RSS(ref_gather_ngeom(ref_grid_node(ref_grid), ref_grid_geom(ref_grid),
                       REF_GEOM_FACE, &ngeom),
      "count ngeom");
  if (ngeom > 0) {
    curvature_constraint = REF_TRUE;
    RSS(ref_geom_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
    RSS(ref_geom_verify_topo(ref_grid), "geom topo");
    RSS(ref_geom_verify_param(ref_grid), "geom param");
  }

  if (ref_geom_curvature_unlimited(ref_grid_geom(ref_grid))) {
    RAS(!curvature_metric,
        "[-r seg_per_radian] must be grater than 0.1 or -m provided");
    curvature_constraint = REF_FALSE;
  }

  if (curvature_metric) {
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
  } else {
    RSS(ref_grid_deep_copy(&background_grid, ref_grid), "import");
    ref_grid_background(ref_grid) = background_grid;
  }

  RSS(ref_gather_tec_movie_record_button(ref_grid_gather(ref_grid),
                                         tecplot_movie),
      "show time");

  RSS(ref_validation_cell_volume(ref_grid), "vol");
  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");

  if (curvature_constraint) {
    if (ref_mpi_once(ref_mpi)) printf("constrain curvature\n");
    RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "crv const");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
  }
  if (sanitize_metric) {
    if (ref_mpi_once(ref_mpi)) printf("sanitizing metric\n");
    RSS(ref_metric_sanitize(ref_grid), "sant metric");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
  }

  RSS(ref_migrate_to_balance(ref_grid), "balance");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "balance");
  RSS(ref_grid_pack(ref_grid), "pack");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pack");

  for (pass = 0; pass < passes; pass++) {
    if (ref_mpi_once(ref_mpi))
      printf("\n pass %d of %d with %d ranks\n", pass + 1, passes,
             ref_mpi_n(ref_grid_mpi(ref_grid)));
    RSS(ref_adapt_parameter(ref_grid, &all_done), "param");
    if (all_done) break;
    RSS(ref_adapt_pass(ref_grid), "pass");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pass");
    if (curvature_metric) {
      RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "curvature");
    }
    if (NULL != background_grid) {
      RSS(ref_metric_interpolate(ref_grid, background_grid), "interp");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "interp");
    }
    if (curvature_constraint) {
      RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "crv const");
    }
    if (sanitize_metric) {
      RSS(ref_metric_sanitize(ref_grid), "sant metric");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "sant");
    }
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    RSS(ref_migrate_to_balance(ref_grid), "balance");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "balance");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pack");
  }

  RSS(ref_geom_verify_param(ref_grid), "final params");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "verify final params");

  if (strcmp(output_project, "") != 0) {
    if (!ref_grid_twod(ref_grid)) {
      snprintf(output_filename, 1024, "%s.meshb", output_project);
      RSS(ref_gather_by_extension(ref_grid, output_filename), "export");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "gather meshb");
    } else {
      /* single core, not mixed element, 2D */
      if (!ref_mpi_para(ref_mpi) && 0 == ref_cell_n(ref_grid_hex(ref_grid))) {
        snprintf(output_filename, 1024, "%s.meshb", output_project);
        RSS(ref_export_twod_meshb(ref_grid, output_filename), "export");
        ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "gather 2D meshb");
      }
    }
    snprintf(output_filename, 1024, "%s.b8.ugrid", output_project);
    RSS(ref_gather_by_extension(ref_grid, output_filename), "b8.ugrid");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "gather b8.ugrid");

    snprintf(output_filename, 1024, "%s-final-metric.solb", output_project);
    RSS(ref_gather_metric(ref_grid, output_filename), "met met");
    if (!ref_mpi_para(ref_mpi)) {
      snprintf(output_filename, 1024, "%s_surf.tec", output_project);
      RSS(ref_export_tec_surf(ref_grid, output_filename), "surf tec");
      snprintf(output_filename, 1024, "%s_geom.tec", output_project);
      RSS(ref_geom_tec(ref_grid, output_filename), "geom tec");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "tec");
    }
    if (debug_verbose && !ref_grid_twod(ref_grid) && !ref_mpi_para(ref_mpi)) {
      snprintf(output_filename, 1024, "%s_metric_ellipse.tec", output_project);
      RSS(ref_export_tec_metric_ellipse(ref_grid, output_project), "al");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "ellipse");
    }
  }

  /* export via -x grid.ext */
  for (opt = 0; opt < argc - 1; opt++) {
    if (strcmp(argv[opt], "-x") == 0) {
      if (ref_mpi_para(ref_mpi)) {
        RSS(ref_gather_by_extension(ref_grid, argv[opt + 1]), "gather -x");
      } else {
        RSS(ref_export_by_extension(ref_grid, argv[opt + 1]), "export -x");
      }
    }
  }

  if (NULL != background_grid) RSS(ref_grid_free(background_grid), "free");
  if (NULL != ref_grid) RSS(ref_grid_free(ref_grid), "free");

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
