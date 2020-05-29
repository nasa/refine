
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

#include "ref_adapt.h"
#include "ref_args.h"
#include "ref_defs.h"
#include "ref_dist.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_histogram.h"
#include "ref_import.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_meshlink.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_split.h"
#include "ref_validation.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define VERSION "not available"
#endif

static void usage(const char *name) {
  printf("usage: \n %s [--help] <subcommand> [<args>]\n", name);
  printf("\n");
  printf("ref subcommands:\n");
  printf("  adapt        Adapt a mesh\n");
  printf("  bootstrap    Create initial mesh from EGADS file\n");
  printf("  fun3d        Extract a scalar from the primitive solution\n");
  printf("  interpolate  Interpolate a field from one mesh to another\n");
  printf("  loop         Multiscale metric, adapt, and interpolation.\n");
  printf("  multiscale   Compute a multiscale metric.\n");
  printf("  surface      Extract mesh surface.\n");
  printf("  translate    Convert mesh formats.\n");
  printf("\n");
  printf("'ref <command> -h' provides details on a specific subcommand.\n");
}
static void adapt_help(const char *name) {
  printf("usage: \n %s adapt input_mesh.extension [<options>]\n", name);
  printf("  -g  geometry.egads\n");
  printf("  -m  metric.solb (geometry feature metric when missing)\n");
  printf("  -x  output_mesh.extension\n");
  printf("  --partioner selects domain decomposition method.\n");
  printf("      2: ParMETIS graph partioning.\n");
  printf("      3: Zoltan graph partioning.\n");
  printf("      4: Zoltan recursive bisection.\n");
  printf("      5: native recursive bisection.\n");
  printf("\n");
}
static void bootstrap_help(const char *name) {
  printf("usage: \n %s bootstrap project.egads [-t]\n", name);
  printf("  -t  tecplot movie of surface curvature adaptation\n");
  printf("        in files ref_gather_movie.tec and ref_gather_histo.tec\n");
  printf("  --mesher {tetgen|aflr} volume mesher\n");
  printf("  --auto-tparams {or combination of options} adjust .tParams\n");
  printf("        1:missing faces, 2:chord violation, 4:face width (-1:all)\n");
  printf("\n");
}
/*
static void fill_help(const char *name) {
  printf("usage: \n %s fill surface.meshb volume.meshb\n", name);
  printf("\n");
}
*/
static void fun3d_help(const char *name) {
  printf("usage: \n %s fun3d mach project.meshb primitive.solb mach.solb\n",
         name);
  printf(" where primitive.solb is [rho,u,v,w,p] or [rho,u,v,w,p,turb...]\n");
  printf("   in fun3d nondimensionalization\n");
  printf("\n");
}
static void interpolate_help(const char *name) {
  printf(
      "usage: \n %s interpolate donor.meshb donor.solb receptor.meshb "
      "receptor.solb\n",
      name);
  printf("\n");
  printf("  options:\n");
  printf("   --extrude receptor.solb data to two planes.\n");
  printf("   --face <face id> <persist>.solb\n");
  printf("       where persist.solb is copied to receptor.solb\n");
  printf("       and face id is replaced with donor.solb.\n");
  printf("\n");
}

static void vertex_help(const char *name) {
  printf("usage: \n %s vertex input.meshb vertex_index vertex_index ...\n",
         name);
  printf("  vertex_index is zero-based\n");
  printf("\n");
}

static void loop_help(const char *name) {
  printf(
      "usage: \n %s loop <input_project_name> <output_project_name>"
      " complexity [<options>]\n",
      name);
  printf("\n");
  printf("  expects:\n");
  printf(
      "   <input_project_name>.meshb is"
      " mesh with geometry association and model.\n");
  printf(
      "   <input_project_name>_volume.solb is"
      " [rho,u,v,w,p] or [rho,u,v,w,p,turb1]\n");
  printf("    in FUN3D nondimensionalization.\n");
  printf("   complexity is half of the target number of vertices.\n");
  printf("\n");
  printf("  creates:\n");
  printf(
      "   <output_project_name>.meshb is"
      " mesh with geometry association and model.\n");
  printf(
      "   <output_project_name>.lb8.ugrid is"
      " FUN3D compatible little-endian mesh.\n");
  printf(
      "   <output_project_name>-restart.solb is"
      " an interpolated solution.\n");
  printf("\n");
  printf("  options:\n");
  printf("   --norm-power <power> multiscale metric norm power (default 2)\n");
  printf("   --gradation <gradation> (default -1)\n");
  printf("       positive: metric-space gradation stretching ratio.\n");
  printf("       negative: mixed-space gradation.\n");
  printf("   --buffer coarsens the metric approaching the x max boundary.\n");
  printf("   --partioner selects domain decomposition method.\n");
  printf("       2: ParMETIS graph partioning.\n");
  printf("       3: Zoltan graph partioning.\n");
  printf("       4: Zoltan recursive bisection.\n");
  printf("       5: native recursive bisection.\n");
  printf("   --mesh-extension output mesh extension (replaces lb8.ugrid).\n");

  printf("\n");
}
static void multiscale_help(const char *name) {
  printf(
      "usage: \n %s multiscale input_mesh.extension scalar.{solb,snap} "
      "complexity metric.solb\n",
      name);
  printf("   complexity is approximately half the target number of vertices\n");
  printf("\n");
  printf("  options:\n");
  printf("   --norm-power <power> multiscale metric norm power (default 2)\n");
  printf("   --gradation <gradation> (default -1)\n");
  printf("       positive: metric-space gradation stretching ratio.\n");
  printf("       negative: mixed-space gradation.\n");
  printf("   --buffer coarsens the metric approaching the x max boundary.\n");
  printf("   --hessian expects hessian.* in place of scalar.{solb,snap}.\n");
  printf("\n");
}
static void surface_help(const char *name) {
  printf("usage: \n %s surface input_mesh.extension [surface_mesh.tec] \n",
         name);
  printf("\n");
}
static void translate_help(const char *name) {
  printf("usage: \n %s translate input_mesh.extension output_mesh.extension \n",
         name);
  printf("\n");
  printf("  options:\n");
  printf("   --extrude a dim=2 meshb to single layer of prisms.\n");
  printf("   --zero-y-face [face id] explicitly set y=0 on face id.\n");
  printf("\n");
}

static REF_STATUS adapt(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_mesh = NULL;
  char *in_metric = NULL;
  char *in_egads = NULL;
  REF_GRID ref_grid = NULL;
  REF_BOOL curvature_metric = REF_TRUE;
  REF_BOOL all_done = REF_FALSE;
  REF_BOOL all_done0 = REF_FALSE;
  REF_BOOL all_done1 = REF_FALSE;
  REF_INT pass, passes = 30;
  REF_INT opt, pos;
  REF_LONG ntet;

  if (argc < 3) goto shutdown;
  in_mesh = argv[2];

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_mesh);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_mesh), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_mesh);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_mesh), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  RXS(ref_args_char(argc, argv, "-g", &in_egads), REF_NOT_FOUND,
      "egads arg search");
  if (NULL != in_egads) {
    if (ref_mpi_once(ref_mpi)) printf("load egads from %s\n", in_egads);
    RSS(ref_egads_load(ref_grid_geom(ref_grid), in_egads), "load egads");
    ref_mpi_stopwatch_stop(ref_mpi, "load egads");
  } else {
    if (0 < ref_geom_cad_data_size(ref_grid_geom(ref_grid))) {
      if (ref_mpi_once(ref_mpi))
        printf("load egadslite from .meshb byte stream\n");
      RSS(ref_egads_load(ref_grid_geom(ref_grid), NULL), "load egads");
      ref_mpi_stopwatch_stop(ref_mpi, "load egads");
    } else {
      THROW("No geometry available via .meshb or -g option");
    }
  }
  ref_grid_surf(ref_grid) = ref_grid_twod(ref_grid);
  RSS(ref_gather_ncell(ref_grid_node(ref_grid), ref_grid_tet(ref_grid), &ntet),
      "global tets");
  if (0 == ntet) ref_grid_surf(ref_grid) = REF_TRUE;
  RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
  RSS(ref_geom_verify_topo(ref_grid), "geom topo");
  RSS(ref_geom_verify_param(ref_grid), "geom param");
  ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");

  RXS(ref_args_find(argc, argv, "-t", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos)
    RSS(ref_gather_tec_movie_record_button(ref_grid_gather(ref_grid), REF_TRUE),
        "movie on");

  RXS(ref_args_find(argc, argv, "-s", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    passes = atoi(argv[pos + 1]);
    printf("-s %d adaptation passes\n", passes);
  }

  RXS(ref_args_find(argc, argv, "--partioner", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    ref_grid_partitioner(ref_grid) = (REF_MIGRATE_PARTIONER)atoi(argv[pos + 1]);
    printf("--partioner %d partitioner\n", (int)ref_grid_partitioner(ref_grid));
  }

  RXS(ref_args_char(argc, argv, "-m", &in_metric), REF_NOT_FOUND,
      "metric arg search");
  if (NULL != in_metric) {
    if (ref_mpi_once(ref_mpi)) printf("part metric %s\n", in_metric);
    RSS(ref_part_metric(ref_grid_node(ref_grid), in_metric), "part metric");
    curvature_metric = REF_FALSE;
    ref_mpi_stopwatch_stop(ref_mpi, "part metric");
  }

  if (curvature_metric) {
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
  } else {
    RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    ref_mpi_stopwatch_stop(ref_mpi, "crv const");
    RSS(ref_grid_cache_background(ref_grid), "cache");
    ref_mpi_stopwatch_stop(ref_mpi, "cache background metric");
  }

  RSS(ref_validation_cell_volume(ref_grid), "vol");
  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");
  ref_mpi_stopwatch_stop(ref_mpi, "histogram");

  RSS(ref_migrate_to_balance(ref_grid), "balance");
  RSS(ref_grid_pack(ref_grid), "pack");
  ref_mpi_stopwatch_stop(ref_mpi, "pack");

  for (pass = 0; !all_done && pass < passes; pass++) {
    if (ref_mpi_once(ref_mpi))
      printf("\n pass %d of %d with %d ranks\n", pass + 1, passes,
             ref_mpi_n(ref_mpi));
    all_done1 = all_done0;
    RSS(ref_adapt_pass(ref_grid, &all_done0), "pass");
    all_done = all_done0 && all_done1 && (pass > MIN(5, passes));
    ref_mpi_stopwatch_stop(ref_mpi, "pass");
    if (curvature_metric) {
      RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
      ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
    } else {
      RSS(ref_metric_synchronize(ref_grid), "sync with background");
      ref_mpi_stopwatch_stop(ref_mpi, "metric sync");
    }
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    ref_mpi_stopwatch_stop(ref_mpi, "histogram");
    RSS(ref_adapt_tattle_faces(ref_grid), "tattle");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "tattle faces");
    RSS(ref_migrate_to_balance(ref_grid), "balance");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_mpi, "pack");
  }

  RSS(ref_node_implicit_global_from_local(ref_grid_node(ref_grid)),
      "implicit global");
  ref_mpi_stopwatch_stop(ref_mpi, "implicit global");

  RSS(ref_geom_verify_param(ref_grid), "final params");
  ref_mpi_stopwatch_stop(ref_mpi, "verify final params");

  /* export via -x grid.ext and -f final-surf.tec*/
  for (opt = 0; opt < argc - 1; opt++) {
    if (strcmp(argv[opt], "-x") == 0) {
      if (ref_mpi_para(ref_mpi)) {
        if (ref_mpi_once(ref_mpi))
          printf("gather " REF_GLOB_FMT " nodes to %s\n",
                 ref_node_n_global(ref_grid_node(ref_grid)), argv[opt + 1]);
        RSS(ref_gather_by_extension(ref_grid, argv[opt + 1]), "gather -x");
      } else {
        if (ref_mpi_once(ref_mpi))
          printf("export " REF_GLOB_FMT " nodes to %s\n",
                 ref_node_n_global(ref_grid_node(ref_grid)), argv[opt + 1]);
        RSS(ref_export_by_extension(ref_grid, argv[opt + 1]), "export -x");
      }
    }
    if (strcmp(argv[opt], "-f") == 0) {
      if (ref_mpi_once(ref_mpi))
        printf("gather final surface status %s\n", argv[opt + 1]);
      RSS(ref_gather_surf_status_tec(ref_grid, argv[opt + 1]), "gather -f");
    }
  }

  if (NULL != ref_grid) RSS(ref_grid_free(ref_grid), "free");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) adapt_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS bootstrap(REF_MPI ref_mpi, int argc, char *argv[]) {
  size_t end_of_string;
  char project[1000];
  char filename[1024];
  REF_GRID ref_grid = NULL;
  REF_INT t_pos = REF_EMPTY;
  REF_INT s_pos = REF_EMPTY;
  REF_INT mesher_pos = REF_EMPTY;
  REF_INT auto_tparams_pos = REF_EMPTY;
  REF_INT auto_tparams = REF_EGADS_MISSING_TPARAM;
  const char *mesher = "tetgen";
  REF_INT passes = 15;
  REF_INT self_intersections;
  REF_INT global_pos = REF_EMPTY;
  REF_DBL *global_params = NULL;

  if (argc < 3) goto shutdown;
  end_of_string = MIN(1023, strlen(argv[2]));
  if (7 > end_of_string ||
      strncmp(&(argv[2][end_of_string - 6]), ".egads", 6) != 0)
    goto shutdown;
  strncpy(project, argv[2], end_of_string - 6);
  project[end_of_string - 6] = '\0';

  RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
  if (ref_mpi_once(ref_mpi)) {
    printf("loading %s.egads\n", project);
    RSS(ref_egads_out_level(ref_grid_geom(ref_grid), 2), "standard info");
  }
  RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[2]), "ld egads");
  ref_mpi_stopwatch_stop(ref_mpi, "egads load");

  RXS(ref_args_find(argc, argv, "--auto-tparams", &auto_tparams_pos),
      REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != auto_tparams_pos && auto_tparams_pos < argc - 1) {
    auto_tparams = atoi(argv[auto_tparams_pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf("--auto-tparams %d requested\n", auto_tparams);
    if (auto_tparams < 0) {
      auto_tparams = REF_EGADS_ALL_TPARAM;
      if (ref_mpi_once(ref_mpi))
        printf("--auto-tparams %d set to all\n", auto_tparams);
    }
  }

  RXS(ref_args_find(argc, argv, "--global", &global_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != global_pos && global_pos < argc - 3) {
    ref_malloc(global_params, 3, REF_DBL);
    global_params[0] = atof(argv[global_pos + 1]);
    global_params[1] = atof(argv[global_pos + 2]);
    global_params[2] = atof(argv[global_pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf("initial tessellation, global param %f %f %f\n", global_params[0],
             global_params[1], global_params[2]);
  } else {
    if (ref_mpi_once(ref_mpi)) printf("initial tessellation, default param\n");
  }
  RSS(ref_egads_tess(ref_grid, auto_tparams, global_params), "tess egads");
  ref_free(global_params);
  ref_mpi_stopwatch_stop(ref_mpi, "egads tess");
  sprintf(filename, "%s-init-geom.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_geom_tec(ref_grid, filename), "geom export");
  ref_mpi_stopwatch_stop(ref_mpi, "export init-geom");
  sprintf(filename, "%s-init-surf.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_export_tec_surf(ref_grid, filename), "dbg surf");
  ref_mpi_stopwatch_stop(ref_mpi, "export init-surf");
  if (REF_FALSE) {
    sprintf(filename, "%s-init-surf.meshb", project);
    if (ref_mpi_once(ref_mpi))
      RSS(ref_export_by_extension(ref_grid, filename), "dbg meshb");
    ref_mpi_stopwatch_stop(ref_mpi, "export init-surf");
  }
  if (ref_mpi_once(ref_mpi)) printf("verify topo\n");
  RSS(ref_geom_verify_topo(ref_grid), "adapt topo");
  ref_mpi_stopwatch_stop(ref_mpi, "verify topo");
  if (ref_mpi_once(ref_mpi)) printf("verify EGADS param\n");
  RSS(ref_geom_verify_param(ref_grid), "egads params");
  ref_mpi_stopwatch_stop(ref_mpi, "verify param");

  /* slow due to edges, maybe wait for parallel? */
  if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
  RSS(ref_geom_constrain_all(ref_grid), "constrain");
  ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
  if (ref_mpi_once(ref_mpi)) printf("verify constrained param\n");
  RSS(ref_geom_verify_param(ref_grid), "constrained params");
  ref_mpi_stopwatch_stop(ref_mpi, "verify param");

  if (ref_geom_manifold(ref_grid_geom(ref_grid))) {
    if (ref_mpi_once(ref_mpi)) printf("verify manifold\n");
    RSS(ref_validation_boundary_manifold(ref_grid), "manifold");
    ref_mpi_stopwatch_stop(ref_mpi, "tess verification");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("manifold not required for wirebody\n");
  }

  RXS(ref_args_find(argc, argv, "-t", &t_pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != t_pos)
    RSS(ref_gather_tec_movie_record_button(ref_grid_gather(ref_grid), REF_TRUE),
        "movie on");

  RXS(ref_args_find(argc, argv, "--mesher", &mesher_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != mesher_pos && mesher_pos < argc - 1) {
    mesher = argv[mesher_pos + 1];
    printf("--mesher %s requested\n", mesher);
  }

  RXS(ref_args_find(argc, argv, "-s", &s_pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != s_pos && s_pos < argc - 1) {
    passes = atoi(argv[s_pos + 1]);
    printf("-s %d surface adaptation passes\n", passes);
  }

  RSS(ref_adapt_surf_to_geom(ref_grid, passes), "ad");

  RSS(ref_geom_report_tri_area_normdev(ref_grid), "tri status");
  if (ref_mpi_once(ref_mpi)) printf("verify topo\n");
  RSS(ref_geom_verify_topo(ref_grid), "adapt topo");
  if (ref_mpi_once(ref_mpi)) printf("verify param\n");
  RSS(ref_geom_verify_param(ref_grid), "adapt params");
  ref_mpi_stopwatch_stop(ref_mpi, "surf verification");

  ref_grid_partitioner(ref_grid) = REF_MIGRATE_SINGLE;
  RSS(ref_migrate_to_balance(ref_grid), "migrate to single part");
  RSS(ref_grid_pack(ref_grid), "pack");
  ref_mpi_stopwatch_stop(ref_mpi, "pack");

  sprintf(filename, "%s-adapt-surf.meshb", project);
  RSS(ref_gather_by_extension(ref_grid, filename), "gather surf meshb");
  sprintf(filename, "%s-adapt-geom.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_geom_tec(ref_grid, filename), "geom export");
  sprintf(filename, "%s-adapt-surf.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_export_tec_surf(ref_grid, filename), "dbg surf");
  sprintf(filename, "%s-adapt-prop.tec", project);
  RSS(ref_gather_surf_status_tec(ref_grid, filename), "gather surf status");
  ref_mpi_stopwatch_stop(ref_mpi, "export adapt surf");

  RSS(ref_geom_feedback(ref_grid), "feedback");
  ref_mpi_stopwatch_stop(ref_mpi, "geom feedback");

  if (ref_geom_manifold(ref_grid_geom(ref_grid))) {
    if (strncmp(mesher, "t", 1) == 0) {
      if (ref_mpi_once(ref_mpi)) {
        printf("fill volume with TetGen\n");
        RSB(ref_geom_tetgen_volume(ref_grid), "tetgen surface to volume", {
          printf("probing adapted tessellation self-intersections\n");
          RSS(ref_dist_collisions(ref_grid, REF_TRUE, &self_intersections),
              "bumps");
          printf("%d segment-triangle intersections detected.\n",
                 self_intersections);
        });
      }
      ref_mpi_stopwatch_stop(ref_mpi, "tetgen volume");
    } else if (strncmp(mesher, "a", 1) == 0) {
      if (ref_mpi_once(ref_mpi)) {
        printf("fill volume with AFLR3\n");
        RSB(ref_geom_aflr_volume(ref_grid), "aflr surface to volume", {
          printf("probing adapted tessellation self-intersections\n");
          RSS(ref_dist_collisions(ref_grid, REF_TRUE, &self_intersections),
              "bumps");
          printf("%d segment-triangle intersections detected.\n",
                 self_intersections);
        });
      }
      ref_mpi_stopwatch_stop(ref_mpi, "aflr volume");
    } else {
      printf("mesher '%s' not implemented\n", mesher);
      goto shutdown;
    }
    ref_grid_surf(ref_grid) = REF_FALSE; /* needed until vol mesher para */

    RSS(ref_split_edge_geometry(ref_grid), "split geom");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "split geom");
  } else {
    ref_grid_twod(ref_grid) = REF_TRUE; /* assume flat facebody */
  }
  RSS(ref_node_synchronize_globals(ref_grid_node(ref_grid)), "sync glob");

  sprintf(filename, "%s-vol.meshb", project);
  if (ref_mpi_once(ref_mpi))
    printf("gather " REF_GLOB_FMT " nodes to %s\n",
           ref_node_n_global(ref_grid_node(ref_grid)), filename);
  RSS(ref_gather_by_extension(ref_grid, filename), "vol export");
  ref_mpi_stopwatch_stop(ref_mpi, "export volume");

  RSS(ref_validation_cell_volume(ref_grid), "vol");

  RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "curvature");
  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "histogram");

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) bootstrap_help(argv[0]);
  return REF_FAILURE;
}

/*
static REF_STATUS fill(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_file;
  char *in_file;
  REF_GRID ref_grid = NULL;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref fill is not parallel");
  }
  if (argc < 4) goto shutdown;
  in_file = argv[2];
  out_file = argv[3];

  printf("import %s\n", in_file);
  RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "load surface");

  RSS(ref_geom_tetgen_volume(ref_grid), "tetgen surface to volume ");

  printf("export %s\n", out_file);
  RSS(ref_export_by_extension(ref_grid, out_file), "vol export");

  RSS(ref_grid_free(ref_grid), "create");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) fill_help(argv[0]);
  return REF_FAILURE;
}
*/

static REF_STATUS fun3d(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *scalar_name;
  char *out_solb;
  char *in_solb;
  char *in_meshb;
  REF_GRID ref_grid = NULL;
  REF_DBL gamma = 1.4;
  REF_INT ldim, node;
  REF_DBL *solution, *scalar;

  if (argc < 6) goto shutdown;
  scalar_name = argv[2];
  in_meshb = argv[3];
  in_solb = argv[4];
  out_solb = argv[5];

  if (strncmp(scalar_name, "mach", 4) != 0) {
    printf("scalar %s not implemented\n", scalar_name);
    goto shutdown;
  }

  if (ref_mpi_once(ref_mpi)) printf("gamma %f\n", gamma);

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_meshb);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_meshb), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_meshb);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_meshb), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  if (ref_mpi_once(ref_mpi)) printf("part solution %s\n", in_solb);
  RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &solution, in_solb),
      "part solution");
  RAS(5 <= ldim, "expected 5 or more variables per vertex");
  ref_mpi_stopwatch_stop(ref_mpi, "part solution");

  if (ref_mpi_once(ref_mpi)) printf("compute %s\n", scalar_name);
  ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    REF_DBL rho, u, v, w, p, temp;
    rho = solution[0 + ldim * node];
    u = solution[1 + ldim * node];
    v = solution[2 + ldim * node];
    w = solution[3 + ldim * node];
    p = solution[4 + ldim * node];
    temp = gamma * p / rho;
    scalar[node] = sqrt((u * u + v * v + w * w) / temp);
  }
  ref_mpi_stopwatch_stop(ref_mpi, "compute scalar");

  if (ref_mpi_once(ref_mpi))
    printf("writing %s to %s\n", scalar_name, out_solb);
  RSS(ref_gather_scalar_by_extension(ref_grid, 1, scalar, NULL, out_solb),
      "export mach");
  ref_mpi_stopwatch_stop(ref_mpi, "gather scalar");

  ref_free(scalar);
  ref_free(solution);
  RSS(ref_grid_free(ref_grid), "create");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) fun3d_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS ref_grid_extrude_field(REF_GRID twod_grid, REF_INT ldim,
                                         REF_DBL *twod_field,
                                         REF_GRID extruded_grid,
                                         REF_DBL *extruded_field) {
  REF_INT node, local, i;
  REF_GLOB twod_nnode, global;
  twod_nnode = ref_node_n_global(ref_grid_node(twod_grid));
  each_ref_node_valid_node(ref_grid_node(extruded_grid), node) {
    global = ref_node_global(ref_grid_node(extruded_grid), node);
    if (global >= twod_nnode) global -= twod_nnode;
    RSS(ref_node_local(ref_grid_node(twod_grid), global, &local),
        "twod global missing");
    for (i = 0; i < ldim; i++) {
      extruded_field[i + ldim * node] = twod_field[i + ldim * local];
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS interpolate(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *receipt_solb;
  char *receipt_meshb;
  char *donor_solb;
  char *donor_meshb;
  char *persist_solb;
  REF_GRID donor_grid = NULL;
  REF_GRID receipt_grid = NULL;
  REF_INT ldim;
  REF_DBL *donor_solution, *receipt_solution;
  REF_INTERP ref_interp;
  REF_INT pos;
  REF_INT faceid;

  if (argc < 6) goto shutdown;
  donor_meshb = argv[2];
  donor_solb = argv[3];
  receipt_meshb = argv[4];
  receipt_solb = argv[5];

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", donor_meshb);
    RSS(ref_part_by_extension(&donor_grid, ref_mpi, donor_meshb), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "donor part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", donor_meshb);
    RSS(ref_import_by_extension(&donor_grid, ref_mpi, donor_meshb), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "donor import");
  }

  if (ref_mpi_once(ref_mpi)) printf("part solution %s\n", donor_solb);
  RSS(ref_part_scalar(ref_grid_node(donor_grid), &ldim, &donor_solution,
                      donor_solb),
      "part solution");
  ref_mpi_stopwatch_stop(ref_mpi, "donor part solution");

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", donor_meshb);
    RSS(ref_part_by_extension(&receipt_grid, ref_mpi, receipt_meshb), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "receptor part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", donor_meshb);
    RSS(ref_import_by_extension(&receipt_grid, ref_mpi, receipt_meshb),
        "import");
    ref_mpi_stopwatch_stop(ref_mpi, "receptor import");
  }

  if (ref_mpi_once(ref_mpi)) {
    printf("%d leading dim from " REF_GLOB_FMT " donor nodes to " REF_GLOB_FMT
           " receptor nodes\n",
           ldim, ref_node_n_global(ref_grid_node(donor_grid)),
           ref_node_n_global(ref_grid_node(receipt_grid)));
  }

  RXS(ref_args_find(argc, argv, "--face", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 2) {
    faceid = atoi(argv[pos + 1]);
    persist_solb = argv[pos + 2];
    if (ref_mpi_once(ref_mpi))
      printf("part persist solution %s\n", persist_solb);
    RSS(ref_part_scalar(ref_grid_node(receipt_grid), &ldim, &receipt_solution,
                        persist_solb),
        "part solution");
    ref_mpi_stopwatch_stop(ref_mpi, "persist part solution");
    if (ref_mpi_once(ref_mpi)) printf("update solution on faceid %d\n", faceid);
    RSS(ref_interp_create(&ref_interp, donor_grid, receipt_grid),
        "make interp");
    RSS(ref_interp_face_only(ref_interp, faceid, ldim, donor_solution,
                             receipt_solution),
        "map");
    ref_mpi_stopwatch_stop(ref_mpi, "update");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("locate receptor nodes\n");
    RSS(ref_interp_create(&ref_interp, donor_grid, receipt_grid),
        "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    ref_mpi_stopwatch_stop(ref_mpi, "locate");
    if (ref_mpi_once(ref_mpi)) printf("interpolate receptor nodes\n");
    ref_malloc(receipt_solution,
               ldim * ref_node_max(ref_grid_node(receipt_grid)), REF_DBL);
    RSS(ref_interp_scalar(ref_interp, ldim, donor_solution, receipt_solution),
        "interp scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "interp");
  }

  RXS(ref_args_find(argc, argv, "--extrude", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_GRID extruded_grid;
    REF_DBL *extruded_solution = NULL;
    if (ref_mpi_once(ref_mpi)) printf("extrude receptor solution\n");
    RSS(ref_grid_extrude_twod(&extruded_grid, receipt_grid), "extrude");
    ref_malloc(extruded_solution,
               ldim * ref_node_max(ref_grid_node(extruded_grid)), REF_DBL);
    RSS(ref_grid_extrude_field(receipt_grid, ldim, receipt_solution,
                               extruded_grid, extruded_solution),
        "extrude solution");
    if (ref_mpi_once(ref_mpi))
      printf("writing interpolated extruded solution %s\n", receipt_solb);
    RSS(ref_gather_scalar_by_extension(extruded_grid, ldim, extruded_solution,
                                       NULL, receipt_solb),
        "gather recept");
    ref_free(extruded_solution);
    RSS(ref_grid_free(extruded_grid), "free");
  } else {
    if (ref_mpi_once(ref_mpi))
      printf("writing receptor solution %s\n", receipt_solb);
    RSS(ref_gather_scalar_by_extension(receipt_grid, ldim, receipt_solution,
                                       NULL, receipt_solb),
        "gather recept");
    ref_mpi_stopwatch_stop(ref_mpi, "gather receptor");
  }

  ref_free(receipt_solution);
  ref_interp_free(ref_interp);
  RSS(ref_grid_free(receipt_grid), "receipt");
  ref_free(donor_solution);
  RSS(ref_grid_free(donor_grid), "donor");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) interpolate_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS vertex(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_file;
  REF_INT pos, global, local;
  REF_GRID ref_grid = NULL;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref vertex is not parallel");
  }
  if (argc < 4) goto shutdown;
  in_file = argv[2];

  printf("import %s\n", in_file);
  RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "load surface");

  for (pos = 3; pos < argc; pos++) {
    global = atoi(argv[pos]);
    printf("global index %d\n", global);
    RSS(ref_node_local(ref_grid_node(ref_grid), global, &local),
        "global node_index not found");
    RSS(ref_node_location(ref_grid_node(ref_grid), local), "location");
  }

  RSS(ref_grid_free(ref_grid), "create");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) vertex_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS loop(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_project = NULL;
  char *out_project = NULL;
  char filename[1024];
  REF_GRID ref_grid = NULL;
  REF_GRID initial_grid = NULL;
  REF_GRID extruded_grid = NULL;
  REF_BOOL all_done = REF_FALSE;
  REF_BOOL all_done0 = REF_FALSE;
  REF_BOOL all_done1 = REF_FALSE;
  REF_INT pass, passes = 30;
  REF_DBL gamma = 1.4;
  REF_INT ldim, node;
  REF_DBL *initial_field, *ref_field, *extruded_field = NULL, *scalar, *metric;
  REF_INT p = 2;
  REF_DBL gradation = -1.0, complexity;
  REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
  REF_BOOL buffer = REF_FALSE;
  REF_INTERP ref_interp;
  REF_INT pos;
  const char *mach_interpolant = "mach";
  const char *interpolant = mach_interpolant;
  const char *lb8_ugrid = "lb8.ugrid";
  const char *b8_ugrid = "b8.ugrid";
  const char *mesh_extension = lb8_ugrid;

  if (argc < 5) goto shutdown;
  in_project = argv[2];
  out_project = argv[3];
  complexity = atof(argv[4]);

  p = 2;
  RXS(ref_args_find(argc, argv, "--norm-power", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (pos >= argc - 1) {
      printf("option missing value: --norm-power <norm power>\n");
      goto shutdown;
    }
    p = atoi(argv[pos + 1]);
  }

  gradation = -1.0;
  RXS(ref_args_find(argc, argv, "--gradation", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (pos >= argc - 1) {
      printf("option missing value: --gradation <gradation>\n");
      goto shutdown;
    }
    gradation = atof(argv[pos + 1]);
  }

  buffer = REF_FALSE;
  RXS(ref_args_find(argc, argv, "--buffer", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    buffer = REF_TRUE;
  }

  RXS(ref_args_find(argc, argv, "--interpolant", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    interpolant = argv[pos + 1];
  }

  RXS(ref_args_find(argc, argv, "--usm3d", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    mesh_extension = b8_ugrid;
  }

  RXS(ref_args_find(argc, argv, "--mesh-extension", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    mesh_extension = argv[pos + 1];
  }

  if (ref_mpi_once(ref_mpi)) {
    printf("complexity %f\n", complexity);
    printf("Lp=%d\n", p);
    printf("gradation %f\n", gradation);
    printf("reconstruction %d\n", (int)reconstruction);
    printf("buffer %d (zero is inactive)\n", buffer);
    printf("interpolant %s\n", interpolant);
  }

  RXS(ref_args_find(argc, argv, "-s", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    passes = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf("-s %d adaptation passes\n", passes);
  }

  sprintf(filename, "%s.meshb", in_project);
  if (ref_mpi_once(ref_mpi)) printf("part mesh %s\n", filename);
  RSS(ref_part_by_extension(&ref_grid, ref_mpi, filename), "part");
  ref_mpi_stopwatch_stop(ref_mpi, "part");

  RXS(ref_args_find(argc, argv, "--partioner", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    ref_grid_partitioner(ref_grid) = (REF_MIGRATE_PARTIONER)atoi(argv[pos + 1]);
    printf("--partioner %d partitioner\n", (int)ref_grid_partitioner(ref_grid));
  }

  RXS(ref_args_find(argc, argv, "--meshlink", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("meshlink with %s\n", argv[pos + 1]);
    RSS(ref_meshlink_open(ref_grid, argv[pos + 1]), "meshlink init");
    if (ref_mpi_once(ref_mpi)) printf("cache geode orientation\n");
    RSS(ref_meshlink_infer_orientation(ref_grid), "meshlink orient");
  } else {
    RXS(ref_args_find(argc, argv, "--egads", &pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY != pos && pos < argc - 1) {
      if (ref_mpi_once(ref_mpi)) printf("load egads from %s\n", argv[pos + 1]);
      RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[pos + 1]), "load egads");
      ref_mpi_stopwatch_stop(ref_mpi, "load egads");
    } else {
      RAS(0 < ref_geom_cad_data_size(ref_grid_geom(ref_grid)),
          "project.meshb is missing the geometry model record");
      if (ref_mpi_once(ref_mpi))
        printf("load egadslite from .meshb byte stream\n");
      RSS(ref_egads_load(ref_grid_geom(ref_grid), NULL), "load egads");
      ref_mpi_stopwatch_stop(ref_mpi, "load egadslite cad data");
    }
  }

  RSS(ref_grid_deep_copy(&initial_grid, ref_grid), "import");

  RXS(ref_args_find(argc, argv, "--usm3d", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY == pos) {
    sprintf(filename, "%s_volume.solb", in_project);
    if (ref_mpi_once(ref_mpi)) printf("part scalar %s\n", filename);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &initial_field,
                        filename),
        "part scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "part scalar");
  } else {
    sprintf(filename, "%s_volume.plt", in_project);
    if (ref_mpi_once(ref_mpi)) printf("reconstruct scalar %s\n", filename);
    RSS(ref_interp_plt(ref_grid, filename, &ldim, &initial_field),
        "part scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "reconstruct scalar");
  }
  if (ref_mpi_once(ref_mpi)) printf("compute %s\n", interpolant);
  ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  if (strcmp(interpolant, "incomp") == 0) {
    RAS(4 <= ldim,
        "expected 4 or more variables per vertex for incompressible");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL u, v, w, u2;
      u = initial_field[0 + ldim * node];
      v = initial_field[1 + ldim * node];
      w = initial_field[2 + ldim * node];
      /* press = initial_field[3 + ldim * node]; */
      u2 = u * u + v * v + w * w;
      scalar[node] = sqrt(u2);
    }
    ref_mpi_stopwatch_stop(ref_mpi, "compute incompressible scalar");
  } else {
    RAS(5 <= ldim, "expected 5 or more variables per vertex for compressible");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, u, v, w, press, temp, u2, mach2;
      rho = initial_field[0 + ldim * node];
      u = initial_field[1 + ldim * node];
      v = initial_field[2 + ldim * node];
      w = initial_field[3 + ldim * node];
      press = initial_field[4 + ldim * node];
      RAB(ref_math_divisible(press, rho), "can not divide by rho", {
        printf("rho = %e  u = %e  v = %e  w = %e  press = %e\n", rho, u, v, w,
               press);
      });
      temp = gamma * (press / rho);
      u2 = u * u + v * v + w * w;
      RAB(ref_math_divisible(u2, temp), "can not divide by temp", {
        printf("rho = %e  u = %e  v = %e  w = %e  press = %e  temp = %e\n", rho,
               u, v, w, press, temp);
      });
      mach2 = u2 / temp;
      RAB(mach2 >= 0, "negative mach2", {
        printf("rho = %e  u = %e  v = %e  w = %e  press = %e  temp = %e\n", rho,
               u, v, w, press, temp);
      });
      if (strcmp(interpolant, "mach") == 0) {
        scalar[node] = sqrt(mach2);
      } else if (strcmp(interpolant, "htot") == 0) {
        scalar[node] = temp * (1.0 / (gamma - 1.0)) + 0.5 * u2;
      } else if (strcmp(interpolant, "pressure") == 0) {
        scalar[node] = press;
      } else if (strcmp(interpolant, "density") == 0) {
        scalar[node] = rho;
      } else if (strcmp(interpolant, "temperature") == 0) {
        scalar[node] = temp;
      } else {
        RSS(REF_INVALID, "unknown scalar interpolant");
      }
    }
    ref_mpi_stopwatch_stop(ref_mpi, "compute compressible scalar");
  }

  if (ref_mpi_once(ref_mpi)) printf("reconstruct Hessian, compute metric\n");
  ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                    gradation, complexity),
      "lp norm");
  ref_mpi_stopwatch_stop(ref_mpi, "compute metric");

  ref_free(scalar);

  if (buffer) {
    if (ref_mpi_once(ref_mpi)) printf("buffer at complexity %e\n", complexity);
    RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
        "buffer at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "buffer");
  }

  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
  ref_free(metric);

  ref_grid_surf(ref_grid) = ref_grid_twod(ref_grid);
  RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
  RSS(ref_geom_verify_topo(ref_grid), "geom topo");
  RSS(ref_geom_verify_param(ref_grid), "geom param");
  ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");

  RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
  RSS(ref_validation_cell_volume(ref_grid), "vol");
  ref_mpi_stopwatch_stop(ref_mpi, "crv const");
  RSS(ref_grid_cache_background(ref_grid), "cache");
  ref_mpi_stopwatch_stop(ref_mpi, "cache background metric");

  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");
  ref_mpi_stopwatch_stop(ref_mpi, "histogram");

  RSS(ref_migrate_to_balance(ref_grid), "balance");
  RSS(ref_grid_pack(ref_grid), "pack");
  ref_mpi_stopwatch_stop(ref_mpi, "pack");

  for (pass = 0; !all_done && pass < passes; pass++) {
    if (ref_mpi_once(ref_mpi))
      printf("\n pass %d of %d with %d ranks\n", pass + 1, passes,
             ref_mpi_n(ref_mpi));
    all_done1 = all_done0;
    RSS(ref_adapt_pass(ref_grid, &all_done0), "pass");
    all_done = all_done0 && all_done1 && (pass > MIN(5, passes));
    ref_mpi_stopwatch_stop(ref_mpi, "pass");
    RSS(ref_metric_synchronize(ref_grid), "sync with background");
    ref_mpi_stopwatch_stop(ref_mpi, "metric sync");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    ref_mpi_stopwatch_stop(ref_mpi, "histogram");
    RSS(ref_adapt_tattle_faces(ref_grid), "tattle");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "tattle faces");
    RSS(ref_migrate_to_balance(ref_grid), "balance");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_mpi, "pack");
  }

  RSS(ref_node_implicit_global_from_local(ref_grid_node(ref_grid)),
      "implicit global");
  ref_mpi_stopwatch_stop(ref_mpi, "implicit global");

  RSS(ref_geom_verify_param(ref_grid), "final params");
  ref_mpi_stopwatch_stop(ref_mpi, "verify final params");

  sprintf(filename, "%s.meshb", out_project);
  if (ref_mpi_once(ref_mpi))
    printf("gather " REF_GLOB_FMT " nodes to %s\n",
           ref_node_n_global(ref_grid_node(ref_grid)), filename);
  RSS(ref_gather_by_extension(ref_grid, filename), "gather .meshb");
  ref_mpi_stopwatch_stop(ref_mpi, "gather meshb");

  sprintf(filename, "%s.%s", out_project, mesh_extension);
  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("extrude twod\n");
    RSS(ref_grid_extrude_twod(&extruded_grid, ref_grid), "extrude");
    if (ref_mpi_once(ref_mpi))
      printf("gather extruded " REF_GLOB_FMT " nodes to %s\n",
             ref_node_n_global(ref_grid_node(extruded_grid)), filename);
    if (ref_mpi_once(ref_mpi)) printf("gather extruded %s\n", filename);
    RSS(ref_gather_by_extension(extruded_grid, filename),
        "gather mesh extension");
  } else {
    if (ref_mpi_once(ref_mpi))
      printf("gather " REF_GLOB_FMT " nodes to %s\n",
             ref_node_n_global(ref_grid_node(ref_grid)), filename);
    RSS(ref_gather_by_extension(ref_grid, filename), "gather mesh extension");
  }
  ref_mpi_stopwatch_stop(ref_mpi, "gather mesh extension");

  if (ref_mpi_once(ref_mpi)) {
    printf("%d leading dim from " REF_GLOB_FMT " donor nodes to " REF_GLOB_FMT
           " receptor nodes\n",
           ldim, ref_node_n_global(ref_grid_node(initial_grid)),
           ref_node_n_global(ref_grid_node(ref_grid)));
  }

  if (ref_mpi_once(ref_mpi)) printf("locate receptor nodes\n");
  RSS(ref_interp_create(&ref_interp, initial_grid, ref_grid), "make interp");
  RSS(ref_interp_locate(ref_interp), "map");
  ref_mpi_stopwatch_stop(ref_mpi, "locate");

  if (ref_mpi_once(ref_mpi)) printf("interpolate receptor nodes\n");
  ref_malloc(ref_field, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  RSS(ref_interp_scalar(ref_interp, ldim, initial_field, ref_field),
      "interp scalar");
  RSS(ref_interp_free(ref_interp), "free");
  ref_mpi_stopwatch_stop(ref_mpi, "interp");

  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("extruding field of %d\n", ldim);
    ref_malloc(extruded_field,
               ldim * ref_node_max(ref_grid_node(extruded_grid)), REF_DBL);
    RSS(ref_grid_extrude_field(ref_grid, ldim, ref_field, extruded_grid,
                               extruded_field),
        "extrude field");
    RXS(ref_args_find(argc, argv, "--usm3d", &pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY == pos) {
      sprintf(filename, "%s-restart.solb", out_project);
      if (ref_mpi_once(ref_mpi))
        printf("writing interpolated extruded field %s\n", filename);
      RSS(ref_gather_scalar_by_extension(extruded_grid, ldim, extruded_field,
                                         NULL, filename),
          "gather recept");
    } else {
      sprintf(filename, "%s.solb", out_project);
      if (ref_mpi_once(ref_mpi))
        printf("writing interpolated field at prism cell centers %s\n",
               filename);
      RSS(ref_gather_scalar_cell_solb(extruded_grid, ldim, extruded_field,
                                      filename),
          "gather cell center");
    }
    ref_free(extruded_field);
    ref_grid_free(extruded_grid);
  } else {
    RXS(ref_args_find(argc, argv, "--usm3d", &pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY == pos) {
      sprintf(filename, "%s-restart.solb", out_project);
      if (ref_mpi_once(ref_mpi))
        printf("writing interpolated field %s\n", filename);
      RSS(ref_gather_scalar_by_extension(ref_grid, ldim, ref_field, NULL,
                                         filename),
          "gather recept");
    } else {
      sprintf(filename, "%s.solb", out_project);
      if (ref_mpi_once(ref_mpi))
        printf("writing interpolated field at tet cell centers %s\n", filename);
      RSS(ref_gather_scalar_cell_solb(ref_grid, ldim, ref_field, filename),
          "gather cell center");
    }
  }
  ref_mpi_stopwatch_stop(ref_mpi, "gather receptor");

  ref_free(ref_field);
  ref_free(initial_field);
  RSS(ref_grid_free(initial_grid), "free");

  /* export via -x grid.ext and -f final-surf.tec*/
  for (pos = 0; pos < argc - 1; pos++) {
    if (strcmp(argv[pos], "-x") == 0) {
      if (ref_mpi_para(ref_mpi)) {
        if (ref_mpi_once(ref_mpi))
          printf("gather " REF_GLOB_FMT " nodes to %s\n",
                 ref_node_n_global(ref_grid_node(ref_grid)), argv[pos + 1]);
        RSS(ref_gather_by_extension(ref_grid, argv[pos + 1]), "gather -x");
      } else {
        if (ref_mpi_once(ref_mpi))
          printf("export " REF_GLOB_FMT " nodes to %s\n",
                 ref_node_n_global(ref_grid_node(ref_grid)), argv[pos + 1]);
        RSS(ref_export_by_extension(ref_grid, argv[pos + 1]), "export -x");
      }
    }
    if (strcmp(argv[pos], "-f") == 0) {
      if (ref_mpi_once(ref_mpi))
        printf("gather final surface status %s\n", argv[pos + 1]);
      RSS(ref_gather_surf_status_tec(ref_grid, argv[pos + 1]), "gather -f");
    }
  }

  RSS(ref_grid_free(ref_grid), "free");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) loop_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS multiscale(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_metric;
  char *in_mesh;
  char *in_scalar;
  REF_GRID ref_grid = NULL;
  REF_INT ldim;
  REF_DBL *scalar = NULL;
  REF_DBL *metric = NULL;
  REF_INT p;
  REF_DBL gradation, complexity, current_complexity;
  REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
  REF_INT pos;
  REF_BOOL buffer;

  if (argc < 6) goto shutdown;
  in_mesh = argv[2];
  in_scalar = argv[3];
  complexity = atof(argv[4]);
  out_metric = argv[5];

  p = 2;
  RXS(ref_args_find(argc, argv, "--norm-power", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (pos >= argc - 1) {
      printf("option missing value: --norm-power <norm power>\n");
      goto shutdown;
    }
    p = atoi(argv[pos + 1]);
  }

  gradation = -1.0;
  RXS(ref_args_find(argc, argv, "--gradation", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (pos >= argc - 1) {
      printf("option missing value: --gradation <gradation>\n");
      goto shutdown;
    }
    gradation = atof(argv[pos + 1]);
  }

  buffer = REF_FALSE;
  RXS(ref_args_find(argc, argv, "--buffer", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    buffer = REF_TRUE;
  }

  if (ref_mpi_once(ref_mpi)) {
    printf("complexity %f\n", complexity);
    printf("Lp=%d\n", p);
    printf("gradation %f\n", gradation);
    printf("reconstruction %d\n", (int)reconstruction);
    printf("buffer %d (zero is inactive)\n", buffer);
  }

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_mesh);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_mesh), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_mesh);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_mesh), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  RXS(ref_args_find(argc, argv, "--hessian", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (ref_mpi_once(ref_mpi)) printf("part hessian %s\n", in_scalar);
    RSS(ref_part_metric(ref_grid_node(ref_grid), in_scalar), "part scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "part metric");
    RSS(ref_metric_from_node(metric, ref_grid_node(ref_grid)), "get node");
    RSS(ref_recon_abs_value_hessian(ref_grid, metric), "abs val");
    RSS(ref_recon_roundoff_limit(metric, ref_grid),
        "floor metric eignvalues based on grid size and solution jitter");
    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local scale lp norm");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "compute metric from hessian");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("part scalar %s\n", in_scalar);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, in_scalar),
        "part scalar");
    REIS(1, ldim, "expected one scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "part scalar");

    if (ref_mpi_once(ref_mpi)) printf("reconstruct Hessian, compute metric\n");
    RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                      gradation, complexity),
        "lp norm");
    ref_free(scalar);
    ref_mpi_stopwatch_stop(ref_mpi, "compute metric");
  }

  if (buffer) {
    if (ref_mpi_once(ref_mpi)) printf("buffer at complexity %e\n", complexity);
    RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
        "buffer at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "buffer");
  }

  RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
  if (ref_mpi_once(ref_mpi))
    printf("actual complexity %e\n", current_complexity);
  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");

  ref_free(metric);

  if (ref_mpi_once(ref_mpi)) printf("gather %s\n", out_metric);
  RSS(ref_gather_metric(ref_grid, out_metric), "gather metric");
  ref_mpi_stopwatch_stop(ref_mpi, "gather metric");

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) multiscale_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS surface(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_file;
  char *in_file;
  char filename[1024];
  REF_GRID ref_grid = NULL;

  if (argc < 3) goto shutdown;
  in_file = argv[2];
  if (argc < 4) {
    RAS(strlen(in_file) < 1014, "input filename too long (>1014)");
    sprintf(filename, "%s-surf.tec", in_file);
    out_file = filename;
  } else {
    out_file = argv[3];
  }

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_file);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_file), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_file);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("gather %s\n", out_file);
    RSS(ref_gather_scalar_surf_tec(ref_grid, 0, NULL, NULL, out_file),
        "gather surf tec");
    ref_mpi_stopwatch_stop(ref_mpi, "gather");
  } else {
    if (ref_mpi_once(ref_mpi))
      printf("export " REF_GLOB_FMT " nodes to %s\n",
             ref_node_n_global(ref_grid_node(ref_grid)), out_file);
    RSS(ref_export_tec_surf(ref_grid, out_file), "export tec surf");
    ref_mpi_stopwatch_stop(ref_mpi, "export");
  }

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) surface_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS translate(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_file;
  char *in_file;
  REF_GRID ref_grid = NULL;
  REF_INT pos;

  if (argc < 4) goto shutdown;
  in_file = argv[2];
  out_file = argv[3];

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_file);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_file), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_file);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  RXS(ref_args_find(argc, argv, "--extrude", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_GRID twod_grid = ref_grid;
    if (ref_mpi_once(ref_mpi)) printf("extrude prims\n");
    RSS(ref_grid_extrude_twod(&ref_grid, twod_grid), "extrude");
    RSS(ref_grid_free(twod_grid), "free");
  }

  RXS(ref_args_find(argc, argv, "--zero-y-face", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_DBL deviation, total_deviation;
    REF_CELL ref_cell;
    REF_NODE ref_node = ref_grid_node(ref_grid);
    REF_INT faceid, cell, node, nodes[REF_CELL_MAX_SIZE_PER];
    if (pos + 1 >= argc) goto shutdown;
    faceid = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf("zero y of face %d\n", faceid);
    deviation = 0.0;
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (faceid == nodes[ref_cell_node_per(ref_cell)]) {
        each_ref_cell_cell_node(ref_cell, node) {
          deviation =
              MAX(deviation, ABS(ref_node_xyz(ref_node, 1, nodes[node])));
          ref_node_xyz(ref_node, 1, nodes[node]) = 0.0;
        }
      }
    }
    RSS(ref_mpi_max(ref_mpi, &deviation, &total_deviation, REF_DBL_TYPE),
        "mpi max");
    printf("max deviation %e\n", deviation);
  }

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi))
      printf("gather " REF_GLOB_FMT " nodes to %s\n",
             ref_node_n_global(ref_grid_node(ref_grid)), out_file);
    RSS(ref_gather_by_extension(ref_grid, out_file), "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "gather");
  } else {
    if (ref_mpi_once(ref_mpi))
      printf("export " REF_GLOB_FMT " nodes to %s\n",
             ref_node_n_global(ref_grid_node(ref_grid)), out_file);
    RSS(ref_export_by_extension(ref_grid, out_file), "export");
    ref_mpi_stopwatch_stop(ref_mpi, "export");
  }

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) translate_help(argv[0]);
  return REF_FAILURE;
}

static void echo_argv(int argc, char *argv[]) {
  int pos;
  printf("\n");
  for (pos = 0; pos < argc; pos++) printf(" %s", argv[pos]);
  printf("\n\n");
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT help_pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");
  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_once(ref_mpi)) {
    printf("refine %s, on or after 1.9.1\n", VERSION);
    echo_argv(argc, argv);
  }

  RXS(ref_args_find(argc, argv, "--help", &help_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY == help_pos) {
    RXS(ref_args_find(argc, argv, "-h", &help_pos), REF_NOT_FOUND,
        "arg search");
  }

  if (1 == argc || 1 == help_pos) {
    if (ref_mpi_once(ref_mpi)) usage(argv[0]);
    goto shutdown;
  }

  if (strncmp(argv[1], "a", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(adapt(ref_mpi, argc, argv), "adapt");
    } else {
      if (ref_mpi_once(ref_mpi)) adapt_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "b", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(bootstrap(ref_mpi, argc, argv), "bootstrap");
    } else {
      if (ref_mpi_once(ref_mpi)) bootstrap_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "f", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(fun3d(ref_mpi, argc, argv), "fun3d");
    } else {
      if (ref_mpi_once(ref_mpi)) fun3d_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "i", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(interpolate(ref_mpi, argc, argv), "interpolate");
    } else {
      if (ref_mpi_once(ref_mpi)) interpolate_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "l", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(loop(ref_mpi, argc, argv), "loop");
    } else {
      if (ref_mpi_once(ref_mpi)) loop_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "m", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(multiscale(ref_mpi, argc, argv), "multiscale");
    } else {
      if (ref_mpi_once(ref_mpi)) multiscale_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "s", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(surface(ref_mpi, argc, argv), "surface");
    } else {
      if (ref_mpi_once(ref_mpi)) surface_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "t", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(translate(ref_mpi, argc, argv), "translate");
    } else {
      if (ref_mpi_once(ref_mpi)) translate_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "v", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(vertex(ref_mpi, argc, argv), "translate");
    } else {
      if (ref_mpi_once(ref_mpi)) vertex_help(argv[0]);
      goto shutdown;
    }
  } else {
    if (ref_mpi_once(ref_mpi)) usage(argv[0]);
    goto shutdown;
  }

  ref_mpi_stopwatch_stop(ref_mpi, "done.");
shutdown:
  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
