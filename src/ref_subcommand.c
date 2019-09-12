
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

#include "ref_defs.h"

#include "ref_args.h"
#include "ref_mpi.h"

#include "ref_geom.h"
#include "ref_grid.h"

#include "ref_histogram.h"
#include "ref_metric.h"
#include "ref_split.h"
#include "ref_validation.h"

#include "ref_export.h"
#include "ref_gather.h"
#include "ref_import.h"
#include "ref_part.h"

static void usage(const char *name) {
  printf("usage: \n %s [--help] <command> [<args>]\n", name);
  printf("\n");
  printf("ref commands:\n");
  printf("  bootstrap Create initial grid from EGADS file\n");
  printf("  fill      Fill a surface shell mesh with a volume.\n");
  printf("  location  Report the locations of verticies in the mesh.\n");
  printf("  translate Convert mesh formats.\n");
}
static void bootstrap_help(const char *name) {
  printf("usage: \n %s boostrap project.egads [-t]\n", name);
  printf("  -t  tecplot movie of surface curvature adaptation\n");
  printf("        in files ref_gather_movie.tec and ref_gather_histo.tec\n");
  printf("\n");
}
static void fill_help(const char *name) {
  printf("usage: \n %s fill surface.meshb volume.meshb\n", name);
  printf("\n");
}
static void location_help(const char *name) {
  printf("usage: \n %s location input.meshb node_index node_index ...\n", name);
  printf("  node_index is zero-based\n");
  printf("\n");
}
static void translate_help(const char *name) {
  printf("usage: \n %s input_mesh.extension output_mesh.extension \n", name);
  printf("\n");
}

static REF_STATUS bootstrap(REF_MPI ref_mpi, int argc, char *argv[]) {
  size_t end_of_string;
  char project[1000];
  char filename[1024];
  REF_GRID ref_grid = NULL;
  REF_DBL params[3];
  REF_INT t_pos = REF_EMPTY;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref bootstrap is not parallel");
  }
  if (argc < 3) goto shutdown;
  end_of_string = MIN(1023, strlen(argv[2]));
  if (7 > end_of_string ||
      strncmp(&(argv[2][end_of_string - 6]), ".egads", 6) != 0)
    goto shutdown;
  strncpy(project, argv[2], end_of_string - 6);
  project[end_of_string - 6] = '\0';

  RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
  printf("loading %s.egads\n", project);
  RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[2]), "ld egads");
  ref_mpi_stopwatch_stop(ref_mpi, "egads load");

  printf("initial tessellation\n");
  RSS(ref_geom_egads_suggest_tess_params(ref_grid, params), "suggest params");
  RSS(ref_geom_egads_tess(ref_grid, params), "tess egads");
  ref_mpi_stopwatch_stop(ref_mpi, "egads tess");
  sprintf(filename, "%s-init.meshb", project);
  RSS(ref_export_by_extension(ref_grid, filename), "tess export");
  sprintf(filename, "%s-init-geom.tec", project);
  RSS(ref_geom_tec(ref_grid, filename), "geom export");
  sprintf(filename, "%s-init-surf.tec", project);
  RSS(ref_export_tec_surf(ref_grid, filename), "dbg surf");
  ref_mpi_stopwatch_stop(ref_mpi, "export init-surf");
  printf("verify topo\n");
  RSS(ref_geom_verify_topo(ref_grid), "adapt topo");
  printf("verify EGADS param\n");
  RSS(ref_geom_verify_param(ref_grid), "egads params");
  printf("constrain all\n");
  RSS(ref_geom_constrain_all(ref_grid), "constrain");
  printf("verify constrained param\n");
  RSS(ref_geom_verify_param(ref_grid), "constrained params");
  printf("verify manifold\n");
  RSS(ref_validation_boundary_manifold(ref_grid), "manifold");
  ref_mpi_stopwatch_stop(ref_mpi, "tess verification");

  RXS(ref_args_find(argc, argv, "-t", &t_pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != t_pos)
    RSS(ref_gather_tec_movie_record_button(ref_grid_gather(ref_grid), REF_TRUE),
        "movie on");

  RSS(ref_adapt_surf_to_geom(ref_grid), "ad");
  RSS(ref_geom_report_tri_area_normdev(ref_grid), "tri status");
  printf("verify topo\n");
  RSS(ref_geom_verify_topo(ref_grid), "adapt topo");
  printf("verify param\n");
  RSS(ref_geom_verify_param(ref_grid), "adapt params");
  ref_mpi_stopwatch_stop(ref_mpi, "surf verification");
  sprintf(filename, "%s-adapt-geom.tec", project);
  RSS(ref_geom_tec(ref_grid, filename), "geom export");
  sprintf(filename, "%s-adapt-surf.tec", project);
  RSS(ref_export_tec_surf(ref_grid, filename), "dbg surf");
  sprintf(filename, "%s-adapt-surf.meshb", project);
  printf("export %s\n", filename);
  RSS(ref_export_by_extension(ref_grid, filename), "surf export");
  ref_mpi_stopwatch_stop(ref_mpi, "export adapt surf");

  RSS(ref_geom_tetgen_volume(ref_grid), "tetgen surface to volume ");
  ref_mpi_stopwatch_stop(ref_mpi, "fill volume");

  RSS(ref_split_edge_geometry(ref_grid), "split geom");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "split geom");

  sprintf(filename, "%s-vol.meshb", project);
  printf("export %s\n", filename);
  RSS(ref_export_by_extension(ref_grid, filename), "vol export");
  ref_mpi_stopwatch_stop(ref_mpi, "export volume");

  RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "curvature");
  RSS(ref_histogram_quality(ref_grid), "gram");
  RSS(ref_histogram_ratio(ref_grid), "gram");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "histogram");

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  bootstrap_help(argv[0]);
  return REF_FAILURE;
}

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
  fill_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS location(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_file;
  REF_INT pos, global, local;
  REF_GRID ref_grid = NULL;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref location is not parallel");
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
  location_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS translate(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_file;
  char *in_file;
  REF_GRID ref_grid = NULL;

  if (argc < 4) goto shutdown;
  in_file = argv[2];
  out_file = argv[3];

  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_file);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_file), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_file);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("gather %s\n", out_file);
    RSS(ref_gather_by_extension(ref_grid, out_file), "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "gather");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("export %s\n", out_file);
    RSS(ref_export_by_extension(ref_grid, out_file), "export");
    ref_mpi_stopwatch_stop(ref_mpi, "export");
  }

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  translate_help(argv[0]);
  return REF_FAILURE;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT help_pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");
  ref_mpi_stopwatch_start(ref_mpi);

  RXS(ref_args_find(argc, argv, "--help", &help_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY == help_pos) {
    RXS(ref_args_find(argc, argv, "-h", &help_pos), REF_NOT_FOUND,
        "arg search");
  }

  if (1 == argc || 1 == help_pos) {
    usage(argv[0]);
    goto shutdown;
  }

  if (strncmp(argv[1], "b", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(bootstrap(ref_mpi, argc, argv), "bootstrap");
    } else {
      bootstrap_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "f", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(fill(ref_mpi, argc, argv), "fill");
    } else {
      fill_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "l", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(location(ref_mpi, argc, argv), "location");
    } else {
      location_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "t", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(translate(ref_mpi, argc, argv), "translate");
    } else {
      translate_help(argv[0]);
      goto shutdown;
    }
  } else {
    usage(argv[0]);
    goto shutdown;
  }

  ref_mpi_stopwatch_stop(ref_mpi, "done.");
shutdown:
  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
