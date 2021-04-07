
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
#include "ref_import.h"
#include "ref_iso.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_meshlink.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_phys.h"
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
  printf("  distance     Calculate wall distance (for turbulence model)\n");
  printf("  examine      Report mesh or solution file meta data.\n");
  /*printf("  grow         Fills surface mesh with volume to debug
   * bootstrap\n");*/
  printf("  interpolate  Interpolate a field from one mesh to another\n");
  printf("  loop         Multiscale metric, adapt, and interpolation.\n");
  printf("  multiscale   Compute a multiscale metric.\n");
  /*printf("  node       Reports location of a node by index\n");*/
  /*printf("  quilt      Construct effective EGADS model.\n");*/
  printf("  surface      Extract mesh surface.\n");
  printf("  translate    Convert mesh formats.\n");
  printf("  visualize    Convert solution formats.\n");
  printf("\n");
  printf("'ref <command> -h' provides details on a specific subcommand.\n");
}

static void option_uniform_help(void) {
  printf(
      "  --uniform box {ceil,floor} h0 decay_distance xmin ymin zmin "
      "xmax ymax zmax\n");
  printf(
      "  --uniform cyl {ceil,floor} h0 decay_distance x1 y1 z1 "
      "x2 y2 z2 r1 r2\n");
  printf("      decay_distance is negative to increase h with distance.\n");
  printf("      decay_distance is positive to decrease h with distance.\n");
}

static void option_auto_tprarms_help(void) {
  printf("  --auto-tparams {or combination of options} adjust .tParams\n");
  printf("        1:single edge, 2:chord violation, 4:face width (-1:all)\n");
}

static void adapt_help(const char *name) {
  printf("usage: \n %s adapt input_mesh.extension [<options>]\n", name);
  printf("  -x  output_mesh.extension\n");
  printf("  -g  geometry.egads\n");
  printf("  -m  metric.solb (geometry feature metric when missing)\n");
  printf("  --implied-complexity [complexity] imply metric from input mesh\n");
  printf("      and scale to complexity\n");
  printf("  --spalding [y+=1] [complexity]\n");
  printf("      construct a multiscale metric to control interpolation\n");
  printf("      error in u+ of Spalding's Law. Requires boundary conditions\n");
  printf("      via the --fun3d-mapbc or --viscous-tags options.\n");
  option_uniform_help();
  printf("  --fun3d-mapbc fun3d_format.mapbc\n");
  printf("  --viscous-tags <comma-separated list of viscous boundary tags>\n");
  printf("  --partitioner selects domain decomposition method.\n");
  printf("      2: ParMETIS graph partitioning.\n");
  printf("      3: Zoltan graph partitioning.\n");
  printf("      4: Zoltan recursive bisection.\n");
  printf("      5: native recursive bisection.\n");
  printf("\n");
}
static void bootstrap_help(const char *name) {
  printf("usage: \n %s bootstrap project.egads [-t]\n", name);
  printf("  -t  tecplot movie of surface curvature adaptation\n");
  printf("        in files ref_gather_movie.tec and ref_gather_histo.tec\n");
  printf("  --mesher {tetgen|aflr} volume mesher\n");
  printf("  --mesher-options \"<options>\" quoted mesher options.\n");
  option_auto_tprarms_help();
  printf("\n");
}
static void distance_help(const char *name) {
  printf("usage: \n %s distance input_mesh.extension distance.solb\n", name);
  printf("  --fun3d-mapbc fun3d_format.mapbc\n");
  printf("  --viscous-tags <comma-separated list of viscous boundary tags>\n");
  printf("\n");
}
static void examine_help(const char *name) {
  printf("usage: \n %s examine input_mesh_or_solb.extension\n", name);
  printf("\n");
}
static void grow_help(const char *name) {
  printf("usage: \n %s grow surface.meshb volume.meshb\n", name);
  printf("  --mesher {tetgen|aflr} volume mesher\n");
  printf("  --mesher-options \"<options>\" quoted mesher options.\n");
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
  printf("   --norm-power <power> multiscale metric norm power.\n");
  printf("       Default power is 2 (1 for goal-based metrics)\n");
  printf("   --gradation <gradation> (default -1)\n");
  printf("       positive: metric-space gradation stretching ratio.\n");
  printf("       negative: mixed-space gradation.\n");
  printf("   --partitioner <id> selects domain decomposition method.\n");
  printf("       2: ParMETIS graph partitioning.\n");
  printf("       3: Zoltan graph partitioning.\n");
  printf("       4: Zoltan recursive bisection.\n");
  printf("       5: native recursive bisection.\n");
  printf("   --mesh-extension <output mesh extension> (replaces lb8.ugrid).\n");
  printf("   --fixed-point <middle-string> \\\n");
  printf("       <first_timestep> <timestep_increment> <last_timestep>\n");
  printf("       where <input_project_name><middle-string>N.solb are\n");
  printf("       scalar fields and N is the timestep index.\n");
  printf("   --interpolant <type or file.solb> multiscale scalar field.\n");
  printf(
      "       Type is mach (default), "
      "incomp (incompressible vel magnitude),\n");
  printf("       htot, pressure, density, or temperature.\n");
  printf("       Read from file.solb, if not a recognized type.\n");
  printf("   --export-metric writes <input_project_name>-metric.solb.\n");
  printf("   --opt-goal metric of Loseille et al. AIAA 2007--4186.\n");
  printf("        Include flow and adjoint information in volume.solb.\n");
  printf("        Use --fun3d-mapbc or --viscous-tags with strong BCs.\n");
  printf("   --cons-visc <mach> <re> <temperature> see AIAA 2019--2947.\n");
  printf("        <mach> is reference Mach nubmer.\n");
  printf("        <re> is reference Reylonds number in grid units.\n");
  printf("        <temperature> is reference temperature in K.\n");
  printf("        Include flow and adjoint information in volume.solb.\n");
  printf("        Use --fun3d-mapbc or --viscous-tags with strong BCs.\n");
  printf("  --fun3d-mapbc fun3d_format.mapbc\n");
  printf("  --viscous-tags <comma-separated list of viscous boundary tags>\n");
  printf("  --deforming mesh flow solve, include xyz in *_volume.solb.\n");
  printf("  --buffer coarsens the metric approaching the x max boundary.\n");
  option_uniform_help();

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
  option_uniform_help();
  printf("   --hessian expects hessian.* in place of scalar.{solb,snap}.\n");
  printf("   --pcd <project.pcd> exports isotropic spacing.\n");
  printf("\n");
}
static void node_help(const char *name) {
  printf("usage: \n %s node input.meshb node_index node_index ...\n", name);
  printf("  node_index is zero-based\n");
  printf("\n");
}
static void quilt_help(const char *name) {
  printf("usage: \n %s quilt original.egads\n", name);
  printf("  originaleff.egads is output EGADS model with EBODY\n");
  option_auto_tprarms_help();
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
  printf("   --extrude a 2D mesh to single layer of prisms.\n");
  printf("       extrusion added implicitly for ugrid output files\n");
  printf("   --planes <N> extrude a 2D mesh to N layers of prisms.\n");
  printf("   --zero-y-face [face id] explicitly set y=0 on face id.\n");
  printf("\n");
}
static void visualize_help(const char *name) {
  printf(
      "usage: \n %s visualize input_mesh.extension input_solution.extension "
      "output_solution.extension\n",
      name);
  printf("\n");
  printf(
      "   --subtract <baseline_solution.extension> "
      "computes (input-baseline).\n");
  printf(
      "   --iso <0-based variable index> <threshold> "
      "extracts an isosurface.\n");
  printf("\n");
}

static REF_STATUS spalding_metric(REF_GRID ref_grid, REF_DICT ref_dict_bcs,
                                  REF_DBL spalding_yplus, REF_DBL complexity,
                                  int argc, char *argv[]) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_DBL *metric;
  REF_DBL *distance, *uplus, yplus;
  REF_INT node;
  REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
  REF_DBL gradation = 10.0;

  ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(distance, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(uplus, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  RSS(ref_phys_wall_distance(ref_grid, ref_dict_bcs, distance), "wall dist");
  ref_mpi_stopwatch_stop(ref_mpi, "wall distance");
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    RAS(ref_math_divisible(distance[node], spalding_yplus),
        "\nare viscous boundarys set with --viscous-tags or --fun3d-mapbc?"
        "\nwall distance not divisible by y+=1");
    yplus = distance[node] / spalding_yplus;
    RSS(ref_phys_spalding_uplus(yplus, &(uplus[node])), "uplus");
  }
  RSS(ref_recon_hessian(ref_grid, uplus, metric, reconstruction), "hess");
  RSS(ref_recon_roundoff_limit(metric, ref_grid),
      "floor metric eigenvalues based on grid size and solution jitter");
  RSS(ref_metric_local_scale(metric, NULL, ref_grid, 4),
      "local lp=4 norm scaling");
  ref_mpi_stopwatch_stop(ref_mpi, "spalding metric");
  RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                         complexity),
      "set complexity");
  RSS(ref_metric_parse(metric, ref_grid, argc, argv), "parse metric");
  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "node metric");
  ref_free(uplus);
  ref_free(distance);
  ref_free(metric);
  ref_mpi_stopwatch_stop(ref_mpi, "spalding gradation");
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid)) ||
      ref_geom_meshlinked(ref_grid_geom(ref_grid))) {
    RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
    ref_mpi_stopwatch_stop(ref_mpi, "crv const");
  }
  return REF_SUCCESS;
}

static REF_STATUS adapt(REF_MPI ref_mpi_orig, int argc, char *argv[]) {
  char *in_mesh = NULL;
  char *in_metric = NULL;
  char *in_egads = NULL;
  REF_GRID ref_grid = NULL;
  REF_MPI ref_mpi = ref_mpi_orig;
  REF_BOOL curvature_metric = REF_TRUE;
  REF_BOOL all_done = REF_FALSE;
  REF_BOOL all_done0 = REF_FALSE;
  REF_BOOL all_done1 = REF_FALSE;
  REF_INT pass, passes = 30;
  REF_INT opt, pos;
  REF_LONG ntet;
  REF_DICT ref_dict_bcs = NULL;
  REF_DBL spalding_yplus = -1.0;
  REF_DBL complexity = -1.0;

  if (argc < 3) goto shutdown;
  in_mesh = argv[2];

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_mesh);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_mesh), "part");
    ref_mpi = ref_grid_mpi(ref_grid); /* ref_grid made a deep copy */
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_mesh);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_mesh), "import");
    ref_mpi = ref_grid_mpi(ref_grid); /* ref_grid made a deep copy */
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

  RXS(ref_args_find(argc, argv, "--meshlink", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("meshlink with %s\n", argv[pos + 1]);
    RSS(ref_meshlink_open(ref_grid, argv[pos + 1]), "meshlink init");
    RSS(ref_meshlink_infer_orientation(ref_grid), "meshlink orient");
  } else {
    RXS(ref_args_char(argc, argv, "-g", &in_egads), REF_NOT_FOUND,
        "egads arg search");
    if (NULL != in_egads) {
      if (ref_mpi_once(ref_mpi)) printf("load egads from %s\n", in_egads);
      RSS(ref_egads_load(ref_grid_geom(ref_grid), in_egads), "load egads");
      if (ref_mpi_once(ref_mpi) && ref_geom_effective(ref_grid_geom(ref_grid)))
        printf("EBody Effective Body loaded\n");
      ref_mpi_stopwatch_stop(ref_mpi, "load egads");
    } else {
      if (0 < ref_geom_cad_data_size(ref_grid_geom(ref_grid))) {
        if (ref_mpi_once(ref_mpi))
          printf("load egadslite from .meshb byte stream\n");
        RSS(ref_egads_load(ref_grid_geom(ref_grid), NULL), "load egads");
        if (ref_mpi_once(ref_mpi) &&
            ref_geom_effective(ref_grid_geom(ref_grid)))
          printf("EBody Effective Body loaded\n");
        ref_mpi_stopwatch_stop(ref_mpi, "load egads");
      } else {
        if (ref_mpi_once(ref_mpi))
          printf("warning: no geometry loaded, assuming planar faces.\n");
      }
    }
  }

  if (ref_geom_model_loaded(ref_grid_geom(ref_grid))) {
    RSS(ref_cell_ncell(ref_grid_tet(ref_grid), ref_grid_node(ref_grid), &ntet),
        "global tets");
    if (0 == ntet) ref_grid_surf(ref_grid) = REF_TRUE;
    RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
  }
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid)) ||
      ref_geom_meshlinked(ref_grid_geom(ref_grid))) {
    RSS(ref_geom_verify_topo(ref_grid), "geom topo");
    RSS(ref_geom_verify_param(ref_grid), "geom param");
    ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");
  }

  RXS(ref_args_find(argc, argv, "--facelift", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("--facelift %s import\n", argv[pos + 1]);
    RSS(ref_facelift_import(ref_grid, argv[pos + 1]), "attach");
    ref_mpi_stopwatch_stop(ref_mpi, "facelift loaded");
  }

  RXS(ref_args_find(argc, argv, "--surrogate", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("--surrogate %s import\n", argv[pos + 1]);
    RSS(ref_facelift_surrogate(ref_grid, argv[pos + 1]), "attach");
    ref_mpi_stopwatch_stop(ref_mpi, "facelift loaded");
    if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
    if (ref_mpi_once(ref_mpi)) printf("verify constrained param\n");
    RSS(ref_geom_verify_param(ref_grid), "constrained params");
    ref_mpi_stopwatch_stop(ref_mpi, "verify param");
  }

  RXS(ref_args_find(argc, argv, "-t", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos)
    RSS(ref_gather_tec_movie_record_button(ref_grid_gather(ref_grid), REF_TRUE),
        "movie on");

  RXS(ref_args_find(argc, argv, "-s", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    passes = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf("-s %d adaptation passes\n", passes);
  }

  RXS(ref_args_find(argc, argv, "--partitioner", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    REF_INT part_int = atoi(argv[pos + 1]);
    ref_grid_partitioner(ref_grid) = (REF_MIGRATE_PARTIONER)part_int;
    if (ref_mpi_once(ref_mpi))
      printf("--partitioner %d partitioner\n",
             (int)ref_grid_partitioner(ref_grid));
  }

  RXS(ref_args_find(argc, argv, "--ratio-method", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    ref_grid_node(ref_grid)->ratio_method = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf("--ratio-method %d\n", ref_grid_node(ref_grid)->ratio_method);
  }

  RXS(ref_args_find(argc, argv, "--topo", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    ref_grid_adapt(ref_grid, watch_topo) = REF_TRUE;
    if (ref_mpi_once(ref_mpi)) printf("--topo checks active\n");
  }

  RXS(ref_args_char(argc, argv, "-m", &in_metric), REF_NOT_FOUND,
      "metric arg search");
  if (NULL != in_metric) {
    if (ref_mpi_once(ref_mpi)) printf("part metric %s\n", in_metric);
    RSS(ref_part_metric(ref_grid_node(ref_grid), in_metric), "part metric");
    curvature_metric = REF_FALSE;
    ref_mpi_stopwatch_stop(ref_mpi, "part metric");
  }

  RSS(ref_dict_create(&ref_dict_bcs), "make dict");

  RXS(ref_args_find(argc, argv, "--fun3d-mapbc", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *mapbc;
    mapbc = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("reading fun3d bc map %s\n", mapbc);
      RSS(ref_phys_read_mapbc(ref_dict_bcs, mapbc),
          "unable to read fun3d formatted mapbc");
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  RXS(ref_args_find(argc, argv, "--viscous-tags", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *tags;
    tags = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("parsing viscous tags\n");
      RSS(ref_phys_parse_tags(ref_dict_bcs, tags),
          "unable to parse viscous tags");
      printf(" %d viscous tags parsed\n", ref_dict_n(ref_dict_bcs));
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  RXS(ref_args_find(argc, argv, "--spalding", &pos), REF_NOT_FOUND,
      "metric arg search");
  if (REF_EMPTY != pos && pos < argc - 2) {
    if (0 == ref_dict_n(ref_dict_bcs)) {
      if (ref_mpi_once(ref_mpi))
        printf(
            "\nset viscous boundaries via --fun3d-mapbc or --viscous-tags "
            "to use --spalding\n\n");
      goto shutdown;
    }

    spalding_yplus = atof(argv[pos + 1]);
    complexity = atof(argv[pos + 2]);
    if (ref_mpi_once(ref_mpi))
      printf(" --spalding %e %f law of the wall metric\n", spalding_yplus,
             complexity);
    curvature_metric = REF_TRUE;
  }

  RXS(ref_args_find(argc, argv, "--implied-complexity", &pos), REF_NOT_FOUND,
      "metric arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    REF_DBL *metric;
    complexity = atof(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf(" --implied-complexity %f implied metric scaled to complexity\n",
             complexity);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_imply_from(metric, ref_grid), "imply metric");
    ref_mpi_stopwatch_stop(ref_mpi, "imply metric");
    RSS(ref_metric_set_complexity(metric, ref_grid, complexity),
        "scale metric");
    RSS(ref_metric_parse(metric, ref_grid, argc, argv), "parse metric");
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "node metric");
    ref_free(metric);
    curvature_metric = REF_FALSE;
  }

  if (curvature_metric) {
    if (spalding_yplus > 0.0) {
      RSS(spalding_metric(ref_grid, ref_dict_bcs, spalding_yplus, complexity,
                          argc, argv),
          "spalding");
    } else {
      RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
      ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
      RXS(ref_args_find(argc, argv, "--facelift-metric", &pos), REF_NOT_FOUND,
          "arg search");
      if (REF_EMPTY != pos && pos < argc - 1) {
        complexity = atof(argv[pos + 1]);
        if (ref_mpi_once(ref_mpi)) printf("--facelift-metric %f\n", complexity);
        RSS(ref_facelift_multiscale(ref_grid, complexity), "metric");
        ref_mpi_stopwatch_stop(ref_mpi, "facelift metric");
      }
    }
  } else {
    if (ref_geom_model_loaded(ref_grid_geom(ref_grid)) ||
        ref_geom_meshlinked(ref_grid_geom(ref_grid))) {
      RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
      RSS(ref_validation_cell_volume(ref_grid), "vol");
      ref_mpi_stopwatch_stop(ref_mpi, "crv const");
    }
    RSS(ref_grid_cache_background(ref_grid), "cache");
    ref_mpi_stopwatch_stop(ref_mpi, "cache background metric");
  }

  RSS(ref_validation_cell_volume(ref_grid), "vol");

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
    if (curvature_metric) {
      if (spalding_yplus > 0.0) {
        RSS(spalding_metric(ref_grid, ref_dict_bcs, spalding_yplus, complexity,
                            argc, argv),
            "spalding");
      } else {
        RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
        ref_mpi_stopwatch_stop(ref_mpi, "curvature metric");
        if (REF_EMPTY != pos && pos < argc - 1) {
          complexity = atof(argv[pos + 1]);
          if (ref_mpi_once(ref_mpi))
            printf("--facelift-metric %f\n", complexity);
          RSS(ref_facelift_multiscale(ref_grid, complexity), "metric");
          ref_mpi_stopwatch_stop(ref_mpi, "facelift metric");
        }
      }
    } else {
      RSS(ref_metric_synchronize(ref_grid), "sync with background");
      ref_mpi_stopwatch_stop(ref_mpi, "metric sync");
    }
    RSS(ref_validation_cell_volume(ref_grid), "vol");
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

  /* export via -x grid.ext and -f final-surf.tec and -q final-vol.plt */
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
    if (strcmp(argv[opt], "-q") == 0) {
      if (ref_mpi_once(ref_mpi))
        printf("gather final volume status %s\n", argv[opt + 1]);
      RSS(ref_gather_volume_status_tec(ref_grid, argv[opt + 1]), "gather -f");
    }
  }

  RSS(ref_dict_free(ref_dict_bcs), "free");
  RSS(ref_grid_free(ref_grid), "free");

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
  REF_INT facelift_pos = REF_EMPTY;
  REF_INT pos = REF_EMPTY;
  REF_INT auto_tparams = REF_EGADS_RECOMMENDED_TPARAM;
  const char *mesher = "tetgen";
  const char *mesher_options = NULL;
  REF_INT passes = 15;
  REF_INT self_intersections;
  REF_DBL *global_params = NULL;

  if (!ref_egads_allows_construction()) {
    if (ref_mpi_once(ref_mpi))
      printf("bootstrap requires EGADS(full) use ref or refmpifull\n");
    goto shutdown;
  }

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
  }
  RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[2]), "ld egads");
  if (ref_mpi_once(ref_mpi) && ref_geom_effective(ref_grid_geom(ref_grid)))
    printf("EBody Effective Body loaded\n");
  ref_mpi_stopwatch_stop(ref_mpi, "egads load");

  if (ref_mpi_once(ref_mpi)) {
    sprintf(filename, "%s-vol.mapbc", project);
    printf("extracting %s from 'bc_name' attributes\n", filename);
    if (REF_SUCCESS ==
        ref_egads_extract_mapbc(ref_grid_geom(ref_grid), filename)) {
      printf("%s extracted\n", filename);
    } else {
      printf("one or more 'bc_name' attributes not set, mapbc not written\n");
      printf(
          " All faces (or edges for 2D) should have bc_name attributes "
          "like so:\n");
      printf("         select face # all faces\n");
      printf("         attribute bc_name $4000_wall\n");
      printf("         select face 5\n");
      printf("         attribute bc_name $5000_farfield\n");
    }
  }

  RXS(ref_args_find(argc, argv, "--auto-tparams", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    auto_tparams = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf("--auto-tparams %d requested\n", auto_tparams);
    if (auto_tparams < 0) {
      auto_tparams = REF_EGADS_ALL_TPARAM;
      if (ref_mpi_once(ref_mpi))
        printf("--auto-tparams %d set to all\n", auto_tparams);
    }
  }

  RXS(ref_args_find(argc, argv, "--global", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 3) {
    ref_malloc(global_params, 3, REF_DBL);
    global_params[0] = atof(argv[pos + 1]);
    global_params[1] = atof(argv[pos + 2]);
    global_params[2] = atof(argv[pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf("initial tessellation, global param %f %f %f\n", global_params[0],
             global_params[1], global_params[2]);
  } else {
    if (ref_mpi_once(ref_mpi)) printf("initial tessellation, default param\n");
  }
  RSS(ref_egads_tess(ref_grid, auto_tparams, global_params), "tess egads");
  ref_free(global_params);
  global_params = NULL;
  ref_mpi_stopwatch_stop(ref_mpi, "egads tess");
  sprintf(filename, "%s-init-surf.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_export_tec_surf(ref_grid, filename), "dbg surf");
  ref_mpi_stopwatch_stop(ref_mpi, "export init-surf");
  sprintf(filename, "%s-init-geom.tec", project);
  if (ref_mpi_once(ref_mpi))
    RSS(ref_geom_tec(ref_grid, filename), "geom export");
  ref_mpi_stopwatch_stop(ref_mpi, "export init-geom");
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

  if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
  RSS(ref_geom_constrain_all(ref_grid), "constrain");
  ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
  if (ref_mpi_once(ref_mpi)) printf("verify constrained param\n");
  RSS(ref_geom_verify_param(ref_grid), "constrained params");
  ref_mpi_stopwatch_stop(ref_mpi, "verify param");

  if (REF_FALSE) {
    sprintf(filename, "%s-const-geom.tec", project);
    if (ref_mpi_once(ref_mpi))
      RSS(ref_geom_tec(ref_grid, filename), "geom export");
    ref_mpi_stopwatch_stop(ref_mpi, "export init-geom");
  }

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

  RXS(ref_args_find(argc, argv, "--mesher", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    mesher = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) printf("--mesher %s requested\n", mesher);
  }

  RXS(ref_args_find(argc, argv, "--mesher-options", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    mesher_options = argv[pos + 1];
    if (ref_mpi_once(ref_mpi))
      printf("--mesher-options %s requested\n", mesher_options);
  }

  RXS(ref_args_find(argc, argv, "-s", &s_pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != s_pos && s_pos < argc - 1) {
    passes = atoi(argv[s_pos + 1]);
    if (ref_mpi_once(ref_mpi))
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

  RXS(ref_args_find(argc, argv, "--facelift", &facelift_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != facelift_pos && facelift_pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) {
      printf("--facelift %s requested\n", argv[facelift_pos + 1]);
      RSS(ref_facelift_attach(ref_grid), "attach");
    }
    ref_mpi_stopwatch_stop(ref_mpi, "facelift attached");
    if (ref_mpi_once(ref_mpi)) {
      REF_FACELIFT ref_facelift = ref_geom_facelift(ref_grid_geom(ref_grid));
      RSS(ref_export_by_extension(ref_facelift_grid(ref_facelift),
                                  argv[facelift_pos + 1]),
          "facelift export");
      sprintf(filename, "%s-facelift-geom.tec", project);
      RSS(ref_facelift_tec(ref_facelift, filename), "facelift viz");
    }
    ref_mpi_stopwatch_stop(ref_mpi, "facelift dumped");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
    RSS(ref_geom_verify_param(ref_grid), "facelift params");
    ref_mpi_stopwatch_stop(ref_mpi, "verify param");
    RSS(ref_adapt_surf_to_geom(ref_grid, 3), "ad");
    ref_mpi_stopwatch_stop(ref_mpi, "untangle");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_mpi, "pack");
  }

  RXS(ref_args_find(argc, argv, "--surrogate", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    REF_FACELIFT ref_facelift;
    REF_GRID surrogate;
    REF_DBL gap;
    REF_GLOB nnode = 0;
    if (ref_mpi_once(ref_mpi)) {
      printf("--surrogate %s requested\n", argv[pos + 1]);
    }
    REIS(REF_MIGRATE_SINGLE, ref_grid_partitioner(ref_grid),
         "parallel implementation is incomplete");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    if (ref_mpi_once(ref_mpi)) printf("original gap %e\n", gap);
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_grid_deep_copy(&surrogate, ref_grid), "free grid");
      RSS(ref_geom_enrich3(surrogate), "enrich3");
      nnode = ref_node_n_global(ref_grid_node(surrogate));
      RSS(ref_mpi_bcast(ref_mpi, &nnode, 1, REF_GLOB_TYPE), "bcast nnode");
    } else {
      RSS(ref_grid_create(&surrogate, ref_mpi), "create grid");
      RSS(ref_mpi_bcast(ref_mpi, &nnode, 1, REF_GLOB_TYPE), "bcast nnode");
      RSS(ref_node_initialize_n_global(ref_grid_node(surrogate), nnode),
          "init nnodesg");
    }
    RSS(ref_migrate_replicate_ghost(surrogate), "replicant");
    RSS(ref_facelift_create(&ref_facelift, surrogate, REF_TRUE), "create");
    ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
    ref_mpi_stopwatch_stop(ref_mpi, "enrich attach surrogate");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    RSS(ref_geom_max_gap(ref_grid, &gap), "geom gap");
    if (ref_mpi_once(ref_mpi)) printf("surrogate gap %e\n", gap);
    if (ref_mpi_once(ref_mpi)) {
      printf("gather %s\n", argv[pos + 1]);
    }
    RSS(ref_gather_by_extension(surrogate, argv[pos + 1]), "gather surrogate");
    ref_mpi_stopwatch_stop(ref_mpi, "gather surrogate");
  }

  if (ref_geom_manifold(ref_grid_geom(ref_grid))) {
    if (strncmp(mesher, "t", 1) == 0) {
      if (ref_mpi_once(ref_mpi)) {
        printf("fill volume with TetGen\n");
        RSB(ref_geom_tetgen_volume(ref_grid, project, mesher_options),
            "tetgen surface to volume", {
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
        RSB(ref_geom_aflr_volume(ref_grid, project, mesher_options),
            "aflr surface to volume", {
              printf("probing adapted tessellation self-intersections\n");
              RSS(ref_dist_collisions(ref_grid, REF_TRUE, &self_intersections),
                  "bumps");
              printf("%d segment-triangle intersections detected.\n",
                     self_intersections);
            });
      }
      ref_mpi_stopwatch_stop(ref_mpi, "aflr volume");
    } else {
      if (ref_mpi_once(ref_mpi))
        printf("mesher '%s' not implemented\n", mesher);
      goto shutdown;
    }
    ref_grid_surf(ref_grid) = REF_FALSE; /* needed until vol mesher para */
    RSS(ref_validation_boundary_face(ref_grid),
        "boundary-interior connectivity");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "boundary-volume check");

    RSS(ref_split_edge_geometry(ref_grid), "split geom");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "split geom");
  } else {
    REF_BOOL flat;
    RSS(ref_egads_twod_flat_z(ref_grid_geom(ref_grid), &flat), "flatness");
    ref_grid_twod(ref_grid) = flat;
    if (ref_mpi_once(ref_mpi)) {
      if (ref_grid_twod(ref_grid)) {
        printf(" 2D mode inferred from model flatness\n");
      } else {
        printf(" model curved, assume 3D surface\n");
      }
    }
  }
  RSS(ref_node_synchronize_globals(ref_grid_node(ref_grid)), "sync glob");

  sprintf(filename, "%s-vol.meshb", project);
  if (ref_mpi_once(ref_mpi))
    printf("gather " REF_GLOB_FMT " nodes to %s\n",
           ref_node_n_global(ref_grid_node(ref_grid)), filename);
  RSS(ref_gather_by_extension(ref_grid, filename), "vol export");
  ref_mpi_stopwatch_stop(ref_mpi, "export volume");

  RSS(ref_validation_cell_volume(ref_grid), "vol");

  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) bootstrap_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS distance(REF_MPI ref_mpi, int argc, char *argv[]) {
  REF_GRID ref_grid;
  REF_DICT ref_dict_bcs;
  REF_DBL *distance;
  char *in_mesh = NULL;
  char *out_file = NULL;
  REF_INT pos;
  if (argc < 4) goto shutdown;
  in_mesh = argv[2];
  out_file = argv[3];

  RSS(ref_dict_create(&ref_dict_bcs), "create");

  RXS(ref_args_find(argc, argv, "--fun3d-mapbc", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *mapbc;
    mapbc = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("reading fun3d bc map %s\n", mapbc);
      RSS(ref_phys_read_mapbc(ref_dict_bcs, mapbc),
          "unable to read fun3d formatted mapbc");
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  /* delete this block when f3d uses --fun3d-mapbc */
  RXS(ref_args_find(argc, argv, "--fun3d", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *mapbc;
    mapbc = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      distance_help(argv[0]);
      printf(" use --fun3d-mapbc, --fun3d no longer supported \n");
      printf("reading fun3d bc map %s\n", mapbc);
      RSS(ref_phys_read_mapbc(ref_dict_bcs, mapbc),
          "unable to read fun3d formatted mapbc");
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  RXS(ref_args_find(argc, argv, "--viscous-tags", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *tags;
    tags = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("parsing viscous tags\n");
      RSS(ref_phys_parse_tags(ref_dict_bcs, tags),
          "unable to parse viscous tags");
      printf(" %d viscous tags parsed\n", ref_dict_n(ref_dict_bcs));
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  if (0 == ref_dict_n(ref_dict_bcs)) {
    if (ref_mpi_once(ref_mpi))
      printf(
          "\nno solid walls specified\n"
          "set viscous boundaries via --fun3d-mapbc or --viscous-tags\n\n");
    goto shutdown;
  }

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf("part %s\n", in_mesh);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, in_mesh), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("import %s\n", in_mesh);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_mesh), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
  }
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

  ref_malloc_init(distance, ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                  -1.0);
  RSS(ref_phys_wall_distance(ref_grid, ref_dict_bcs, distance), "store");
  ref_mpi_stopwatch_stop(ref_mpi, "wall distance");

  if (ref_mpi_once(ref_mpi)) printf("gather %s\n", out_file);
  RSS(ref_gather_scalar_by_extension(ref_grid, 1, distance, NULL, out_file),
      "gather");
  ref_mpi_stopwatch_stop(ref_mpi, "gather");

  ref_free(distance);
  ref_dict_free(ref_dict_bcs);
  ref_grid_free(ref_grid);

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) distance_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS examine(REF_MPI ref_mpi, int argc, char *argv[]) {
  if (argc < 3) goto shutdown;

  RSS(ref_import_examine_header(argv[2]), "examine header");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) examine_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS grow(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *out_file;
  char *in_file;
  char project[1000];
  size_t end_of_string;
  REF_GRID ref_grid = NULL;
  const char *mesher = "tetgen";
  const char *mesher_options = NULL;
  REF_INT pos;
  REF_INT self_intersections;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref grow is not parallel");
  }
  if (argc < 4) goto shutdown;
  in_file = argv[2];
  out_file = argv[3];
  end_of_string = MIN(1023, strlen(argv[2]));
  if (7 > end_of_string ||
      strncmp(&(argv[2][end_of_string - 6]), ".meshb", 6) != 0)
    goto shutdown;
  strncpy(project, argv[2], end_of_string - 6);
  project[end_of_string - 6] = '\0';

  printf("import %s\n", in_file);
  RSS(ref_import_by_extension(&ref_grid, ref_mpi, in_file), "load surface");

  RXS(ref_args_find(argc, argv, "--mesher", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    mesher = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) printf("--mesher %s requested\n", mesher);
  }

  RXS(ref_args_find(argc, argv, "--mesher-options", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    mesher_options = argv[pos + 1];
    if (ref_mpi_once(ref_mpi))
      printf("--mesher-options %s requested\n", mesher_options);
  }

  if (strncmp(mesher, "t", 1) == 0) {
    if (ref_mpi_once(ref_mpi)) {
      printf("fill volume with TetGen\n");
      RSB(ref_geom_tetgen_volume(ref_grid, project, mesher_options),
          "tetgen surface to volume", {
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
      RSB(ref_geom_aflr_volume(ref_grid, project, mesher_options),
          "aflr surface to volume", {
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
  RSS(ref_validation_boundary_face(ref_grid), "boundary-interior connectivity");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "boundary-volume check");

  RSS(ref_split_edge_geometry(ref_grid), "split geom");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "split geom");

  RSS(ref_node_synchronize_globals(ref_grid_node(ref_grid)), "sync glob");

  printf("export %s\n", out_file);
  RSS(ref_export_by_extension(ref_grid, out_file), "vol export");

  RSS(ref_validation_cell_volume(ref_grid), "vol");

  RSS(ref_grid_free(ref_grid), "create");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) grow_help(argv[0]);
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
  REF_INT ldim, persist_ldim;
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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(donor_grid)));

  if (ref_mpi_once(ref_mpi)) printf("part solution %s\n", donor_solb);
  RSS(ref_part_scalar(donor_grid, &ldim, &donor_solution, donor_solb),
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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(receipt_grid)));

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
    RSS(ref_part_scalar(receipt_grid, &persist_ldim, &receipt_solution,
                        persist_solb),
        "part solution");
    ref_mpi_stopwatch_stop(ref_mpi, "persist part solution");
    REIS(ldim, persist_ldim, "persist leading dimension different than donor");

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
    RSS(ref_grid_extrude_twod(&extruded_grid, receipt_grid, 2), "extrude");
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

static REF_STATUS initial_field_scalar(REF_GRID ref_grid, REF_INT ldim,
                                       REF_DBL *initial_field,
                                       const char *interpolant,
                                       REF_DBL *scalar) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_INT node;
  REF_DBL gamma = 1.4;
  REF_BOOL recognized = REF_FALSE;

  RSS(ref_validation_finite(ref_grid, ldim, initial_field), "init field");
  if (ref_mpi_once(ref_mpi)) printf("compute %s\n", interpolant);
  if (strcmp(interpolant, "incomp") == 0) {
    recognized = REF_TRUE;
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
  }
  if ((strcmp(interpolant, "mach") == 0) ||
      (strcmp(interpolant, "htot") == 0) ||
      (strcmp(interpolant, "pressure") == 0) ||
      (strcmp(interpolant, "density") == 0) ||
      (strcmp(interpolant, "temperature") == 0)) {
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
        recognized = REF_TRUE;
        scalar[node] = sqrt(mach2);
      } else if (strcmp(interpolant, "htot") == 0) {
        recognized = REF_TRUE;
        scalar[node] = temp * (1.0 / (gamma - 1.0)) + 0.5 * u2;
      } else if (strcmp(interpolant, "pressure") == 0) {
        recognized = REF_TRUE;
        scalar[node] = press;
      } else if (strcmp(interpolant, "density") == 0) {
        recognized = REF_TRUE;
        scalar[node] = rho;
      } else if (strcmp(interpolant, "temperature") == 0) {
        recognized = REF_TRUE;
        scalar[node] = temp;
      }
    }
    if (recognized)
      ref_mpi_stopwatch_stop(ref_mpi, "compute compressible scalar");
  }

  if (!recognized) {
    REF_INT solb_ldim;
    REF_DBL *solb_scalar;
    if (ref_mpi_once(ref_mpi))
      printf("opening %s as solb multiscale interpolant\n", interpolant);
    RSS(ref_part_scalar(ref_grid, &solb_ldim, &solb_scalar, interpolant),
        "unable to load interpolant scalar");
    REIS(1, solb_ldim, "expected one interpolant scalar");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = solb_scalar[node];
    }
    ref_free(solb_scalar);
    ref_mpi_stopwatch_stop(ref_mpi, "read interpolant from file");
  }

  return REF_SUCCESS;
}

static REF_STATUS fixed_point_metric(
    REF_DBL *metric, REF_GRID ref_grid, REF_INT first_timestep,
    REF_INT last_timestep, REF_INT timestep_increment, const char *in_project,
    const char *solb_middle, REF_RECON_RECONSTRUCTION reconstruction, REF_INT p,
    REF_DBL gradation, REF_DBL complexity) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_DBL *hess, *scalar;
  REF_INT timestep, total_timesteps;
  char solb_filename[1024];
  REF_DBL inv_total;
  REF_INT im, node;
  REF_INT fixed_point_ldim;

  ref_malloc(hess, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  total_timesteps = 0;
  for (timestep = first_timestep; timestep <= last_timestep;
       timestep += timestep_increment) {
    snprintf(solb_filename, 1024, "%s%s%d.solb", in_project, solb_middle,
             timestep);
    if (ref_mpi_once(ref_mpi))
      printf("read and hess recon for %s\n", solb_filename);
    RSS(ref_part_scalar(ref_grid, &fixed_point_ldim, &scalar, solb_filename),
        "unable to load scalar");
    REIS(1, fixed_point_ldim, "expected one scalar");
    RSS(ref_recon_hessian(ref_grid, scalar, hess, reconstruction), "hess");
    ref_free(scalar);
    total_timesteps++;
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (im = 0; im < 6; im++) {
        metric[im + 6 * node] += hess[im + 6 * node];
      }
    }
  }
  free(hess);
  ref_mpi_stopwatch_stop(ref_mpi, "all timesteps processed");

  RAS(0 < total_timesteps, "expected one or more timesteps");
  inv_total = 1.0 / (REF_DBL)total_timesteps;
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (im = 0; im < 6; im++) {
      metric[im + 6 * node] *= inv_total;
    }
  }
  RSS(ref_recon_roundoff_limit(metric, ref_grid),
      "floor metric eigenvalues based on grid size and solution jitter");
  RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
      "local lp norm scaling");
  ref_mpi_stopwatch_stop(ref_mpi, "local scale metric");
  RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                         complexity),
      "gradation at complexity");
  ref_mpi_stopwatch_stop(ref_mpi, "metric gradation and complexity");
  return REF_SUCCESS;
}

static REF_STATUS extract_displaced_xyz(REF_NODE ref_node, REF_INT *ldim,
                                        REF_DBL **initial_field,
                                        REF_DBL **displaced) {
  REF_INT i, node;

  ref_malloc(*displaced, 3 * ref_node_max(ref_node), REF_DBL);
  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < 3; i++) {
      (*displaced)[i + 3 * node] = (*initial_field)[i + (*ldim) * node];
    }
  }
  (*ldim) -= 3;
  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < (*ldim); i++) {
      (*initial_field)[i + (*ldim) * node] =
          (*initial_field)[i + 3 + ((*ldim) + 3) * node];
    }
  }
  ref_realloc(*initial_field, (*ldim) * ref_node_max(ref_node), REF_DBL);
  return REF_SUCCESS;
}

static REF_STATUS moving_fixed_point_metric(
    REF_DBL *metric, REF_GRID ref_grid, REF_INT first_timestep,
    REF_INT last_timestep, REF_INT timestep_increment, const char *in_project,
    const char *solb_middle, REF_RECON_RECONSTRUCTION reconstruction, REF_INT p,
    REF_DBL gradation, REF_DBL complexity) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *hess, *scalar;
  REF_DBL *jac, *x, *grad, *this_metric, *xyz, det;
  REF_INT timestep, total_timesteps;
  char solb_filename[1024];
  REF_DBL inv_total;
  REF_INT im, node;
  REF_INT fixed_point_ldim;
  REF_DBL *displaced;
  REF_INT i, j;

  ref_malloc(hess, 6 * ref_node_max(ref_node), REF_DBL);
  ref_malloc(this_metric, 6 * ref_node_max(ref_node), REF_DBL);
  ref_malloc(jac, 9 * ref_node_max(ref_node), REF_DBL);
  ref_malloc(x, ref_node_max(ref_node), REF_DBL);
  ref_malloc(grad, 3 * ref_node_max(ref_node), REF_DBL);
  ref_malloc(xyz, 3 * ref_node_max(ref_node), REF_DBL);

  total_timesteps = 0;
  for (timestep = first_timestep; timestep <= last_timestep;
       timestep += timestep_increment) {
    snprintf(solb_filename, 1024, "%s%s%d.solb", in_project, solb_middle,
             timestep);
    if (ref_mpi_once(ref_mpi))
      printf("read and hess recon for %s\n", solb_filename);
    RSS(ref_part_scalar(ref_grid, &fixed_point_ldim, &scalar, solb_filename),
        "unable to load scalar");
    REIS(4, fixed_point_ldim, "expected x,y,z and one scalar");
    RSS(extract_displaced_xyz(ref_node, &fixed_point_ldim, &scalar, &displaced),
        "disp");
    if (ref_grid_twod(ref_grid)) {
      each_ref_node_valid_node(ref_node, node) {
        displaced[1 + 3 * node] = displaced[2 + 3 * node];
        displaced[2 + 3 * node] = 0.0;
      }
    }
    for (j = 0; j < 3; j++) {
      each_ref_node_valid_node(ref_node, node) {
        x[node] = displaced[j + 3 * node];
      }
      RSS(ref_recon_gradient(ref_grid, x, grad, reconstruction), "recon x");
      if (ref_grid_twod(ref_grid)) {
        each_ref_node_valid_node(ref_node, node) { grad[2 + 3 * node] = 1.0; }
      }
      each_ref_node_valid_node(ref_node, node) {
        for (i = 0; i < 3; i++) {
          jac[i + 3 * j + 9 * node] = grad[i + 3 * node];
        }
      }
    }

    each_ref_node_valid_node(ref_node, node) {
      for (i = 0; i < 3; i++) {
        xyz[i + 3 * node] = ref_node_xyz(ref_node, i, node);
        ref_node_xyz(ref_node, i, node) = displaced[i + 3 * node];
      }
    }
    RSS(ref_recon_hessian(ref_grid, scalar, hess, reconstruction), "hess");
    RSS(ref_recon_roundoff_limit(hess, ref_grid),
        "floor metric eigenvalues based on grid size and solution jitter");
    each_ref_node_valid_node(ref_node, node) {
      for (i = 0; i < 3; i++) {
        ref_node_xyz(ref_node, i, node) = xyz[i + 3 * node];
      }
    }

    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_matrix_jac_m_jact(&(jac[9 * node]), &(hess[6 * node]),
                                &(this_metric[6 * node])),
          "J M J^t");

      RSS(ref_matrix_det_gen(3, &(jac[9 * node]), &det), "gen det");
      for (i = 0; i < 6; i++) {
        this_metric[i + 6 * node] *= pow(ABS(det), 1.0 / (REF_DBL)p);
      }
    }

    total_timesteps++;
    each_ref_node_valid_node(ref_node, node) {
      for (im = 0; im < 6; im++) {
        metric[im + 6 * node] += this_metric[im + 6 * node];
      }
    }

    ref_free(displaced);
    ref_free(scalar);
  }
  free(xyz);
  free(grad);
  free(x);
  free(jac);
  free(this_metric);
  free(hess);
  ref_mpi_stopwatch_stop(ref_mpi, "all timesteps processed");

  RAS(0 < total_timesteps, "expected one or more timesteps");
  inv_total = 1.0 / (REF_DBL)total_timesteps;
  each_ref_node_valid_node(ref_node, node) {
    for (im = 0; im < 6; im++) {
      metric[im + 6 * node] *= inv_total;
    }
  }
  RSS(ref_recon_roundoff_limit(metric, ref_grid),
      "floor metric eigenvalues based on grid size and solution jitter");
  RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
      "local lp norm scaling");
  ref_mpi_stopwatch_stop(ref_mpi, "local scale metric");
  RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                         complexity),
      "gradation at complexity");
  ref_mpi_stopwatch_stop(ref_mpi, "metric gradation and complexity");
  return REF_SUCCESS;
}

static REF_STATUS remove_initial_field_adjoint(REF_NODE ref_node, REF_INT *ldim,
                                               REF_DBL **initial_field) {
  REF_INT i, node;
  RAS((*ldim) % 2 == 0, "volume field should have a even leading dimension");
  (*ldim) /= 2;
  each_ref_node_valid_node(ref_node, node) {
    if (0 != node) {
      for (i = 0; i < (*ldim); i++) {
        (*initial_field)[i + (*ldim) * node] =
            (*initial_field)[i + 2 * (*ldim) * node];
      }
    }
  }
  ref_realloc(*initial_field, (*ldim) * ref_node_max(ref_node), REF_DBL);
  return REF_SUCCESS;
}

static REF_STATUS mask_strong_bc_adjoint(REF_GRID ref_grid,
                                         REF_DICT ref_dict_bcs, REF_INT ldim,
                                         REF_DBL *prim_dual) {
  REF_BOOL *replace;
  ref_malloc(replace, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_BOOL);
  RSS(ref_phys_mask_strong_bcs(ref_grid, ref_dict_bcs, replace, ldim), "mask");
  RSS(ref_recon_extrapolate_kexact(ref_grid, prim_dual, replace, ldim),
      "extrapolate kexact");
  ref_free(replace);

  return REF_SUCCESS;
}

static REF_STATUS flip_twod_yz(REF_NODE ref_node, REF_INT ldim,
                               REF_DBL *field) {
  REF_INT node;
  REF_DBL temp;
  REF_INT nequ;

  nequ = 0;
  if (ldim > 5 && 0 == ldim % 5) nequ = 5;
  if (ldim > 6 && 0 == ldim % 6) nequ = 6;

  each_ref_node_valid_node(ref_node, node) {
    if (ldim >= 5) {
      temp = field[2 + ldim * node];
      field[2 + ldim * node] = field[3 + ldim * node];
      field[3 + ldim * node] = temp;
    }
    if (nequ > 0) {
      temp = field[nequ + 2 + ldim * node];
      field[nequ + 2 + ldim * node] = field[nequ + 3 + ldim * node];
      field[nequ + 3 + ldim * node] = temp;
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS loop(REF_MPI ref_mpi_orig, int argc, char *argv[]) {
  char *in_project = NULL;
  char *out_project = NULL;
  char filename[1024];
  REF_GRID ref_grid = NULL;
  REF_MPI ref_mpi = ref_mpi_orig;
  REF_GRID extruded_grid = NULL;
  REF_BOOL all_done = REF_FALSE;
  REF_BOOL all_done0 = REF_FALSE;
  REF_BOOL all_done1 = REF_FALSE;
  REF_INT pass, passes = 30;
  REF_INT ldim;
  REF_DBL *initial_field, *ref_field, *extruded_field = NULL, *scalar, *metric;
  REF_DBL *displaced = NULL;
  REF_INT p = 2;
  REF_DBL gradation = -1.0, complexity;
  REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
  REF_BOOL buffer = REF_FALSE;
  REF_BOOL multiscale_metric;
  REF_DICT ref_dict_bcs = NULL;
  REF_INT pos;
  REF_INT fixed_point_pos, deforming_pos;
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
  RXS(ref_args_find(argc, argv, "--opt-goal", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    p = 1;
  }
  RXS(ref_args_find(argc, argv, "--cons-euler", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    p = 1;
  }
  RXS(ref_args_find(argc, argv, "--cons-visc", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos + 3 < argc) {
    p = 1;
  }
  RXS(ref_args_find(argc, argv, "--norm-power", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    if (pos >= argc - 1) {
      if (ref_mpi_once(ref_mpi))
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
      if (ref_mpi_once(ref_mpi))
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

  RXS(ref_args_find(argc, argv, "-s", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    passes = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf("-s %d adaptation passes\n", passes);
  }

  RSS(ref_dict_create(&ref_dict_bcs), "make dict");

  RXS(ref_args_find(argc, argv, "--fun3d-mapbc", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *mapbc;
    mapbc = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("reading fun3d bc map %s\n", mapbc);
      RSS(ref_phys_read_mapbc(ref_dict_bcs, mapbc),
          "unable to read fun3d formatted mapbc");
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  RXS(ref_args_find(argc, argv, "--viscous-tags", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    const char *tags;
    tags = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) {
      printf("parsing viscous tags\n");
      RSS(ref_phys_parse_tags(ref_dict_bcs, tags),
          "unable to parse viscous tags");
      printf(" %d viscous tags parsed\n", ref_dict_n(ref_dict_bcs));
    }
    RSS(ref_dict_bcast(ref_dict_bcs, ref_mpi), "bcast");
  }

  sprintf(filename, "%s.meshb", in_project);
  if (ref_mpi_once(ref_mpi)) printf("part mesh %s\n", filename);
  RSS(ref_part_by_extension(&ref_grid, ref_mpi, filename), "part");
  ref_mpi = ref_grid_mpi(ref_grid); /* ref_grid made a deep copy */
  ref_mpi_stopwatch_stop(ref_mpi, "part");
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

  RXS(ref_args_find(argc, argv, "--partitioner", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    REF_INT part_int = atoi(argv[pos + 1]);
    ref_grid_partitioner(ref_grid) = (REF_MIGRATE_PARTIONER)part_int;
    if (ref_mpi_once(ref_mpi))
      printf("--partitioner %d partitioner\n",
             (int)ref_grid_partitioner(ref_grid));
  }

  RXS(ref_args_find(argc, argv, "--ratio-method", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    ref_grid_node(ref_grid)->ratio_method = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf("--ratio-method %d\n", ref_grid_node(ref_grid)->ratio_method);
  }

  RXS(ref_args_find(argc, argv, "--topo", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    ref_grid_adapt(ref_grid, watch_topo) = REF_TRUE;
    if (ref_mpi_once(ref_mpi)) printf("--topo checks active\n");
  }

  RXS(ref_args_find(argc, argv, "--meshlink", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("meshlink with %s\n", argv[pos + 1]);
    RSS(ref_meshlink_open(ref_grid, argv[pos + 1]), "meshlink init");
    RSS(ref_meshlink_infer_orientation(ref_grid), "meshlink orient");
  } else {
    RXS(ref_args_find(argc, argv, "--egads", &pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY != pos && pos < argc - 1) {
      if (ref_mpi_once(ref_mpi)) printf("load egads from %s\n", argv[pos + 1]);
      RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[pos + 1]), "load egads");
      if (ref_mpi_once(ref_mpi) && ref_geom_effective(ref_grid_geom(ref_grid)))
        printf("EBody Effective Body loaded\n");
      ref_mpi_stopwatch_stop(ref_mpi, "load egads");
    } else {
      if (0 < ref_geom_cad_data_size(ref_grid_geom(ref_grid))) {
        if (ref_mpi_once(ref_mpi))
          printf("load egadslite from .meshb byte stream\n");
        RSS(ref_egads_load(ref_grid_geom(ref_grid), NULL), "load egads");
        if (ref_mpi_once(ref_mpi) &&
            ref_geom_effective(ref_grid_geom(ref_grid)))
          printf("EBody Effective Body loaded\n");
        ref_mpi_stopwatch_stop(ref_mpi, "load egadslite cad data");
      } else {
        if (ref_mpi_once(ref_mpi))
          printf("warning: no geometry loaded, assuming planar faces.\n");
      }
    }
  }

  RXS(ref_args_find(argc, argv, "--facelift", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("--facelift %s import\n", argv[pos + 1]);
    RSS(ref_facelift_import(ref_grid, argv[pos + 1]), "attach");
    ref_mpi_stopwatch_stop(ref_mpi, "facelift loaded");
  }

  RXS(ref_args_find(argc, argv, "--surrogate", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    if (ref_mpi_once(ref_mpi)) printf("--surrogate %s import\n", argv[pos + 1]);
    RSS(ref_facelift_surrogate(ref_grid, argv[pos + 1]), "attach");
    ref_mpi_stopwatch_stop(ref_mpi, "facelift loaded");
    if (ref_mpi_once(ref_mpi)) printf("constrain all\n");
    RSS(ref_geom_constrain_all(ref_grid), "constrain");
    ref_mpi_stopwatch_stop(ref_mpi, "constrain param");
    if (ref_mpi_once(ref_mpi)) printf("verify constrained param\n");
    RSS(ref_geom_verify_param(ref_grid), "constrained params");
    ref_mpi_stopwatch_stop(ref_mpi, "verify param");
  }

  RXS(ref_args_find(argc, argv, "--usm3d", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY == pos) {
    sprintf(filename, "%s_volume.solb", in_project);
    if (ref_mpi_once(ref_mpi)) printf("part scalar %s\n", filename);
    RSS(ref_part_scalar(ref_grid, &ldim, &initial_field, filename),
        "part scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "part scalar");
  } else {
    sprintf(filename, "%s_volume.plt", in_project);
    if (ref_mpi_once(ref_mpi)) printf("reconstruct scalar %s\n", filename);
    RSS(ref_part_scalar(ref_grid, &ldim, &initial_field, filename),
        "part scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "reconstruct scalar");
  }

  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("flip initial_field v-w for twod\n");
    RSS(flip_twod_yz(ref_grid_node(ref_grid), ldim, initial_field), "flip");
  }

  RXS(ref_args_find(argc, argv, "--fixed-point", &fixed_point_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--deforming", &deforming_pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != deforming_pos && REF_EMPTY == fixed_point_pos) {
    if (ref_mpi_once(ref_mpi)) printf("extract xyz displacement\n");
    RSS(extract_displaced_xyz(ref_grid_node(ref_grid), &ldim, &initial_field,
                              &displaced),
        "extract displacments");
  }

  if (ref_mpi_once(ref_mpi)) {
    printf("complexity %f\n", complexity);
    printf("Lp=%d\n", p);
    printf("gradation %f\n", gradation);
    printf("reconstruction %d\n", (int)reconstruction);
  }

  ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                  0.0);

  multiscale_metric = REF_TRUE;
  RXS(ref_args_find(argc, argv, "--opt-goal", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    multiscale_metric = REF_FALSE;
    if (ref_mpi_once(ref_mpi)) printf("--opt-goal metric construction\n");
    RSS(mask_strong_bc_adjoint(ref_grid, ref_dict_bcs, ldim, initial_field),
        "maks");
    RSS(ref_metric_belme_gfe(metric, ref_grid, ldim, initial_field,
                             reconstruction),
        "add nonlinear terms");
    RSS(ref_recon_roundoff_limit(metric, ref_grid),
        "floor metric eigenvalues based on grid size and solution jitter");
    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local scale lp norm");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    RSS(remove_initial_field_adjoint(ref_grid_node(ref_grid), &ldim,
                                     &initial_field),
        "rm adjoint");
  }
  RXS(ref_args_find(argc, argv, "--cons-euler", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_DBL *g;
    multiscale_metric = REF_FALSE;
    if (ref_mpi_once(ref_mpi)) printf("--cons-euler metric construction\n");
    RSS(mask_strong_bc_adjoint(ref_grid, ref_dict_bcs, ldim, initial_field),
        "maks");
    ref_malloc_init(g, 5 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    RSS(ref_metric_cons_euler_g(g, ref_grid, ldim, initial_field,
                                reconstruction),
        "cons euler g weights");
    RSS(ref_metric_cons_assembly(metric, g, ref_grid, ldim, initial_field,
                                 reconstruction),
        "cons metric assembly");
    ref_free(g);
    RSS(ref_recon_roundoff_limit(metric, ref_grid),
        "floor metric eigenvalues based on grid size and solution jitter");
    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local scale lp norm");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    RSS(remove_initial_field_adjoint(ref_grid_node(ref_grid), &ldim,
                                     &initial_field),
        "rm adjoint");
  }
  RXS(ref_args_find(argc, argv, "--cons-visc", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos + 3 < argc) {
    REF_DBL *g;
    REF_DBL mach, re, temperature;
    multiscale_metric = REF_FALSE;
    mach = atof(argv[pos + 1]);
    re = atof(argv[pos + 2]);
    temperature = atof(argv[pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf(
          "--cons-visc %.3f Mach %.2e Re %.2f temperature metric "
          "construction\n",
          mach, re, temperature);
    RSS(mask_strong_bc_adjoint(ref_grid, ref_dict_bcs, ldim, initial_field),
        "maks");
    ref_malloc_init(g, 5 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    RSS(ref_metric_cons_euler_g(g, ref_grid, ldim, initial_field,
                                reconstruction),
        "cons euler g weights");
    RSS(ref_metric_cons_viscous_g(g, ref_grid, ldim, initial_field, mach, re,
                                  temperature, reconstruction),
        "cons viscous g weights");
    RSS(ref_metric_cons_assembly(metric, g, ref_grid, ldim, initial_field,
                                 reconstruction),
        "cons metric assembly");
    ref_free(g);
    RSS(ref_recon_roundoff_limit(metric, ref_grid),
        "floor metric eigenvalues based on grid size and solution jitter");
    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local scale lp norm");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    RSS(remove_initial_field_adjoint(ref_grid_node(ref_grid), &ldim,
                                     &initial_field),
        "rm adjoint");
  }
  RXS(ref_args_find(argc, argv, "--fixed-point", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos + 4 < argc) {
    REF_INT first_timestep, last_timestep, timestep_increment;
    const char *solb_middle;
    multiscale_metric = REF_FALSE;
    solb_middle = argv[pos + 1];
    first_timestep = atoi(argv[pos + 2]);
    timestep_increment = atoi(argv[pos + 3]);
    last_timestep = atoi(argv[pos + 4]);
    if (ref_mpi_once(ref_mpi)) {
      printf("--fixed-point\n");
      printf("    %s%s solb project\n", in_project, solb_middle);
      printf("    timesteps [%d ... %d ... %d]\n", first_timestep,
             timestep_increment, last_timestep);
    }
    RXS(ref_args_find(argc, argv, "--deforming", &deforming_pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY == deforming_pos) {
      RSS(fixed_point_metric(metric, ref_grid, first_timestep, last_timestep,
                             timestep_increment, in_project, solb_middle,
                             reconstruction, p, gradation, complexity),
          "fixed point");
    } else {
      RSS(moving_fixed_point_metric(metric, ref_grid, first_timestep,
                                    last_timestep, timestep_increment,
                                    in_project, solb_middle, reconstruction, p,
                                    gradation, complexity),
          "fixed point");
    }
  }
  if (multiscale_metric) {
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    if (ref_mpi_once(ref_mpi))
      printf("computing interpolant %s for multiscale metric\n", interpolant);
    RSS(initial_field_scalar(ref_grid, ldim, initial_field, interpolant,
                             scalar),
        "field metric");

    RXS(ref_args_find(argc, argv, "--deforming", &pos), REF_NOT_FOUND,
        "arg search");
    if (REF_EMPTY != pos) {
      if (ref_mpi_once(ref_mpi))
        printf("reconstruct Hessian, compute metric\n");
      RSS(ref_metric_moving_multiscale(metric, ref_grid, displaced, scalar,
                                       reconstruction, p, gradation,
                                       complexity),
          "lp norm");
      ref_mpi_stopwatch_stop(ref_mpi, "deforming metric");
    } else {
      if (ref_mpi_once(ref_mpi))
        printf("reconstruct Hessian, compute metric\n");
      RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                        gradation, complexity),
          "lp norm");
      ref_mpi_stopwatch_stop(ref_mpi, "multiscale metric");
    }
    ref_free(scalar);
  }

  if (buffer) {
    if (ref_mpi_once(ref_mpi)) printf("buffer at complexity %e\n", complexity);
    RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
        "buffer at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "buffer");
  }

  RXS(ref_args_find(argc, argv, "--uniform", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    RSS(ref_metric_parse(metric, ref_grid, argc, argv), "parse uniform");
  }

  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
  ref_free(metric);

  RXS(ref_args_find(argc, argv, "--export-metric", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    sprintf(filename, "%s-metric.solb", in_project);
    if (ref_mpi_once(ref_mpi)) printf("export metric to %s\n", filename);
    RSS(ref_gather_metric(ref_grid, filename), "export metric");
    ref_mpi_stopwatch_stop(ref_mpi, "export metric");
  }

  ref_grid_surf(ref_grid) = ref_grid_twod(ref_grid);
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid))) {
    RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps; UV degen");
  }
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid)) ||
      ref_geom_meshlinked(ref_grid_geom(ref_grid))) {
    RSS(ref_geom_verify_topo(ref_grid), "geom topo");
    RSS(ref_geom_verify_param(ref_grid), "geom param");
    ref_mpi_stopwatch_stop(ref_mpi, "geom assoc");
    RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
    ref_mpi_stopwatch_stop(ref_mpi, "crv const");
  }
  RSS(ref_grid_cache_background(ref_grid), "cache");
  RSS(ref_node_store_aux(ref_grid_node(ref_grid_background(ref_grid)), ldim,
                         initial_field),
      "store init field with background");
  ref_free(initial_field);
  ref_mpi_stopwatch_stop(ref_mpi, "cache background metric and field");

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
    RSS(ref_metric_synchronize(ref_grid), "sync with background");
    ref_mpi_stopwatch_stop(ref_mpi, "metric sync");
    RSS(ref_validation_cell_volume(ref_grid), "vol");
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

  RXS(ref_args_find(argc, argv, "--export-metric", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    sprintf(filename, "%s-final-metric.solb", out_project);
    if (ref_mpi_once(ref_mpi)) printf("export metric to %s\n", filename);
    RSS(ref_gather_metric(ref_grid, filename), "export metric");
    ref_mpi_stopwatch_stop(ref_mpi, "export metric");
  }

  sprintf(filename, "%s.meshb", out_project);
  if (ref_mpi_once(ref_mpi))
    printf("gather " REF_GLOB_FMT " nodes to %s\n",
           ref_node_n_global(ref_grid_node(ref_grid)), filename);
  RSS(ref_gather_by_extension(ref_grid, filename), "gather .meshb");
  ref_mpi_stopwatch_stop(ref_mpi, "gather meshb");

  sprintf(filename, "%s.%s", out_project, mesh_extension);
  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("extrude twod\n");
    RSS(ref_grid_extrude_twod(&extruded_grid, ref_grid, 2), "extrude");
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
           ldim,
           ref_node_n_global(ref_grid_node(ref_grid_background(ref_grid))),
           ref_node_n_global(ref_grid_node(ref_grid)));
  }

  if (ref_mpi_once(ref_mpi)) printf("interpolate receptor nodes\n");
  ref_malloc_init(ref_field, ldim * ref_node_max(ref_grid_node(ref_grid)),
                  REF_DBL, 0.0);
  RSS(ref_node_extract_aux(ref_grid_node(ref_grid_background(ref_grid)), &ldim,
                           &initial_field),
      "store init field with background");
  RSS(ref_validation_finite(ref_grid_background(ref_grid), ldim, initial_field),
      "recall background field");

  RSS(ref_interp_scalar(ref_grid_interp(ref_grid), ldim, initial_field,
                        ref_field),
      "interp scalar");
  RSS(ref_validation_finite(ref_grid, ldim, ref_field), "interp field");
  ref_free(initial_field);
  /* free interp and background grid */
  RSS(ref_grid_free(ref_grid_background(ref_grid)),
      "free cached background grid");
  RSS(ref_interp_free(ref_grid_interp(ref_grid)), "interp free");
  ref_grid_interp(ref_grid) = NULL;
  ref_mpi_stopwatch_stop(ref_mpi, "interp");

  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("flip ref_field v-w for twod\n");
    RSS(flip_twod_yz(ref_grid_node(ref_grid), ldim, ref_field), "flip");
  }

  if (ref_grid_twod(ref_grid)) {
    if (ref_mpi_once(ref_mpi)) printf("extruding field of %d\n", ldim);
    ref_malloc(extruded_field,
               ldim * ref_node_max(ref_grid_node(extruded_grid)), REF_DBL);
    RSS(ref_grid_extrude_field(ref_grid, ldim, ref_field, extruded_grid,
                               extruded_field),
        "extrude field");
    RSS(ref_validation_finite(extruded_grid, ldim, extruded_field),
        "extruded field");
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

  RSS(ref_dict_free(ref_dict_bcs), "free");
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
      if (ref_mpi_once(ref_mpi))
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
      if (ref_mpi_once(ref_mpi))
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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

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
        "floor metric eigenvalues based on grid size and solution jitter");
    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local scale lp norm");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "compute metric from hessian");
  } else {
    if (ref_mpi_once(ref_mpi)) printf("part scalar %s\n", in_scalar);
    RSS(ref_part_scalar(ref_grid, &ldim, &scalar, in_scalar), "part scalar");
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

  RXS(ref_args_find(argc, argv, "--uniform", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    RSS(ref_metric_parse(metric, ref_grid, argc, argv), "parse uniform");
  }

  RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
  if (ref_mpi_once(ref_mpi))
    printf("actual complexity %e\n", current_complexity);
  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");

  RXS(ref_args_find(argc, argv, "--pcd", &pos), REF_NOT_FOUND, "arg search");
  printf("pos %d arc %d\n", pos, argc);
  if (REF_EMPTY != pos && pos + 1 < argc) {
    REF_DBL *hh;
    const char *title[] = {"spacing", "decay"};
    ref_malloc(hh, 2 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_isotropic(metric, ref_grid, hh), "iso");
    ref_mpi_stopwatch_stop(ref_mpi, "isotropic");
    if (ref_mpi_once(ref_mpi)) printf("gather %s\n", argv[pos + 1]);
    RSS(ref_gather_scalar_by_extension(ref_grid, 2, hh, title, argv[pos + 1]),
        "dump hh");
    ref_free(hh);
  }

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

static REF_STATUS node(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_file;
  REF_INT pos, global, local;
  REF_GRID ref_grid = NULL;

  if (ref_mpi_para(ref_mpi)) {
    RSS(REF_IMPLEMENT, "ref node is not parallel");
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
  if (ref_mpi_once(ref_mpi)) quilt_help(argv[0]);
  return REF_FAILURE;
}

static REF_STATUS quilt(REF_MPI ref_mpi, int argc, char *argv[]) {
  REF_GEOM ref_geom;
  char *input_egads;
  size_t end_of_string;
  char project[1000];
  char output_egads[1024];
  REF_INT pos;
  REF_DBL *global_params = NULL;
  REF_INT auto_tparams = REF_EGADS_RECOMMENDED_TPARAM;

  if (argc < 3) goto shutdown;
  input_egads = argv[2];

  RAS(ref_egads_allows_construction(),
      "EGADS not linked with OpenCASCADE, required to load model")
  RAS(ref_egads_allows_effective(), "EGADS does not support Effective Geometry")

  end_of_string = MIN(1023, strlen(input_egads));
  RAS((7 < end_of_string &&
       strncmp(&(input_egads[end_of_string - 6]), ".egads", 6) == 0),
      ".egads extension missing");
  strncpy(project, input_egads, end_of_string - 6);
  project[end_of_string - 6] = '\0';
  sprintf(output_egads, "%s-eff.egads", project);

  RXS(ref_args_find(argc, argv, "--global", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 3) {
    ref_malloc(global_params, 3, REF_DBL);
    global_params[0] = atof(argv[pos + 1]);
    global_params[1] = atof(argv[pos + 2]);
    global_params[2] = atof(argv[pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf("initial tessellation, global param %f %f %f\n", global_params[0],
             global_params[1], global_params[2]);
  } else {
    if (ref_mpi_once(ref_mpi)) printf("initial tessellation, default param\n");
  }

  RXS(ref_args_find(argc, argv, "--auto-tparams", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    auto_tparams = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi))
      printf("--auto-tparams %d requested\n", auto_tparams);
    if (auto_tparams < 0) {
      auto_tparams = REF_EGADS_ALL_TPARAM;
      if (ref_mpi_once(ref_mpi))
        printf("--auto-tparams %d set to all\n", auto_tparams);
    }
  }

  RSS(ref_geom_create(&ref_geom), "create geom");
  RSS(ref_egads_load(ref_geom, input_egads), "load");
  if (ref_mpi_once(ref_mpi) && ref_geom_effective(ref_geom))
    printf("EBody Effective Body loaded\n");
  RSS(ref_egads_quilt(ref_geom, auto_tparams, global_params), "quilt");
  RSS(ref_egads_save(ref_geom, output_egads), "save");
  RSS(ref_geom_free(ref_geom), "free geom/context");

  ref_free(global_params);
  global_params = NULL;

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) quilt_help(argv[0]);
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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

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
  REF_BOOL extrude = REF_FALSE;
  size_t end_of_string;

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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

  RXS(ref_args_find(argc, argv, "--extrude", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    extrude = REF_TRUE;
  }

  RXS(ref_args_find(argc, argv, "--planes", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && extrude) {
    if (ref_mpi_once(ref_mpi)) printf("--extrude and --planes exclusive\n");
    goto shutdown;
  }
  if (REF_EMPTY != pos) {
    REF_GRID twod_grid = ref_grid;
    REF_INT n_planes;
    if (pos + 1 >= argc) {
      if (ref_mpi_once(ref_mpi)) printf("--planes missing N\n");
      goto shutdown;
    }
    n_planes = atoi(argv[pos + 1]);
    if (n_planes < 2) {
      if (ref_mpi_once(ref_mpi))
        printf("--planes %d must be 2 or more\n", n_planes);
      goto shutdown;
    }
    if (ref_mpi_once(ref_mpi))
      printf("extrude %d layers of prisms\n", n_planes);
    RSS(ref_grid_extrude_twod(&ref_grid, twod_grid, n_planes), "extrude");
    RSS(ref_grid_free(twod_grid), "free");
  } else {
    end_of_string = strlen(out_file);
    if (ref_grid_twod(ref_grid) && (end_of_string >= 6) &&
        (strncmp(&out_file[end_of_string - 6], ".ugrid", 6)) == 0) {
      extrude = REF_TRUE;
      if (ref_mpi_once(ref_mpi))
        printf("  --extrude implicitly added to ugrid output of 2D input.\n");
    }
  }

  if (extrude) {
    REF_GRID twod_grid = ref_grid;
    if (ref_mpi_once(ref_mpi)) printf("extrude prisms\n");
    RSS(ref_grid_extrude_twod(&ref_grid, twod_grid, 2), "extrude");
    RSS(ref_grid_free(twod_grid), "free");
  }

  RXS(ref_args_find(argc, argv, "--zero-y-face", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_DBL deviation, total_deviation;
    REF_CELL ref_cell;
    REF_NODE ref_node = ref_grid_node(ref_grid);
    REF_INT faceid, group, cell, node, nodes[REF_CELL_MAX_SIZE_PER];
    if (pos + 1 >= argc) {
      if (ref_mpi_once(ref_mpi)) printf("--zero-y-face missing faceid\n");
      goto shutdown;
    }
    faceid = atoi(argv[pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf("zero y of face %d\n", faceid);
    deviation = 0.0;
    each_ref_grid_2d_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (faceid == nodes[ref_cell_node_per(ref_cell)]) {
          each_ref_cell_cell_node(ref_cell, node) {
            deviation =
                MAX(deviation, ABS(ref_node_xyz(ref_node, 1, nodes[node])));
            ref_node_xyz(ref_node, 1, nodes[node]) = 0.0;
          }
        }
      }
    }
    RSS(ref_mpi_max(ref_mpi, &deviation, &total_deviation, REF_DBL_TYPE),
        "mpi max");
    if (ref_mpi_once(ref_mpi)) printf("max deviation %e\n", deviation);
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

static REF_STATUS visualize(REF_MPI ref_mpi, int argc, char *argv[]) {
  char *in_mesh;
  char *in_sol;
  char *out_sol;
  REF_GRID ref_grid = NULL;
  REF_INT ldim;
  REF_DBL *field;
  REF_INT pos;

  if (argc < 5) goto shutdown;
  in_mesh = argv[2];
  in_sol = argv[3];
  out_sol = argv[4];

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
  if (ref_mpi_once(ref_mpi))
    printf("  read " REF_GLOB_FMT " vertices\n",
           ref_node_n_global(ref_grid_node(ref_grid)));

  if (ref_mpi_once(ref_mpi)) printf("read solution %s\n", in_sol);
  RSS(ref_part_scalar(ref_grid, &ldim, &field, in_sol), "scalar");
  if (ref_mpi_once(ref_mpi)) printf("  with leading dimension %d\n", ldim);
  ref_mpi_stopwatch_stop(ref_mpi, "read solution");

  RXS(ref_args_find(argc, argv, "--boom", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos + 4 < argc) {
    REF_INT node, i;
    REF_DBL center[3], aoa, phi, h;
    REF_DBL *dp_pinf;
    FILE *file;
    const char *vars[] = {"dp/pinf"};
    ref_malloc(dp_pinf, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_INT pressure_index = 4;
      REF_DBL gamma = 1.4;
      dp_pinf[node] =
          (field[pressure_index + ldim * node] - 1.0 / gamma) * gamma;
    }
    center[0] = atof(argv[pos + 1]);
    center[1] = atof(argv[pos + 2]);
    center[2] = atof(argv[pos + 3]);
    aoa = atof(argv[pos + 4]);
    if (ref_mpi_once(ref_mpi))
      printf("  center %f %f %f\n", center[0], center[1], center[2]);
    if (ref_mpi_once(ref_mpi)) printf("  angle of attack %f\n", aoa);
    file = NULL;
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_iso_boom_header(&file, 1, vars, out_sol), "boom header");
    }
    if (ref_mpi_once(ref_mpi)) printf(" open %s\n", out_sol);
    for (i = pos + 5; i + 1 < argc; i += 2) {
      phi = atof(argv[i]);
      h = atof(argv[i + 1]);
      if (ref_mpi_once(ref_mpi)) printf("   phi %f h %f\n", phi, h);
      RSS(ref_iso_boom_zone(file, ref_grid, dp_pinf, 1, center, aoa, phi, h),
          " boom zone");
      ref_mpi_stopwatch_stop(ref_mpi, "export ray");
    }
    ref_free(dp_pinf);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free grid");
    return REF_SUCCESS;
  }

  RXS(ref_args_find(argc, argv, "--subtract", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos && pos < argc - 1) {
    char *in_diff;
    REF_INT diff_ldim;
    REF_DBL *diff_field;
    REF_INT node, i;
    in_diff = argv[pos + 1];
    if (ref_mpi_once(ref_mpi)) printf("read diff solution %s\n", in_diff);
    RSS(ref_part_scalar(ref_grid, &diff_ldim, &diff_field, in_diff), "diff");
    ref_mpi_stopwatch_stop(ref_mpi, "read diff solution");
    REIS(ldim, diff_ldim, "difference field must have same leading dimension");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < ldim; i++) {
        field[i + ldim * node] -= diff_field[i + ldim * node];
      }
    }
    ref_free(diff_field);
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "diff field");
    for (i = 0; i < ldim; i++) {
      REF_DBL max_diff = 0.0;
      REF_DBL master_diff = 0.0;
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        max_diff = MAX(max_diff, ABS(field[i + ldim * node]));
      }
      RSS(ref_mpi_max(ref_mpi, &max_diff, &master_diff, REF_DBL_TYPE),
          "mpi max");
      if (ref_mpi_once(ref_mpi)) printf("%d max diff %e\n", i, max_diff);
    }
  }

  RXS(ref_args_find(argc, argv, "--overfun", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) {
    REF_DBL *overflow;
    REF_INT node, i, ldim_overflow;
    ldim_overflow = ldim;
    ref_malloc(overflow, ldim_overflow * ref_node_max(ref_grid_node(ref_grid)),
               REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < ldim_overflow; i++) {
        overflow[i + ldim_overflow * node] = field[i + ldim_overflow * node];
      }
    }
    ldim = ldim_overflow - 1;
    ref_free(field);
    ref_malloc(field, (ldim)*ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, u, v, w, e_0, gamma, e_i, p;

      rho = overflow[0 + ldim_overflow * node];
      u = overflow[1 + ldim_overflow * node] / rho;
      v = overflow[2 + ldim_overflow * node] / rho;
      w = overflow[3 + ldim_overflow * node] / rho;
      e_0 = overflow[4 + ldim_overflow * node] / rho;
      gamma = overflow[5 + ldim_overflow * node];
      e_i = e_0 - 0.5 * (u * u + v * v + w * w);
      p = (gamma - 1.0) * rho * e_i;

      field[0 + ldim * node] = rho;
      field[1 + ldim * node] = u;
      field[2 + ldim * node] = v;
      field[3 + ldim * node] = w;
      field[4 + ldim * node] = p;

      for (i = 5; i < ldim; i++) {
        field[i + ldim * node] = overflow[(i + 1) + ldim_overflow * node];
      }
    }
    ref_free(overflow);
  }

  RXS(ref_args_find(argc, argv, "--iso", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos < argc - 2) {
    REF_DBL *scalar;
    REF_DBL threshold;
    REF_GRID iso_grid;
    REF_INT var;
    REF_INT node;
    var = atoi(argv[pos + 1]);
    threshold = atof(argv[pos + 2]);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = field[var + ldim * node] - threshold;
    }
    RSS(ref_iso_insert(&iso_grid, ref_grid, scalar), "iso");
    if (ref_mpi_once(ref_mpi))
      printf("write isosurface geometry %s\n", out_sol);
    RSS(ref_gather_by_extension(iso_grid, out_sol), "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "write isosurface geometry");

    ref_grid_free(iso_grid);
    ref_free(scalar);
  } else {
    if (ref_mpi_once(ref_mpi))
      printf("write %d ldim solution %s\n", ldim, out_sol);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, field, NULL, out_sol),
        "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "write solution");
  }
  ref_free(field);
  RSS(ref_grid_free(ref_grid), "free grid");

  return REF_SUCCESS;
shutdown:
  if (ref_mpi_once(ref_mpi)) visualize_help(argv[0]);
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
    printf("refine %s, on or after 1.9.3\n", VERSION);
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
  } else if (strncmp(argv[1], "d", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(distance(ref_mpi, argc, argv), "distance");
    } else {
      if (ref_mpi_once(ref_mpi)) distance_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "e", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(examine(ref_mpi, argc, argv), "examine");
    } else {
      if (ref_mpi_once(ref_mpi)) examine_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "g", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(grow(ref_mpi, argc, argv), "grow");
    } else {
      if (ref_mpi_once(ref_mpi)) grow_help(argv[0]);
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
  } else if (strncmp(argv[1], "n", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(node(ref_mpi, argc, argv), "translate");
    } else {
      if (ref_mpi_once(ref_mpi)) node_help(argv[0]);
      goto shutdown;
    }
  } else if (strncmp(argv[1], "q", 1) == 0) {
    if (REF_EMPTY == help_pos) {
      RSS(quilt(ref_mpi, argc, argv), "quilt");
    } else {
      if (ref_mpi_once(ref_mpi)) quilt_help(argv[0]);
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
      RSS(visualize(ref_mpi, argc, argv), "visualize");
    } else {
      if (ref_mpi_once(ref_mpi)) visualize_help(argv[0]);
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
