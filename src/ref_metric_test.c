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

#include "ref_metric.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_histogram.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_validation.h"

/*
./test.sh ref_metric && ./ref_metric_test \
  ../acceptance/2d/linear/two/accept-2d-two-09.b8.ugrid \
  --parent ../acceptance/2d/linear/two/accept-2d-two-08.b8.ugrid \
  ../acceptance/2d/linear/two/accept-2d-two-08.metric
*/

int main(int argc, char *argv[]) {
  REF_INT fixed_point_pos = REF_EMPTY;
  REF_INT curve_limit_pos = REF_EMPTY;
  REF_INT parent_pos = REF_EMPTY;
  REF_INT wlp_pos = REF_EMPTY;
  REF_INT moving_pos = REF_EMPTY;
  REF_INT explore_pos = REF_EMPTY;
  REF_INT lp_pos = REF_EMPTY;
  REF_INT opt_goal_pos = REF_EMPTY;
  REF_INT no_goal_pos = REF_EMPTY;
  REF_INT venditti_pos = REF_EMPTY;
  REF_INT belme_pos = REF_EMPTY;
  REF_INT euler_opt_goal_pos = REF_EMPTY;
  REF_INT euler_cons_pos = REF_EMPTY;
  REF_INT viscous_cons_pos = REF_EMPTY;
  REF_INT hmax_pos = REF_EMPTY;
  REF_INT buffer_pos = REF_EMPTY;
  REF_INT kexact_pos = REF_EMPTY;
  REF_INT complexity_pos = REF_EMPTY;
  REF_INT gradation_pos = REF_EMPTY;
  REF_INT cloud_pos = REF_EMPTY;
  REF_INT wake_pos = REF_EMPTY;

  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");
  ref_mpi_stopwatch_start(ref_mpi);

  RXS(ref_args_find(argc, argv, "--fixed-point", &fixed_point_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--curve-limit", &curve_limit_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--parent", &parent_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--wlp", &wlp_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--lp", &lp_pos), REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--moving", &moving_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--explore", &explore_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--opt-goal", &opt_goal_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--no-goal", &no_goal_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--venditti", &venditti_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--belme", &belme_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--euler-opt-goal", &euler_opt_goal_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--euler-cons", &euler_cons_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--viscous-cons", &viscous_cons_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--kexact", &kexact_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--hmax", &hmax_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--buffer", &buffer_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--complexity", &complexity_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--gradation", &gradation_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--cloud", &cloud_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--wake", &wake_pos), REF_NOT_FOUND,
      "arg search");

  if (curve_limit_pos != REF_EMPTY) {
    REF_GRID ref_grid;

    REIS(1, curve_limit_pos,
         "required args: --curve-limit grid.ext input.metric geom.egads "
         "[assoc.gas]");
    REIS(5, argc,
         "required args: --curve-limit grid.ext input.metric geom.egads "
         "[assoc.gas]");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 1");
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load parent metric in position 2");
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[4]),
        "unable to load egads in position 3");

    RSS(ref_metric_constrain_curvature(ref_grid), "crv const");
    RSS(ref_gather_metric(ref_grid, "ref_metric_test_curve_limit.metric"),
        "export curve limit metric");

    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_metric_test_curve_limit"),
        "al");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (buffer_pos == 1) {
    REF_GRID ref_grid;
    REF_DBL complexity, *metric;
    REIS(6, argc,
         "required args: --buffer grid.ext input-metric.solb complexity "
         "output-metric.solb");

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    if (ref_mpi_once(ref_mpi)) printf("reading metric %s\n", argv[3]);
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load parent metric in position 3");
    ref_mpi_stopwatch_stop(ref_mpi, "read metric");

    complexity = atof(argv[4]);
    if (ref_mpi_once(ref_mpi))
      printf("buffering at complexity %f\n", complexity);

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_from_node(metric, ref_grid_node(ref_grid)), "set node");
    RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
        "buffer at complexity");
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_mpi_stopwatch_stop(ref_mpi, "buffer");

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[5]);
    RSS(ref_gather_metric(ref_grid, argv[5]), "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (wlp_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, *weight, *metric;
    REF_INT p;
    REF_DBL gradation, complexity, current_complexity, hmin, hmax;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim, wdim;
    REIS(1, wlp_pos,
         "required args: --wlp grid.meshb scalar.solb weight.solb "
         "p gradation complexity output-metric.solb");
    if (9 > argc) {
      printf(
          "required args: --wlp grid.meshb scalar.solb weight.solb "
          "p gradation complexity output-metric.solb\n");
      return REF_FAILURE;
    }
    hmin = -1.0;
    hmax = -1.0;
    if (REF_EMPTY != hmax_pos) {
      if (hmax_pos >= argc - 1) {
        printf("option missing value: --hmax max_edge_length\n");
        return REF_FAILURE;
      }
      hmax = atof(argv[hmax_pos + 1]);
    }

    p = atoi(argv[5]);
    gradation = atof(argv[6]);
    complexity = atof(argv[7]);
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("hmin %f hmax %f (negative is inactive)\n", hmin, hmax);
      printf("buffer %d (negative is inactive)\n", buffer_pos);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    if (ref_mpi_para(ref_mpi)) {
      if (ref_mpi_once(ref_mpi)) printf("part %s\n", argv[2]);
      RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]), "part");
      ref_mpi_stopwatch_stop(ref_mpi, "part mesh");
    } else {
      if (ref_mpi_once(ref_mpi)) printf("import %s\n", argv[2]);
      RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "import");
      ref_mpi_stopwatch_stop(ref_mpi, "import mesh");
    }

    if (ref_mpi_once(ref_mpi)) printf("reading scalar %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load scalar in position 3");
    REIS(1, ldim, "expected one scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "read scalar");

    if (ref_mpi_once(ref_mpi)) printf("reading weight %s\n", argv[4]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &wdim, &weight, argv[4]),
        "unable to load weight in position 4");
    REIS(1, wdim, "expected one weight");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_lp(metric, ref_grid, scalar, weight, reconstruction, p,
                      gradation, complexity),
        "lp norm");
    ref_mpi_stopwatch_stop(ref_mpi, "compute metric");
    if (REF_EMPTY != buffer_pos) {
      RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
          "buffer at complexity");
    }
    if (hmin > 0.0 || hmax > 0.0) {
      RSS(ref_metric_limit_h_at_complexity(metric, ref_grid, hmin, hmax,
                                           complexity),
          "limit at complexity");
    }
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_mpi))
      printf("actual complexity %e\n", current_complexity);
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(weight);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[8]);
    RSS(ref_gather_metric(ref_grid, argv[8]), "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (lp_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, *metric;
    REF_INT p;
    REF_DBL gradation, complexity, current_complexity, hmin, hmax;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REIS(1, lp_pos,
         "required args: --lp grid.meshb scalar-mach.solb p gradation "
         "complexity output-metric.solb");
    if (8 > argc) {
      printf(
          "required args: --lp grid.meshb scalar-mach.solb p gradation "
          "complexity output-metric.solb\n");
      return REF_FAILURE;
    }
    hmin = -1.0;
    hmax = -1.0;
    if (REF_EMPTY != hmax_pos) {
      if (hmax_pos >= argc - 1) {
        printf("option missing value: --hmax max_edge_length\n");
        return REF_FAILURE;
      }
      hmax = atof(argv[hmax_pos + 1]);
    }

    p = atoi(argv[4]);
    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("hmin %f hmax %f (negative is inactive)\n", hmin, hmax);
      printf("buffer %d (negative is inactive)\n", buffer_pos);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    if (ref_mpi_once(ref_mpi)) printf("reading scalar %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load scalar in position 3");
    REIS(1, ldim, "expected one scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "read scalar");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                      gradation, complexity),
        "lp norm");
    ref_mpi_stopwatch_stop(ref_mpi, "compute metric");
    if (REF_EMPTY != buffer_pos) {
      RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
          "buffer at complexity");
    }
    if (hmin > 0.0 || hmax > 0.0) {
      RSS(ref_metric_limit_h_at_complexity(metric, ref_grid, hmin, hmax,
                                           complexity),
          "limit at complexity");
    }
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_mpi))
      printf("actual complexity %e\n", current_complexity);
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[7]);
    RSS(ref_gather_metric(ref_grid, argv[7]), "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (moving_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *displaced, *scalar, *metric;
    REF_INT p;
    REF_DBL gradation, complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REIS(1, moving_pos,
         "required args: --moving grid.meshb displaced.solb scalar.solb p "
         "gradation "
         "complexity output-metric.solb");
    if (9 > argc) {
      printf(
          "required args: --moving grid.meshb displaced.solb scalar.solb p "
          "gradation "
          "complexity output-metric.solb\n");
      return REF_FAILURE;
    }

    p = atoi(argv[5]);
    gradation = atof(argv[6]);
    complexity = atof(argv[7]);
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    if (ref_mpi_once(ref_mpi)) printf("reading displaced %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &displaced, argv[3]),
        "unable to load dispaced in position 3");
    REIS(3, ldim, "expected 3 [x,y,z]");
    ref_mpi_stopwatch_stop(ref_mpi, "read scalar");

    if (ref_mpi_once(ref_mpi)) printf("reading scalar %s\n", argv[4]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[4]),
        "unable to load scalar in position 4");
    REIS(1, ldim, "expected 1 scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "read scalar");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_moving_multiscale(metric, ref_grid, displaced, scalar,
                                     reconstruction, p, gradation, complexity),
        "moving multiscale norm");

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(displaced);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[8]);
    RSS(ref_gather_metric(ref_grid, argv[8]), "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (explore_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *field, *scalar, *metric, *output;
    REF_INT p;
    REF_DBL gradation, complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim, node, var;
    REF_DBL h0, multiscale_system[12];

    REIS(1, explore_pos,
         "required args: --explore grid.meshb scalars.solb p gradation "
         "complexity metric-h.tec");
    if (8 > argc) {
      printf(
          "required args: --explore grid.meshb scalars.solb p gradation "
          "complexity metric-h.tec\n");
      return REF_FAILURE;
    }

    p = atoi(argv[4]);
    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("buffer %d (negative is inactive)\n", buffer_pos);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    if (ref_mpi_para(ref_mpi)) {
      if (ref_mpi_once(ref_mpi)) printf("part %s\n", argv[2]);
      RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]), "part");
      ref_mpi_stopwatch_stop(ref_mpi, "part mesh");
    } else {
      if (ref_mpi_once(ref_mpi)) printf("import %s\n", argv[2]);
      RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "import");
      ref_mpi_stopwatch_stop(ref_mpi, "import mesh");
    }

    if (ref_mpi_once(ref_mpi))
      printf("reading field with scalars %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load scalar in position 3");
    RAS(ldim > 0, "expected at least one scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "read scalar");

    ref_malloc(output, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    for (var = 0; var < ldim; var++) {
      if (ref_mpi_once(ref_mpi)) printf("scalar %d of %d\n", var, ldim);
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        scalar[node] = field[node + ldim * var];
      }
      RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                        gradation, complexity),
          "lp norm");
      ref_mpi_stopwatch_stop(ref_mpi, "compute metric");
      if (REF_EMPTY != buffer_pos) {
        RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
            "buffer at complexity");
      }
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        RSS(ref_matrix_diag_m(&(metric[6 * node]), multiscale_system),
            "decomp");
        RSS(ref_matrix_ascending_eig(multiscale_system), "sort eig");
        if (multiscale_system[0] < 0.0) RSS(REF_DIV_ZERO, "sqrt(-1)");
        h0 = sqrt(multiscale_system[0]);
        if (!ref_math_divisible(1.0, h0)) RSS(REF_DIV_ZERO, "inf h0");
        output[node + var * ldim] = 1.0 / h0;
      }
    }
    ref_free(metric);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing sizes %s\n", argv[7]);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, output, NULL, argv[7]),
        "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");
    ref_free(output);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (fixed_point_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, *hess, *metric;
    REF_INT p, n, timestep, timestep_increment, node, im;
    REF_DBL gradation, complexity, current_complexity, hmin, hmax;
    REF_RECON_RECONSTRUCTION reconstruction;
    char solb[1024];
    REF_INT ldim;
    REIS(1, fixed_point_pos,
         "required args: --fixed-point grid.meshb scalar-mach-root Ntimesteps "
         "timestep_increment p gradation complexity output-metric.solb");
    if (10 > argc) {
      printf(
          "required args: --fixed-point grid.meshb scalar-mach-root Ntimesteps "
          "timestep_increment p gradation complexity output-metric.solb");
      return REF_FAILURE;
    }
    hmin = -1.0;
    hmax = -1.0;
    if (REF_EMPTY != hmax_pos) {
      if (hmax_pos >= argc - 1) {
        printf("option missing value: --hmax max_edge_length\n");
        return REF_FAILURE;
      }
      hmax = atof(argv[hmax_pos + 1]);
    }

    n = atoi(argv[4]);
    timestep_increment = atoi(argv[5]);
    p = atoi(argv[6]);
    gradation = atof(argv[7]);
    complexity = atof(argv[8]);
    reconstruction = REF_RECON_KEXACT;

    if (ref_mpi_once(ref_mpi)) {
      printf("N=%d\n", n);
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("hmin %f hmax %f (negative is inactive)\n", hmin, hmax);
      printf("buffer %d (negative is inactive)\n", buffer_pos);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc(hess, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    for (timestep = timestep_increment; timestep <= n;
         timestep += timestep_increment) {
      snprintf(solb, 1024, "%s%d.solb", argv[3], timestep);
      if (ref_mpi_once(ref_mpi))
        printf("reading and reconstructing hessian for  %s\n", solb);
      RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, solb),
          "unable to load scalar in position 3");
      REIS(1, ldim, "expected one scalar");
      RSS(ref_recon_hessian(ref_grid, scalar, hess, reconstruction), "hess");
      ref_free(scalar);
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        for (im = 0; im < 6; im++) {
          metric[im + 6 * node] += hess[im + 6 * node];
        }
      }
      ref_mpi_stopwatch_stop(ref_mpi, "timestep processed");
    }

    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p),
        "local lp norm scaling");
    ref_mpi_stopwatch_stop(ref_mpi, "local scale metric");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");
    ref_mpi_stopwatch_stop(ref_mpi, "metric gradation and complexity");

    if (REF_EMPTY != buffer_pos) {
      RSS(ref_metric_buffer_at_complexity(metric, ref_grid, complexity),
          "buffer at complexity");
      ref_mpi_stopwatch_stop(ref_mpi, "buffer metric");
    }
    if (hmin > 0.0 || hmax > 0.0) {
      RSS(ref_metric_limit_h_at_complexity(metric, ref_grid, hmin, hmax,
                                           complexity),
          "limit at complexity");
      ref_mpi_stopwatch_stop(ref_mpi, "h-limit metric");
    }
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("actual complexity %e\n", current_complexity);
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(hess);
    ref_free(metric);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[9]);
    RSS(ref_gather_metric(ref_grid, argv[9]), "export curve limit metric");
    ref_mpi_stopwatch_stop(ref_mpi, "write metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (opt_goal_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, *metric;
    REF_INT p;
    REF_DBL gradation, complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_DBL current_complexity, hmin, hmax;

    REIS(1, opt_goal_pos,
         "required args: --opt-goal grid.meshb solution.solb p "
         "gradation complexity output-metric.solb");
    if (8 > argc) {
      printf(
          "required args: --opt-goal grid.meshb solution.solb p "
          "gradation complexity output-metric.solb\n");
      return REF_FAILURE;
    }
    hmin = -1.0;
    hmax = -1.0;
    if (REF_EMPTY != hmax_pos) {
      if (hmax_pos >= argc - 1) {
        printf("option missing value: --hmax max_edge_length\n");
        return REF_FAILURE;
      }
      hmax = atof(argv[hmax_pos + 1]);
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    p = atoi(argv[4]);
    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("hmin %f hmax %f (negative is inactive)\n", hmin, hmax);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading solution %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load scalar in position 3");
    REIS(20, ldim, "expected 20 (5*adj,5*xflux,5*yflux,5*zflux) scalar");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_opt_goal(metric, ref_grid, 5, scalar, reconstruction, p,
                            gradation, complexity),
        "opt goal");
    if (hmin > 0.0 || hmax > 0.0) {
      RSS(ref_metric_limit_h_at_complexity(metric, ref_grid, hmin, hmax,
                                           complexity),
          "limit at complexity");
    }
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_mpi))
      printf("actual complexity %e\n", current_complexity);

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[7]);
    RSS(ref_gather_metric(ref_grid, argv[7]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (no_goal_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_DBL *scalar, *metric;
    REF_INT p;
    REF_DBL gradation, complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_DBL current_complexity, hmin, hmax;
    REF_INT i, node;

    REIS(1, no_goal_pos,
         "required args: --no-goal grid.meshb solution.solb complexity p "
         "gradation output-metric.solb");
    if (8 > argc) {
      printf(
          "required args: --no-goal grid.meshb solution.solb complexity p "
          "gradation output-metric.solb\n");
      return REF_FAILURE;
    }
    hmin = -1.0;
    hmax = -1.0;
    if (REF_EMPTY != hmax_pos) {
      if (hmax_pos >= argc - 1) {
        printf("option missing value: --hmax max_edge_length\n");
        return REF_FAILURE;
      }
      hmax = atof(argv[hmax_pos + 1]);
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    p = atoi(argv[4]);
    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    if (ref_mpi_once(ref_mpi)) {
      printf("Lp=%d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
      printf("hmin %f hmax %f (negative is inactive)\n", hmin, hmax);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_node = ref_grid_node(ref_grid);

    if (ref_mpi_once(ref_mpi)) printf("reading solution %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load scalar in position 3");
    REIS(20, ldim, "expected 20 (5*adj,5*xflux,5*yflux,5*zflux) scalar");

    /* linear function evaluates to unit adjoint weights */
    each_ref_node_valid_node(ref_node, node) {
      for (i = 0; i < 5; i++) {
        scalar[i + 20 * node] = ref_node_xyz(ref_node, 0, node) +
                                ref_node_xyz(ref_node, 1, node) +
                                ref_node_xyz(ref_node, 2, node);
      }
    }

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_opt_goal(metric, ref_grid, 5, scalar, reconstruction, p,
                            gradation, complexity),
        "opt goal");
    if (hmin > 0.0 || hmax > 0.0) {
      RSS(ref_metric_limit_h_at_complexity(metric, ref_grid, hmin, hmax,
                                           complexity),
          "limit at complexity");
    }
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_mpi))
      printf("actual complexity %e\n", current_complexity);

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[7]);
    RSS(ref_gather_metric(ref_grid, argv[7]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (venditti_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_DBL *weight, *scalar, *metric, *implied;
    REF_INT p = 2;
    REF_DBL gradation, complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_DBL h, h0, h_h0, scale, h_ms;
    REF_INT i, node;
    REF_DBL multiscale_system[12];
    REF_DBL *system;
    REF_INT nsystem;

    REIS(1, venditti_pos,
         "required args: --venditti grid.meshb scalar.solb weight.solb "
         "gradation complexity output-metric.solb");
    if (8 > argc) {
      printf(
          "required args: --venditti grid.meshb scalar.solb weight.solb "
          "gradation complexity output-metric.solb");
      return REF_FAILURE;
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    if (ref_mpi_once(ref_mpi)) {
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_node = ref_grid_node(ref_grid);

    if (ref_mpi_once(ref_mpi)) printf("reading scalar %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load scalar in position 3");
    REIS(1, ldim, "expected one scalar");

    if (ref_mpi_once(ref_mpi)) printf("reading weight %s\n", argv[4]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &weight, argv[4]),
        "unable to load scalar in position 4");
    REIS(1, ldim, "expected one weight");

    if (ref_mpi_once(ref_mpi)) printf("multiscale metric\n");
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, reconstruction, p,
                      gradation, complexity),
        "lp");

    if (ref_mpi_once(ref_mpi)) printf("imply current metric\n");
    ref_malloc(implied, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_imply_from(implied, ref_grid), "imply");

    nsystem = 6;
    ref_malloc_init(system, nsystem * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 0.0);

    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_matrix_diag_m(&(metric[6 * node]), multiscale_system), "decomp");
      if (ref_grid_twod(ref_grid)) {
        RSS(ref_matrix_ascending_eig_twod(multiscale_system), "sort eig twod");
      } else {
        RSS(ref_matrix_ascending_eig(multiscale_system), "sort eig");
      }
      h0 = ref_matrix_sqrt_vt_m_v(&(implied[6 * node]),
                                  &(ref_matrix_vec(multiscale_system, 0, 0)));
      if (!ref_math_divisible(1.0, h0)) RSS(REF_DIV_ZERO, "inf h0");
      h0 = 1.0 / h0;
      h_h0 = weight[node];
      h_h0 = MAX(0.1, MIN(10.0, h_h0));
      h = h_h0 * h0;
      h_ms = ref_matrix_eig(multiscale_system, 0);
      if (!ref_math_divisible(1.0, sqrt(h_ms))) RSS(REF_DIV_ZERO, "inf h_ms");
      h_ms = 1.0 / sqrt(h_ms);
      if (!ref_math_divisible((h_ms * h_ms), (h * h)))
        RSS(REF_DIV_ZERO, "inf scale");
      scale = (h_ms * h_ms) / (h * h);
      for (i = 0; i < 6; i++) metric[i + 6 * node] *= scale;
      system[0 + nsystem * node] = h0;
      system[1 + nsystem * node] = h_h0;
      system[2 + nsystem * node] = h;
      system[3 + nsystem * node] = h_ms;
      system[4 + nsystem * node] = scale;
      RSS(ref_matrix_diag_m(&(metric[6 * node]), multiscale_system), "decomp");
      RSS(ref_matrix_ascending_eig(multiscale_system), "sort eig");
      h_ms = ref_matrix_eig(multiscale_system, 0);
      if (!ref_math_divisible(1.0, sqrt(h_ms))) RSS(REF_DIV_ZERO, "post h_ms");
      h_ms = 1.0 / sqrt(h_ms);
      system[5 + nsystem * node] = h_ms / h;
    }

    if (ref_mpi_once(ref_mpi))
      printf("global scaling and gradation limiting\n");

    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation at complexity");

    if (ref_mpi_once(ref_mpi))
      printf("writing res,dual,weight ref_vend_system.tec\n");
    RSS(ref_gather_scalar_by_extension(ref_grid, nsystem, system, NULL,
                                       "ref_vend_system.tec"),
        "export primitive_dual");
    ref_free(system);

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(implied);
    ref_free(metric);
    ref_free(weight);
    ref_free(scalar);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[7]);
    RSS(ref_gather_metric(ref_grid, argv[7]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (belme_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *prim_dual, *metric;
    REF_DBL gradation = -1.0;
    REF_DBL mach, re, temperature;
    REF_DBL complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_INT p = 1;

    REIS(1, belme_pos,
         "required args: --belme grid.meshb prim_dual.solb "
         "Mach Re Temperature(Kelvin) "
         "complexity output-metric.solb");
    if (9 > argc) {
      printf(
          "required args: --belme grid.meshb prim_dual.solb "
          "Mach Re Temperature(Kelvin) "
          "complexity output-metric.solb");
      return REF_FAILURE;
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    mach = atof(argv[4]);
    re = atof(argv[5]);
    temperature = atof(argv[6]);
    complexity = atof(argv[7]);
    if (ref_mpi_once(ref_mpi)) {
      printf("p-norm %d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading prim_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &prim_dual, argv[3]),
        "unable to load scalar in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected rho,u,v,w,p,5*adj "
        "or rho,u,v,w,p,turb,6*adj");

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);

    RSS(ref_metric_belme_gfe(metric, ref_grid, ldim, prim_dual, reconstruction),
        "gfe");
    RSS(ref_metric_belme_gu(metric, ref_grid, ldim, prim_dual, mach, re,
                            temperature, reconstruction),
        "gu");

    RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6),
        "update ghosts");

    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p), "local scale");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation");

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[8]);
    RSS(ref_gather_metric(ref_grid, argv[8]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (euler_opt_goal_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *prim_dual, *metric;
    REF_DBL gradation = -1.0;
    REF_DBL complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_INT p = 1;

    REIS(1, euler_opt_goal_pos,
         "required args: --euler-opt-goal grid.meshb prim_dual.solb "
         "complexity output-metric.solb");
    if (6 > argc) {
      printf(
          "required args: --euler-opt-goal grid.meshb prim_dual.solb "
          "complexity output-metric.solb");
      return REF_FAILURE;
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    complexity = atof(argv[4]);
    if (ref_mpi_once(ref_mpi)) {
      printf("p-norm %d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading prim_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &prim_dual, argv[3]),
        "unable to load scalar in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected rho,u,v,w,p,5*adj "
        "or rho,u,v,w,p,turb,6*adj");

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);

    RSS(ref_metric_belme_gfe(metric, ref_grid, ldim, prim_dual, reconstruction),
        "gfe");
    RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6),
        "update ghosts");

    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p), "local scale");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation");

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[5]);
    RSS(ref_gather_metric(ref_grid, argv[5]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (euler_cons_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *prim_dual, *g, *metric;
    REF_DBL gradation = -1.0;
    REF_DBL complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_INT p = 1;

    REIS(1, euler_cons_pos,
         "required args: --euler-cons grid.meshb prim_dual.solb "
         "gradation complexity output-metric.solb");
    if (7 > argc) {
      printf(
          "required args: --euler-cons grid.meshb prim_dual.solb "
          "gradation complexity output-metric.solb");
      return REF_FAILURE;
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    gradation = atof(argv[4]);
    complexity = atof(argv[5]);
    if (ref_mpi_once(ref_mpi)) {
      printf("p-norm %d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading prim_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &prim_dual, argv[3]),
        "unable to load scalar in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected rho,u,v,w,p,5*adj "
        "or rho,u,v,w,p,turb,6*adj");

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc_init(g, 5 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    RSS(ref_metric_cons_euler_g(g, ref_grid, ldim, prim_dual, reconstruction),
        "cons euler g weights");

    RSS(ref_metric_cons_assembly(metric, g, ref_grid, ldim, prim_dual,
                                 reconstruction),
        "cons metric assembly");
    ref_free(g);
    RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6),
        "update ghosts");

    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p), "local scale");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation");

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[6]);
    RSS(ref_gather_metric(ref_grid, argv[6]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (viscous_cons_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *prim_dual, *g, *metric;
    REF_DBL gradation = -1.0;
    REF_DBL mach, re, temperature;
    REF_DBL complexity;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_INT ldim;
    REF_INT p = 1;

    REIS(1, viscous_cons_pos,
         "required args: --viscous-cons grid.meshb prim_dual.solb "
         "Mach Re Temperature(Kelvin) "
         "gradation complexity output-metric.solb");
    if (10 > argc) {
      printf(
          "required args: --viscous-cons grid.meshb prim_dual.solb "
          "Mach Re Temperature(Kelvin) "
          "gradation complexity output-metric.solb");
      return REF_FAILURE;
    }
    if (REF_EMPTY != kexact_pos) {
      reconstruction = REF_RECON_KEXACT;
    }

    mach = atof(argv[4]);
    re = atof(argv[5]);
    temperature = atof(argv[6]);
    gradation = atof(argv[7]);
    complexity = atof(argv[8]);
    if (ref_mpi_once(ref_mpi)) {
      printf("Mach %f\n", mach);
      printf("Re %e\n", re);
      printf("Temp(K) %f\n", temperature);
      printf("p-norm %d\n", p);
      printf("gradation %f\n", gradation);
      printf("complexity %f\n", complexity);
      printf("reconstruction %d\n", (int)reconstruction);
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading prim_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &prim_dual, argv[3]),
        "unable to load scalar in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected rho,u,v,w,p,5*adj "
        "or rho,u,v,w,p,turb,6*adj");

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc_init(g, 5 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    RSS(ref_metric_cons_euler_g(g, ref_grid, ldim, prim_dual, reconstruction),
        "cons viscous g weights");
    RSS(ref_metric_cons_viscous_g(g, ref_grid, ldim, prim_dual, mach, re,
                                  temperature, reconstruction),
        "cons viscous g weights");

    RSS(ref_metric_cons_assembly(metric, g, ref_grid, ldim, prim_dual,
                                 reconstruction),
        "cons metric assembly");
    ref_free(g);
    RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6),
        "update ghosts");

    RSS(ref_metric_local_scale(metric, NULL, ref_grid, p), "local scale");
    RSS(ref_metric_gradation_at_complexity(metric, ref_grid, gradation,
                                           complexity),
        "gradation");

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_mpi)) printf("writing metric %s\n", argv[9]);
    RSS(ref_gather_metric(ref_grid, argv[9]), "export opt goal metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (complexity_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL target_complexity, current_complexity;
    REF_DBL *metric;
    REF_INT i, node;

    REIS(1, complexity_pos,
         "required args: --complexity grid.ext input-metric.solb complexity "
         "output-metric.solb");
    REIS(6, argc,
         "required args: --complexity grid.ext input-metric.solb complexity "
         "output-metric.solb");
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load grid in position 2");
    if (ref_mpi_once(ref_mpi)) printf("reading metric %s\n", argv[3]);
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load metric in position 3");
    target_complexity = atof(argv[4]);
    if (ref_mpi_once(ref_mpi))
      printf("desired complexity %e\n", target_complexity);

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_from_node(metric, ref_grid_node(ref_grid)), "get node");
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_mpi))
      printf("actual complexity %e\n", current_complexity);
    if (!ref_math_divisible(target_complexity, current_complexity)) {
      ref_free(metric);
      RSS(ref_grid_free(ref_grid), "free");
      RSS(ref_mpi_free(ref_mpi), "free");
      RSS(ref_mpi_stop(), "stop");
      return REF_DIV_ZERO;
    }
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 6; i++) {
        metric[i + 6 * node] *=
            pow(target_complexity / current_complexity, 2.0 / 3.0);
      }
    }
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("writing metric %s\n", argv[5]);
    RSS(ref_gather_metric(ref_grid, argv[5]), "export scaled metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (cloud_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_DBL h, m[6], d[12];
    REF_INT node;
    FILE *file;
    char filename[] = "ref_metric_cloud.dat";

    REIS(1, cloud_pos, "required args: --cloud grid.ext input-metric.solb");
    REIS(4, argc, "required args: --cloud grid.ext input-metric.solb");
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load grid in position 2");
    ref_node = ref_grid_node(ref_grid);
    if (ref_mpi_once(ref_mpi)) printf("reading metric %s\n", argv[3]);
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load metric in position 3");

    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RSS(ref_node_metric_get(ref_node, node, m), "get");
      RSS(ref_matrix_diag_m(m, d), "diag");
      h = MAX(ref_matrix_eig(d, 0), ref_matrix_eig(d, 1));
      h = MAX(h, ref_matrix_eig(d, 2));
      h = 1.0 / sqrt(h);
      fprintf(file, "%e %e %e %e\n", ref_node_xyz(ref_node, 0, node),
              ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
              h);
    }

    fclose(file);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (gradation_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL complexity, gradation, t;
    REF_DBL *metric;
    REF_INT npass = 20, pass;
    char *gradation_type;

    REIS(1, gradation_pos,
         "required args: --gradation grid.ext input-metric.solb "
         "output-metric.solb  {metric beta|mixed beta t}");

    if (7 > argc) {
      printf(
          "required args: --gradation grid.ext input-metric.solb "
          "output-metric.solb {metric beta|mixed beta t}");
      return REF_FAILURE;
    }
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load grid in position 2");
    if (ref_mpi_once(ref_mpi)) printf("reading metric %s\n", argv[3]);
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load metric in position 3");
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_from_node(metric, ref_grid_node(ref_grid)), "get node");
    gradation_type = argv[5];

    if (ref_mpi_once(ref_mpi)) printf("gradation type %s\n", gradation_type);
    if (strcmp(gradation_type, "metric") == 0) {
      REIS(7, argc,
           "required args: --gradation grid.ext input-metric.solb "
           "output-metric.solb metric beta");
      gradation = atof(argv[6]);
      if (ref_mpi_once(ref_mpi))
        printf("metric-space gradation %e\n", gradation);
      for (pass = 0; pass < npass; pass++) {
        RSS(ref_metric_complexity(metric, ref_grid, &complexity), "cmp");
        if (ref_mpi_once(ref_mpi))
          printf("pass %d complexity %.5e\n", pass, complexity);
        RSS(ref_metric_metric_space_gradation(metric, ref_grid, gradation),
            "metric_space");
      }
      RSS(ref_metric_complexity(metric, ref_grid, &complexity), "cmp");
      if (ref_mpi_once(ref_mpi))
        printf("pass %d complexity %.5e\n", npass, complexity);
    } else if (strcmp(gradation_type, "mixed") == 0) {
      REIS(8, argc,
           "required args: --gradation grid.ext input-metric.solb "
           "output-metric.solb mixed beta t");
      gradation = atof(argv[6]);
      t = atof(argv[7]);
      if (ref_mpi_once(ref_mpi))
        printf("mixed-space gradation %e %e\n", gradation, t);
      for (pass = 0; pass < npass; pass++) {
        RSS(ref_metric_complexity(metric, ref_grid, &complexity), "cmp");
        if (ref_mpi_once(ref_mpi))
          printf("pass %d complexity %.5e\n", pass, complexity);
        RSS(ref_metric_mixed_space_gradation(metric, ref_grid, gradation, t),
            "metric_space");
      }
      RSS(ref_metric_complexity(metric, ref_grid, &complexity), "cmp");
      if (ref_mpi_once(ref_mpi))
        printf("pass %d complexity %.5e\n", npass, complexity);
    } else {
      printf("%s: %d: %s %s\n", __FILE__, __LINE__, "unknown gradation",
             gradation_type);
      return REF_NOT_FOUND;
    }

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("writing metric %s\n", argv[4]);
    RSS(ref_gather_metric(ref_grid, argv[4]), "export scaled metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (parent_pos != REF_EMPTY) {
    REF_GRID ref_grid, parent_grid;
    REF_INTERP ref_interp;

    REIS(2, parent_pos,
         "required args: grid.ext --parent pgrid.ext pgrid.metric");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]),
        "unable to load target grid in position 1");
    RSS(ref_import_by_extension(&parent_grid, ref_mpi, argv[3]),
        "unable to load parent grid in position 3");
    RSS(ref_part_metric(ref_grid_node(parent_grid), argv[4]),
        "unable to load parent grid in position 4");

    RSS(ref_interp_create(&ref_interp, ref_grid, parent_grid), "map");
    RSS(ref_interp_locate(ref_interp), "map");
    RSS(ref_metric_interpolate(ref_interp), "interp");

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "qual");
    RSS(ref_histogram_ratio(ref_grid), "rat");

    RSS(ref_gather_metric(ref_grid, "ref_metric_interp.solb"),
        "unable to gather metric");

    RSS(ref_interp_free(ref_interp), "free");
    RSS(ref_grid_free(parent_grid), "free");
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (wake_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_DBL *dist, *field, *metric, m[6], m0[6];
    REF_INT ldim, i, node, gradation;

    REIS(1, wake_pos,
         "required args: --wake grid.ext distance.solb volume.solb "
         "metric.solb");
    REIS(6, argc,
         "required args: --wake grid.ext distance.solb volume.solb "
         "metric.solb");
    if (ref_mpi_once(ref_mpi)) printf("part grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to part grid in position 2");
    ref_node = ref_grid_node(ref_grid);
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");
    if (ref_mpi_once(ref_mpi)) printf("reading distance %s\n", argv[3]);
    RSS(ref_part_scalar(ref_node, &ldim, &dist, argv[3]),
        "unable to load distance in position 3");
    if (ref_mpi_once(ref_mpi)) printf("distance ldim %d\n", ldim);
    ref_mpi_stopwatch_stop(ref_mpi, "read dist");
    REIS(1, ldim, "expect [distance]");
    if (ref_mpi_once(ref_mpi)) printf("reading solution %s\n", argv[4]);
    RSS(ref_part_scalar(ref_node, &ldim, &field, argv[4]),
        "unable to load solution in position 4");
    if (ref_mpi_once(ref_mpi)) printf("ldim %d\n", ldim);
    ref_mpi_stopwatch_stop(ref_mpi, "read vol");
    REIS(6, ldim, "expect [rho,u,v,w,p,turb1]");

    if (ref_mpi_once(ref_mpi)) printf("imply current metric\n");
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_imply_from(metric, ref_grid), "imply");
    ref_mpi_stopwatch_stop(ref_mpi, "imply");

    if (ref_grid_twod(ref_grid)) {
      each_ref_node_valid_node(ref_node, node) {
        REF_DBL x0 = -0.1;
        REF_DBL x1 = 2;
        REF_DBL y0 = -2;
        REF_DBL y1 = 2;
        REF_DBL h = 1.0 / 100.0;
        REF_DBL slen = dist[node];
        REF_DBL turb1 = field[5 + ldim * node];
        if (x0 <= ref_node_xyz(ref_node, 0, node) &&
            ref_node_xyz(ref_node, 0, node) <= x1 &&
            y0 <= ref_node_xyz(ref_node, 1, node) &&
            ref_node_xyz(ref_node, 1, node) <= y1 && (4 <= turb1 || h > slen)) {
          m[0] = 1.0 / (h * h);
          m[1] = 0.0;
          m[2] = 0.0;
          m[3] = 1.0 / (h * h);
          m[4] = 0.0;
          m[5] = 1.0;
          for (i = 0; i < 6; i++) m0[i] = metric[i + 6 * node];
          RSS(ref_matrix_intersect(m0, m, &(metric[6 * node])), "intersect");
        }
      }
    } else {
      each_ref_node_valid_node(ref_node, node) {
        REF_DBL x0 = 530;
        REF_DBL x1 = 670;
        REF_DBL y0 = -572.8;
        REF_DBL y1 = -507.4;
        REF_DBL z0 = 100.0;
        REF_DBL z1 = 151.2;
        REF_DBL h = 0.25;
        REF_DBL slen = dist[node];
        REF_DBL turb1 = field[5 + ldim * node];
        if (x0 <= ref_node_xyz(ref_node, 0, node) &&
            ref_node_xyz(ref_node, 0, node) <= x1 &&
            y0 <= ref_node_xyz(ref_node, 1, node) &&
            ref_node_xyz(ref_node, 1, node) <= y1 &&
            z0 <= ref_node_xyz(ref_node, 2, node) &&
            ref_node_xyz(ref_node, 2, node) <= z1 && (4 <= turb1 || h > slen)) {
          m[0] = 1.0 / (h * h);
          m[1] = 0.0;
          m[2] = 0.0;
          m[3] = 1.0 / (h * h);
          m[4] = 0.0;
          m[5] = 1.0 / (h * h);
          for (i = 0; i < 6; i++) m0[i] = metric[i + 6 * node];
          RSS(ref_matrix_intersect(m0, m, &(metric[6 * node])), "intersect");
        }
      }
    }
    RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");
    ref_mpi_stopwatch_stop(ref_mpi, "intersect");

    for (gradation = 0; gradation < 20; gradation++) {
      RSS(ref_metric_mixed_space_gradation(metric, ref_grid, -1.0, -1.0),
          "grad");
      ref_mpi_stopwatch_stop(ref_mpi, "gradation");
    }

    RSS(ref_metric_to_node(metric, ref_node), "set node");
    ref_free(metric);
    ref_free(field);
    ref_free(dist);

    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("writing metric %s\n", argv[5]);
    RSS(ref_gather_metric(ref_grid, argv[5]), "export scaled metric");
    ref_mpi_stopwatch_stop(ref_mpi, "dump metric");

    RSS(ref_grid_free(ref_grid), "free");
    ref_mpi_stopwatch_stop(ref_mpi, "done.");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (argc == 3) {
    REF_GRID ref_grid;

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");

    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[2]), "get metric");

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "qual");
    RSS(ref_histogram_ratio(ref_grid), "rat");

    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_metric_test_s00"), "al");

    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* imply metric right tet */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_from(metric, ref_grid), "imply");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(1.0, metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.5, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.5, metric[2 + 6 * node], tol, "m[2]");
      RWDS(1.0, metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.5, metric[4 + 6 * node], tol, "m[4]");
      RWDS(1.0, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* imply metric right prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_from(metric, ref_grid), "imply");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(1.259596e+00, metric[0 + 6 * node], tol, "m[0]");
      RWDS(-6.394720e-01, metric[1 + 6 * node], tol, "m[1]");
      RWDS(6.394720e-01, metric[2 + 6 * node], tol, "m[2]");
      RWDS(9.546890e-01, metric[3 + 6 * node], tol, "m[3]");
      RWDS(-3.247285e-01, metric[4 + 6 * node], tol, "m[4]");
      RWDS(9.546890e-01, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* imply metric non tet prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_non_tet(metric, ref_grid), "imply");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(1.259596e+00, metric[0 + 6 * node], tol, "m[0]");
      RWDS(-6.394720e-01, metric[1 + 6 * node], tol, "m[1]");
      RWDS(6.394720e-01, metric[2 + 6 * node], tol, "m[2]");
      RWDS(9.546890e-01, metric[3 + 6 * node], tol, "m[3]");
      RWDS(-3.247285e-01, metric[4 + 6 * node], tol, "m[4]");
      RWDS(9.546890e-01, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* imply metric tet prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_pri_tet_cap_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_from(metric, ref_grid), "imply");

    node = 6;
    RWDS(1.00, metric[0 + 6 * node], tol, "m[0]");
    RWDS(0.50, metric[1 + 6 * node], tol, "m[1]");
    RWDS(0.05, metric[2 + 6 * node], tol, "m[2]");
    RWDS(1.00, metric[3 + 6 * node], tol, "m[3]");
    RWDS(0.05, metric[4 + 6 * node], tol, "m[4]");
    RWDS(0.67, metric[5 + 6 * node], tol, "m[5]");

    RSS(ref_metric_imply_non_tet(metric, ref_grid), "imply");

    node = 6;
    RWDS(1.00, metric[0 + 6 * node], tol, "m[0]");
    RWDS(0.50, metric[1 + 6 * node], tol, "m[1]");
    RWDS(0.05, metric[2 + 6 * node], tol, "m[2]");
    RWDS(1.00, metric[3 + 6 * node], tol, "m[3]");
    RWDS(0.05, metric[4 + 6 * node], tol, "m[4]");
    RWDS(0.67, metric[5 + 6 * node], tol, "m[5]");

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* imply metric pyr */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_pyr_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_from(metric, ref_grid), "imply");

    node = 2;
    RWDS(9.690913e-01, metric[0 + 6 * node], tol, "m[0]");
    RWDS(2.525674e-01, metric[1 + 6 * node], tol, "m[1]");
    RWDS(-4.689688e-01, metric[2 + 6 * node], tol, "m[2]");
    RWDS(9.503095e-01, metric[3 + 6 * node], tol, "m[3]");
    RWDS(2.525674e-01, metric[4 + 6 * node], tol, "m[4]");
    RWDS(9.690913e-01, metric[5 + 6 * node], tol, "m[5]");

    RSS(ref_metric_imply_non_tet(metric, ref_grid), "imply");

    node = 2;
    RWDS(9.690913e-01, metric[0 + 6 * node], tol, "m[0]");
    RWDS(2.525674e-01, metric[1 + 6 * node], tol, "m[1]");
    RWDS(-4.689688e-01, metric[2 + 6 * node], tol, "m[2]");
    RWDS(9.503095e-01, metric[3 + 6 * node], tol, "m[3]");
    RWDS(2.525674e-01, metric[4 + 6 * node], tol, "m[4]");
    RWDS(9.690913e-01, metric[5 + 6 * node], tol, "m[5]");

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* imply metric hex */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_hex_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    RSS(ref_metric_imply_from(metric, ref_grid), "imply");

    node = 0;
    RWDS(7.854306e-01, metric[0 + 6 * node], tol, "m[0]");
    RWDS(1.331808e-01, metric[1 + 6 * node], tol, "m[1]");
    RWDS(1.331808e-01, metric[2 + 6 * node], tol, "m[2]");
    RWDS(9.960985e-01, metric[3 + 6 * node], tol, "m[3]");
    RWDS(-5.352162e-01, metric[4 + 6 * node], tol, "m[4]");
    RWDS(9.960985e-01, metric[5 + 6 * node], tol, "m[5]");

    RSS(ref_metric_imply_non_tet(metric, ref_grid), "imply");

    node = 0;
    RWDS(7.854306e-01, metric[0 + 6 * node], tol, "m[0]");
    RWDS(1.331808e-01, metric[1 + 6 * node], tol, "m[1]");
    RWDS(1.331808e-01, metric[2 + 6 * node], tol, "m[2]");
    RWDS(9.960985e-01, metric[3 + 6 * node], tol, "m[3]");
    RWDS(-5.352162e-01, metric[4 + 6 * node], tol, "m[4]");
    RWDS(9.960985e-01, metric[5 + 6 * node], tol, "m[5]");

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* smr octave */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric_file;
    REF_DBL *metric_imply;
    REF_DBL *metric;
    REF_INT node;

    /*
  clear all
  format long
  a1 = [ 1 0 0 ;
         0 1 0 ;
         0 0 1 ]
  a2 = [ 30 -25   -2.5;
        -25  25   1;
        -2.5   1  1.5];
  a3 = inv(a1)*a2
  [vector3, value3] = eig(a3)

  val1 = vector3'*a1*vector3
  val2 = vector3'*a2*vector3
  for ii=1:3
   h1(ii,ii) = sqrt(1/val1(ii,ii));
   h2(ii,ii) = sqrt(1/val2(ii,ii));
   h(ii,ii) = max(0.25*h1(ii,ii),min(4.0*h1(ii,ii),h2(ii,ii)));
   val(ii,ii)=1.0 / (h(ii,ii)* h(ii,ii));
  end
  h1
  h2
  h
  val
  vector = inv(vector3)

  smr = vector'*val*vector
  [sv,se]= eig(smr)
  se.^-0.5
    */

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_node_add(ref_grid_node(ref_grid), 0, &node), "add");

    ref_malloc(metric_file, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(metric_imply, 6 * ref_node_max(ref_grid_node(ref_grid)),
               REF_DBL);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    metric_imply[0 + 6 * node] = 1.0;
    metric_imply[1 + 6 * node] = 0.0;
    metric_imply[2 + 6 * node] = 0.0;
    metric_imply[3 + 6 * node] = 1.0;
    metric_imply[4 + 6 * node] = 0.0;
    metric_imply[5 + 6 * node] = 1.0;

    metric_file[0 + 6 * node] = 30.0;
    metric_file[1 + 6 * node] = -25.0;
    metric_file[2 + 6 * node] = -2.5;
    metric_file[3 + 6 * node] = 25.0;
    metric_file[4 + 6 * node] = 1.0;
    metric_file[5 + 6 * node] = 1.5;

    RSS(ref_metric_smr(metric_imply, metric_file, metric, ref_grid), "smr");

    RWDS(9.812655244359012, metric[0 + 6 * node], tol, "m[0]");
    RWDS(-6.764301991389758, metric[1 + 6 * node], tol, "m[0]");
    RWDS(-1.159409438169853, metric[2 + 6 * node], tol, "m[0]");
    RWDS(8.527269886828027, metric[3 + 6 * node], tol, "m[0]");
    RWDS(-0.210986632201670, metric[4 + 6 * node], tol, "m[0]");
    RWDS(1.410974767795262, metric[5 + 6 * node], tol, "m[0]");

    ref_free(metric);
    ref_free(metric_imply);
    ref_free(metric_file);
    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_GRID ref_grid, parent_grid;
    REF_INTERP ref_interp;
    REF_INT node, im;
    REF_DBL tol = -1.0;
    REF_DBL parent_m[6];
    REF_DBL child_m[6];

    RSS(ref_fixture_tet_brick_grid(&parent_grid, ref_mpi), "brick");
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    RSS(ref_metric_olympic_node(ref_grid_node(parent_grid), 0.001), "oly");

    RSS(ref_interp_create(&ref_interp, ref_grid, parent_grid), "map");
    RSS(ref_interp_locate(ref_interp), "map");
    RSS(ref_metric_interpolate(ref_interp), "interp");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RSS(ref_node_metric_get(ref_grid_node(parent_grid), node, parent_m),
          "get parent m");
      RSS(ref_node_metric_get(ref_grid_node(ref_grid), node, child_m),
          "get child m");
      for (im = 0; im < 6; im++) {
        RWDS(parent_m[im], child_m[im], tol, "interpolant");
      }
    }

    RSS(ref_interp_free(ref_interp), "free");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(parent_grid), "free");
  }

  {
    REF_GRID ref_grid, truth;
    REF_INT node, im;
    REF_DBL tol = -1.0;
    char meshb[] = "ref_metric_test.meshb";
    char solb[] = "ref_metric_test-metric.solb";
    REF_DBL truth_m[6];
    REF_DBL m[6];

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, meshb), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }

    RSS(ref_part_by_extension(&truth, ref_mpi, meshb), "import");
    RSS(ref_metric_ugawg_node(ref_grid_node(truth), 1), "m");
    RSS(ref_gather_metric(truth, solb), "export");

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, meshb), "import");
    RSS(ref_part_metric(ref_grid_node(ref_grid), solb), "export");

    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove(meshb), "test meshb clean up");
      REIS(0, remove(solb), "test solb clean up");
    }

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RSS(ref_node_metric_get(ref_grid_node(truth), node, truth_m),
          "get truth m");
      RSS(ref_node_metric_get(ref_grid_node(ref_grid), node, m), "get m");
      for (im = 0; im < 6; im++) {
        RWDS(truth_m[im], m[im], tol, "interpolant");
      }
    }

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(truth), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* gradation */
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 1.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 1.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 1.0;
    }
    node = 0;
    metric[0 + 6 * node] = 1.0;
    metric[1 + 6 * node] = 0.0;
    metric[2 + 6 * node] = 0.0;
    metric[3 + 6 * node] = 1.0;
    metric[4 + 6 * node] = 0.0;
    metric[5 + 6 * node] = 4.0;

    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    RSS(ref_metric_metric_space_gradation(metric, ref_grid, 1.1), "grad");

    node = 0;
    RWDS(1.0, metric[0 + 6 * node], tol, "m[0]");
    RWDS(1.0, metric[3 + 6 * node], tol, "m[3]");
    RWDS(4.0, metric[5 + 6 * node], tol, "m[5]");

    node = 3;
    RWDS(1.0, metric[0 + 6 * node], tol, "m[0]");
    RWDS(1.0, metric[3 + 6 * node], tol, "m[3]");
    RWDS(2.821716527185583, metric[5 + 6 * node], tol, "m[5]");

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* limit hmin */
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;
    REF_DBL hmax, hmin;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 1.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 1.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 1.0;
    }

    hmin = 2.0;
    hmax = -1.0;
    RSS(ref_metric_limit_h(metric, ref_grid, hmin, hmax), "h limit");

    /* m = hmin**-2  */
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.25, metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.25, metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.25, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* limit hmax */
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;
    REF_DBL hmax, hmin;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 1.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 1.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 1.0;
    }

    hmin = -1.0;
    hmax = 0.5;
    RSS(ref_metric_limit_h(metric, ref_grid, hmin, hmax), "h limit");

    /* m = hmax**-2  */
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(4.0, metric[0 + 6 * node], tol, "m[0]");
      RWDS(4.0, metric[3 + 6 * node], tol, "m[3]");
      RWDS(4.0, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* limit no-op */
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;
    REF_DBL hmax, hmin;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 1.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 1.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 1.0;
    }

    hmin = 0.5;
    hmax = 2.0;
    RSS(ref_metric_limit_h(metric, ref_grid, hmin, hmax), "h limit");

    /* unlimited  */
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(1.0, metric[0 + 6 * node], tol, "m[0]");
      RWDS(1.0, metric[3 + 6 * node], tol, "m[3]");
      RWDS(1.0, metric[5 + 6 * node], tol, "m[5]");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* lp for small variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *metric;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = 0.5 + 0.01 * pow(ref_node_xyz(ref_node, 0, node), 2) +
                     0.02 * pow(ref_node_xyz(ref_node, 1, node), 2) +
                     0.03 * pow(ref_node_xyz(ref_node, 2, node), 2);
    }
    RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, REF_RECON_L2PROJECTION, 2,
                      1.5, 1000.0),
        "lp norm");
    ref_free(metric);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* lp for no variation */
    REF_GRID ref_grid;
    REF_INT node;
    REF_DBL *scalar, *metric;
    REF_DBL current_complexity;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = 0.5;
    }
    RSS(ref_metric_lp(metric, ref_grid, scalar, NULL, REF_RECON_L2PROJECTION, 2,
                      1.5, 1000.0),
        "const metric");
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    RWDS(1000.0, current_complexity, -1.0, "complexity");
    ref_free(metric);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* lp for no variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node, ldim = 10;
    REF_DBL *prim_dual, *metric;
    REF_DBL t, ei0, et0;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    REF_DBL mach = 0.5;
    REF_DBL re = 1.0e6;
    REF_DBL reference_temp = 273.11;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(prim_dual, ldim * ref_node_max(ref_grid_node(ref_grid)),
               REF_DBL);
    each_ref_node_valid_node(ref_node, node) {
      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[0 + ldim * node] = 1.0 + 0.01 * cos(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[1 + ldim * node] = 0.5 + 0.1 * sin(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 1, node);
      prim_dual[2 + ldim * node] = 0.0 + 0.1 * cos(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 2, node);
      prim_dual[3 + ldim * node] = 0.1 + 0.1 * sin(t);
      ei0 = (1.0 / 1.4) / ((1.4 - 1.0) * 1.0);
      et0 = 1.0 * (ei0 + 0.5 * (0.5 * 0.5 + 0.1 * 0.1));
      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[4 + ldim * node] = et0 + 0.01 * sin(t);

      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[5 + ldim * node] = 1.0 * cos(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[6 + ldim * node] = 2.0 * sin(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 1, node);
      prim_dual[7 + ldim * node] = 2.0 * sin(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 2, node);
      prim_dual[8 + ldim * node] = 2.0 * sin(t);
      t = ref_math_pi * ref_node_xyz(ref_node, 0, node);
      prim_dual[9 + ldim * node] = 5.0 * cos(t);
    }

    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);

    RSS(ref_metric_belme_gfe(metric, ref_grid, ldim, prim_dual, reconstruction),
        "gfe");
    RSS(ref_metric_belme_gu(metric, ref_grid, ldim, prim_dual, mach, re,
                            reference_temp, reconstruction),
        "gu");

    ref_free(metric);
    ref_free(prim_dual);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* parse interior box floor spacing */
    char *args[] = {
        "--uniform", "box", "floor", "2", "-1", "0", "0", "0", "1", "1", "1",
    };
    int narg = 11;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.25, metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(0.25, metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(0.25, metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* parse outside box spacing */
    char *args[] = {
        "--uniform", "box", "floor", "3", "-0.5", "-1",
        "-1",        "-1",  "0",     "0", "0",
    };
    int narg = 11;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL r, h;
      r = sqrt(pow(ref_node_xyz(ref_grid_node(ref_grid), 0, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 1, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 2, node), 2));
      h = 3.0 * pow(2.0, -r / -0.5);
      RWDS(1.0 / (h * h), metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(1.0 / (h * h), metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(1.0 / (h * h), metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* parse interior box ceil spacing */
    char *args[] = {
        "--uniform", "box", "ceil", "0.25", "-1", "0", "0", "0", "1", "1", "1",
    };
    int narg = 11;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(16.0, metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(16.0, metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(16.0, metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* parse outside box spacing */
    char *args[] = {
        "--uniform", "box", "ceil", "3", "-0.5", "-1",
        "-1",        "-1",  "0",    "0", "0",
    };
    int narg = 11;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL r, h;
      r = sqrt(pow(ref_node_xyz(ref_grid_node(ref_grid), 0, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 1, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 2, node), 2));
      h = MIN(3.0 * pow(2.0, -r / -0.5), 0.5);
      RWDS(1.0 / (h * h), metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(1.0 / (h * h), metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(1.0 / (h * h), metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* parse interior box floor and box ceil spacing, last wins */
    char *args[] = {
        "--uniform", "box", "floor", "4", "-1", "0", "0", "0", "1", "1", "1",
        "--uniform", "box", "ceil",  "2", "-1", "0", "0", "0", "1", "1", "1",
    };
    int narg = 22;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.25, metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(0.25, metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(0.25, metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* distance to truncated cone (cylinder) */
    REF_DBL cone_geom[] = {0, 0, 0, 1, 0, 0, 1, 1};
    REF_DBL dist;
    REF_DBL xyz[3];
    REF_DBL tol = -1.0;
    /* inside */
    xyz[0] = 0;
    xyz[1] = 0;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.0, dist, tol, "inside");
    /* negative axis */
    xyz[0] = -1;
    xyz[1] = 0.2;
    xyz[2] = 0.3;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(1.0, dist, tol, "neg axis");
    /* positive axis */
    xyz[0] = 3;
    xyz[1] = 0.4;
    xyz[2] = 0.5;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(2.0, dist, tol, "pos axis");
    /* ra end circle position */
    xyz[0] = -1;
    xyz[1] = 2;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(2), dist, tol, "a ra circle");
    /* inside middle radial position */
    xyz[0] = 0.5;
    xyz[1] = 0.5;
    xyz[2] = 0.5;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.0, dist, tol, "inside");
    /* above ra radial position */
    xyz[0] = 0;
    xyz[1] = 2;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(1, dist, tol, "a ra circle");
    /* rb end circle position */
    xyz[0] = 2;
    xyz[1] = 2;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(2), dist, tol, "a rb circle");
    /* outside middle radial position y */
    xyz[0] = 0.5;
    xyz[1] = 2;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(1, dist, tol, "outside middle radial y");
    /* outside middle radial position z */
    xyz[0] = 0.5;
    xyz[1] = 0;
    xyz[2] = 2;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(1, dist, tol, "outside middle radial z");
    /* outside middle radial position yz */
    xyz[0] = 0.5;
    xyz[1] = 2;
    xyz[2] = 2;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(2 * 2 + 2 * 2) - 1, dist, tol, "outside middle radial yz");
  }

  { /* distance to truncated cone (sphere) */
    REF_DBL cone_geom[] = {0, 0, 0, 0, 0, 0, 1, 1.5};
    REF_DBL dist;
    REF_DBL xyz[3];
    REF_DBL tol = -1.0;
    /* inside */
    xyz[0] = 0;
    xyz[1] = 0;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.0, dist, tol, "inside");
    /* outside */
    xyz[0] = 2;
    xyz[1] = 0;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.5, dist, tol, "outside");
  }

  { /* distance to truncated cone (core) */
    REF_DBL cone_geom[] = {0, 0, 0, 1, 0, 0, 0, 1};
    REF_DBL dist;
    REF_DBL xyz[3];
    REF_DBL tol = -1.0;
    /* inside */
    xyz[0] = 0.5;
    xyz[1] = 0.2;
    xyz[2] = 0.3;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.0, dist, tol, "inside");
    /* tip */
    xyz[0] = -2;
    xyz[1] = -1;
    xyz[2] = -1;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(6), dist, tol, "tip");
    /* side */
    xyz[0] = 0;
    xyz[1] = -1;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(0.5 * sqrt(2), dist, tol, "side");
    /* side 11 */
    xyz[0] = -1;
    xyz[1] = -1;
    xyz[2] = 0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(2), dist, tol, "side 11");
  }

  { /* distance to truncated cone (core) */
    REF_DBL cone_geom[] = {-1, -1, -1, -2, -2, -2, 0, 1};
    REF_DBL dist;
    REF_DBL xyz[3];
    REF_DBL tol = -1.0;
    /* inside */
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    RSS(ref_metric_truncated_cone_dist(cone_geom, xyz, &dist), "d");
    RWDS(sqrt(3), dist, tol, "inside");
  }

  { /* parse outside box spacing */
    char *args[] = {
        "--uniform", "cyl", "ceil", "3",  "-0.5", "-1", "-1",
        "-1",        "-2",  "-2",   "-2", "0",    "1",
    };
    int narg = 13;
    REF_DBL tol = -1.0;
    REF_DBL *metric;
    REF_GRID ref_grid;
    REF_INT node;
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_malloc_init(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 4.0;
      metric[1 + 6 * node] = 0.0;
      metric[2 + 6 * node] = 0.0;
      metric[3 + 6 * node] = 4.0;
      metric[4 + 6 * node] = 0.0;
      metric[5 + 6 * node] = 4.0;
    }
    RSS(ref_metric_parse(metric, ref_grid, narg, args), "parse");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL r, h;
      r = sqrt(pow(ref_node_xyz(ref_grid_node(ref_grid), 0, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 1, node), 2) +
               pow(ref_node_xyz(ref_grid_node(ref_grid), 2, node), 2));
      h = MIN(3.0 * pow(2.0, -r / -0.5), 0.5);
      RWDS(1.0 / (h * h), metric[0 + 6 * node], tol, "m[0]");
      RWDS(0.00, metric[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, metric[2 + 6 * node], tol, "m[2]");
      RWDS(1.0 / (h * h), metric[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, metric[4 + 6 * node], tol, "m[4]");
      RWDS(1.0 / (h * h), metric[5 + 6 * node], tol, "m[5]");
    }
    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
