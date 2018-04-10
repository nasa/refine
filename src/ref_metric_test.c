
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

#include "ref_metric.h"

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_node.h"
#include "ref_sort.h"

#include "ref_fixture.h"
#include "ref_malloc.h"

#include "ref_mpi.h"

#include "ref_dict.h"
#include "ref_export.h"
#include "ref_import.h"
#include "ref_migrate.h"
#include "ref_part.h"

#include "ref_adapt.h"
#include "ref_edge.h"

#include "ref_args.h"
#include "ref_collapse.h"
#include "ref_face.h"
#include "ref_gather.h"
#include "ref_histogram.h"
#include "ref_math.h"
#include "ref_smooth.h"
#include "ref_split.h"
#include "ref_twod.h"
#include "ref_validation.h"

#include "ref_geom.h"

#include "ref_clump.h"

/*
./test.sh ref_metric && ./ref_metric_test \
  ../acceptance/2d/linear/two/accept-2d-two-09.b8.ugrid \
  --parent ../acceptance/2d/linear/two/accept-2d-two-08.b8.ugrid \
  ../acceptance/2d/linear/two/accept-2d-two-08.metric
*/

int main(int argc, char *argv[]) {
  REF_INT curve_limit_pos = REF_EMPTY;
  REF_INT curvature_pos = REF_EMPTY;
  REF_INT parent_pos = REF_EMPTY;
  REF_INT xyzdirlen_pos = REF_EMPTY;
  REF_INT lp_pos = REF_EMPTY;

  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--curve-limit", &curve_limit_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--curvature", &curvature_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--parent", &parent_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--xyzdirlen", &xyzdirlen_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--lp", &lp_pos), REF_NOT_FOUND, "arg search");

  if (curve_limit_pos != REF_EMPTY) {
    REF_GRID ref_grid;

    REIS(1, curve_limit_pos,
         "required args: --curve-limit grid.ext input.metric geom.egads "
         "[assoc.gas]");
    RAS(argc == 5 || argc == 6,
        "required args: --curve-limit grid.ext input.metric geom.egads "
        "[assoc.gas]");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 1");
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load parent metric in position 2");
    RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[4]),
        "unable to load egads in position 3");
    if (argc == 6)
      RSS(ref_geom_load(ref_grid, argv[5]),
          "unable to load geom assoc in position 4");

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

  if (lp_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, *metric;
    REF_INT p;
    REF_DBL gradation, complexity, current_complexity;
    REIS(1, lp_pos,
         "required args: --lp grid.meshb scalar-mach.solb p gradation "
         "complexity output-metric.solb");
    REIS(8, argc,
         "required args: --lp grid.meshb scalar-mach.solb p gradation "
         "complexity output-metric.solb");

    p = atoi(argv[4]);
    gradation = atof(argv[5]);
    complexity = atof(argv[6]);
    printf("Lp=%d\n", p);
    printf("gradation %f\n", gradation);
    printf("complexity %f\n", complexity);

    printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    printf("reading scalar %s\n", argv[3]);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), scalar, argv[3]),
        "unable to load scalar in position 3");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_lp(metric, ref_grid, scalar, p, gradation, complexity),
        "lp norm");
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("actual complexity %e\n", current_complexity);
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);
    ref_free(scalar);

    printf("writing metric %s\n", argv[7]);
    RSS(ref_gather_metric(ref_grid, argv[7]), "export curve limit metric");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (curvature_pos != REF_EMPTY) {
    REF_GRID ref_grid;

    REIS(1, curvature_pos,
         "required args: --curvature grid.ext geom.egads assoc.gas");
    REIS(5, argc, "required args: --curvature grid.ext geom.egads assoc.gas");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 1");
    RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[3]),
        "unable to load egads in position 2");
    RSS(ref_geom_load(ref_grid, argv[4]),
        "unable to load geom assoc in position 3");

    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");

    RSS(ref_gather_metric(ref_grid, "ref_metric_test_curve.metric"),
        "export curve metric");

    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_metric_test_curve"), "al");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (parent_pos != REF_EMPTY) {
    REF_GRID ref_grid, parent_grid;

    REIS(2, parent_pos,
         "required args: grid.ext --parent pgrid.ext pgrid.metric");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]),
        "unable to load target grid in position 1");
    RSS(ref_import_by_extension(&parent_grid, ref_mpi, argv[3]),
        "unable to load parent grid in position 3");
    RSS(ref_part_metric(ref_grid_node(parent_grid), argv[4]),
        "unable to load parent grid in position 4");

    RSS(ref_metric_interpolate(ref_grid, parent_grid), "interp");

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "qual");
    RSS(ref_histogram_ratio(ref_grid), "rat");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(parent_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (xyzdirlen_pos != REF_EMPTY) {
    REF_GRID ref_grid;

    REIS(1, xyzdirlen_pos,
         "required args: --xyzdirlen grid.ext input.metric output");
    RAS(argc == 5, "required args: --xyzdirlen grid.ext input.metric output");

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 1");
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[3]),
        "unable to load parent grid in position 2");
    RSS(ref_export_metric_xyzdirlen(ref_grid, argv[4]),
        "export metric in xyzdirlen");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (argc == 2) {
    REF_GRID ref_grid;
    REF_DBL *metric;

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_metric_imply_from(metric, ref_grid), "imply");
    RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "set node");
    ref_free(metric);

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "qual");
    RSS(ref_histogram_ratio(ref_grid), "rat");

    RSS(ref_export_tec_metric_ellipse(ref_grid, "ref_metric_test_implied_met"),
        "al");

    RSS(ref_grid_free(ref_grid), "free");

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

  if (argc == 4) {
    REF_GRID ref_grid;

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "read grid");

    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[2]), "get metric");

    RSS(ref_metric_sanitize(ref_grid), "sant metric");

    RSS(ref_gather_metric(ref_grid, argv[3]), "in");

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

  if (!ref_mpi_para(ref_mpi)) {
    REF_GRID ref_grid, parent_grid;
    REF_INT node, im;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_twod_brick_grid(&parent_grid, ref_mpi), "brick");
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "brick");

    RSS(ref_metric_olympic_node(ref_grid_node(parent_grid), 0.001), "oly");
    RSS(ref_metric_twod_node(ref_grid_node(parent_grid)), "2d");

    RSS(ref_metric_interpolate(ref_grid, parent_grid), "interp");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) for (im = 0; im < 6;
                                                                 im++)
        RWDS(ref_node_metric(ref_grid_node(parent_grid), im, node),
             ref_node_metric(ref_grid_node(ref_grid), im, node), tol,
             "interpolant");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(parent_grid), "free");
  }

  {
    REF_GRID ref_grid, parent_grid;
    REF_INT node, im;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_brick_grid(&parent_grid, ref_mpi), "brick");
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    RSS(ref_metric_olympic_node(ref_grid_node(parent_grid), 0.001), "oly");

    RSS(ref_metric_interpolate(ref_grid, parent_grid), "interp");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) for (im = 0; im < 6;
                                                                 im++)
        RWDS(ref_node_metric(ref_grid_node(parent_grid), im, node),
             ref_node_metric(ref_grid_node(ref_grid), im, node), tol,
             "interpolant");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(parent_grid), "free");
  }

  {
    REF_GRID ref_grid, truth;
    REF_INT node, im;
    REF_DBL tol = -1.0;
    char meshb[] = "ref_metric_test.meshb";
    char solb[] = "ref_metric_test-metric.solb";

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

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) for (im = 0; im < 6;
                                                                 im++)
        RWDS(ref_node_metric(ref_grid_node(truth), im, node),
             ref_node_metric(ref_grid_node(ref_grid), im, node), tol,
             "metric wrong");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_grid_free(truth), "free");
  }

  { /* l2-projection grad */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *scalar, *grad;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(grad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = 1.3 * ref_node_xyz(ref_grid_node(ref_grid), 0, node) +
                     3.5 * ref_node_xyz(ref_grid_node(ref_grid), 1, node) +
                     7.2 * ref_node_xyz(ref_grid_node(ref_grid), 2, node) +
                     15.0;
    }

    RSS(ref_metric_l2_projection_grad(ref_grid, scalar, grad), "l2 grad");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(1.3, grad[0 + 3 * node], tol, "gradx");
      RWDS(3.5, grad[1 + 3 * node], tol, "grady");
      RWDS(7.2, grad[2 + 3 * node], tol, "gradz");
    }

    ref_free(grad);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* l2-projection hessian zero, constant gradient */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *scalar, *hessian;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(hessian, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = 1.3 * ref_node_xyz(ref_grid_node(ref_grid), 0, node) +
                     3.5 * ref_node_xyz(ref_grid_node(ref_grid), 1, node) +
                     7.2 * ref_node_xyz(ref_grid_node(ref_grid), 2, node) +
                     15.0;
    }

    RSS(ref_metric_l2_projection_hessian(ref_grid, scalar, hessian), "l2 hess");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.0, hessian[0 + 6 * node], tol, "m11");
      RWDS(0.0, hessian[1 + 6 * node], tol, "m12");
      RWDS(0.0, hessian[2 + 6 * node], tol, "m13");
      RWDS(0.0, hessian[3 + 6 * node], tol, "m22");
      RWDS(0.0, hessian[4 + 6 * node], tol, "m23");
      RWDS(0.0, hessian[5 + 6 * node], tol, "m33");
    }

    ref_free(hessian);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* boundary averaging constant metric */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 100.0;
      metric[1 + 6 * node] = 7.0;
      metric[2 + 6 * node] = 22.0;
      metric[3 + 6 * node] = 200.0;
      metric[4 + 6 * node] = 15.0;
      metric[5 + 6 * node] = 300.0;
    }

    RSS(ref_metric_extrapolate_boundary(metric, ref_grid), "bound extrap");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(100.0, metric[0 + 6 * node], tol, "m11");
      RWDS(7.0, metric[1 + 6 * node], tol, "m12");
      RWDS(22.0, metric[2 + 6 * node], tol, "m13");
      RWDS(200.0, metric[3 + 6 * node], tol, "m22");
      RWDS(15.0, metric[4 + 6 * node], tol, "m23");
      RWDS(300.0, metric[5 + 6 * node], tol, "m33");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* multipass boundary averaging constant metric */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      metric[0 + 6 * node] = 100.0;
      metric[1 + 6 * node] = 7.0;
      metric[2 + 6 * node] = 22.0;
      metric[3 + 6 * node] = 200.0;
      metric[4 + 6 * node] = 15.0;
      metric[5 + 6 * node] = 300.0;
    }

    RSS(ref_metric_extrapolate_boundary_multipass(metric, ref_grid),
        "bound extrap");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(100.0, metric[0 + 6 * node], tol, "m11");
      RWDS(7.0, metric[1 + 6 * node], tol, "m12");
      RWDS(22.0, metric[2 + 6 * node], tol, "m13");
      RWDS(200.0, metric[3 + 6 * node], tol, "m22");
      RWDS(15.0, metric[4 + 6 * node], tol, "m23");
      RWDS(300.0, metric[5 + 6 * node], tol, "m33");
    }

    ref_free(metric);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) {
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
    RSS(ref_metric_gradation(metric, ref_grid, 1.1), "grad");

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

  if ( ref_mpi_once(ref_mpi) ) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *metric;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      scalar[node] = 0.5 + 0.01*pow(ref_node_xyz(ref_node,0,node),2) + 
        0.02*pow(ref_node_xyz(ref_node,1,node),2) + 
        0.03*pow(ref_node_xyz(ref_node,2,node),2);
    }
    RSS(ref_metric_lp(metric, ref_grid, scalar, 2, 1.5, 1000.0),
        "lp norm");
    ref_free(metric);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
