
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_interp.h"

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cavity.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_twod.h"

#include "ref_args.h"
#include "ref_malloc.h"

static REF_STATUS ref_interp_setup(REF_INTERP *ref_interp_ptr,
                                   REF_MPI ref_mpi) {
  REF_GRID from, to;
  RSS(ref_grid_create(&from, ref_mpi), "create");
  RSS(ref_grid_create(&to, ref_mpi), "create");
  RSS(ref_interp_create(ref_interp_ptr, from, to), "make interp");
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_teardown(REF_INTERP ref_interp) {
  RSS(ref_grid_free(ref_interp_from_grid(ref_interp)), "free");
  RSS(ref_grid_free(ref_interp_to_grid(ref_interp)), "free");
  RSS(ref_interp_free(ref_interp), "interp free");
  return REF_SUCCESS;
}

static REF_STATUS ref_interp_shift_cube_interior(REF_NODE ref_node) {
  REF_INT node;
  each_ref_node_valid_node(ref_node, node) {
    if ((0.01 < ref_node_xyz(ref_node, 0, node) &&
         0.99 > ref_node_xyz(ref_node, 0, node)) ||
        (0.01 < ref_node_xyz(ref_node, 1, node) &&
         0.99 > ref_node_xyz(ref_node, 1, node)) ||
        (0.01 < ref_node_xyz(ref_node, 2, node) &&
         0.99 > ref_node_xyz(ref_node, 2, node))) {
      ref_node_xyz(ref_node, 0, node) += 1.0e-2;
      ref_node_xyz(ref_node, 1, node) += 2.0e-2;
      ref_node_xyz(ref_node, 2, node) += 4.0e-2;
    }
  }
  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_INT pair_pos = REF_EMPTY;
  REF_INT rate_pos = REF_EMPTY;
  REF_INT error_pos = REF_EMPTY;
  REF_INT subset_pos = REF_EMPTY;
  REF_INT field_pos = REF_EMPTY;
  REF_INT mach_pos = REF_EMPTY;
  REF_INT cust_pos = REF_EMPTY;
  REF_INT entropy_pos = REF_EMPTY;
  REF_INT entropyadj_pos = REF_EMPTY;
  REF_INT heat_pos = REF_EMPTY;
  REF_INT gamma_pos = REF_EMPTY;
  REF_INT plt_pos = REF_EMPTY;

  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");
  ref_mpi_stopwatch_start(ref_mpi);

  RXS(ref_args_find(argc, argv, "--pair", &pair_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--rate", &rate_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--error", &error_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--subset", &subset_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--field", &field_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--mach", &mach_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--cust", &cust_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--entropy", &entropy_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--entropyadj", &entropyadj_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--heat", &heat_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--gamma", &gamma_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--plt", &plt_pos), REF_NOT_FOUND,
      "arg search");

  if (REF_EMPTY != plt_pos) {
    REF_GRID ref_grid;
    REF_DBL *scalar;
    REF_INT ldim;
    REIS(1, plt_pos,
         "required args: --plt grid.ext usm3d-solution.plt "
         "stitched-soln.solb\n");
    if (5 > argc) {
      printf(
          "required args: --plt grid.ext usm3d-solution.plt "
          "stitched-soln.solb\n");
      return REF_FAILURE;
    }
    RSS(ref_mpi_stopwatch_start(ref_mpi), "sw start");

    if (ref_mpi_once(ref_mpi)) printf("read grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]), "import");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "read grid"), "sw start");

    if (ref_mpi_once(ref_mpi)) printf("read/stitch plt %s\n", argv[3]);
    RSS(ref_iterp_plt(ref_grid, argv[3], &ldim, &scalar), "plt zone");
    ref_mpi_stopwatch_stop(ref_mpi, "read/stitch zones");

    if (ref_mpi_once(ref_mpi)) printf("write/gather stitched %s\n", argv[4]);
    RSS(ref_gather_scalar(ref_grid, ldim, scalar, argv[4]),
        "write/gather stitched");
    ref_mpi_stopwatch_stop(ref_mpi, "write stitched");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != pair_pos) {
    REF_GRID from, to;
    REF_INTERP ref_interp;

    REIS(1, pair_pos,
         "required args: --pair donor_mesh.ext reciver_mesh.ext\n");
    if (4 > argc) {
      printf("required args: --pair donor_mesh.ext reciver_mesh.ext\n");
      return REF_FAILURE;
    }

    RSS(ref_mpi_stopwatch_start(ref_mpi), "sw start");
    RSS(ref_part_by_extension(&from, ref_mpi, argv[2]), "import");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "from grid"), "sw start");
    RSS(ref_part_by_extension(&to, ref_mpi, argv[3]), "import");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "to grid"), "sw start");

    RSS(ref_export_tec_surf(to, "ref_interp_test_to.tec"), "export");
    RSS(ref_export_tec_surf(from, "ref_interp_test_from.tec"), "export");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "export viz"), "sw start");

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    ref_interp->instrument = REF_TRUE;
    RSS(ref_interp_locate(ref_interp), "map");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "locate"), "sw start");
    RSS(ref_interp_tec(ref_interp, "ref_interp_test_exhaust.tec"), "export");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "tec"), "sw start");
    RSS(ref_interp_stats(ref_interp), "err");
    RSS(ref_mpi_stopwatch_stop(ref_mpi, "stats"), "sw start");

    RSS(ref_interp_free(ref_interp), "interp free");
    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != rate_pos) {
    REF_GRID grid3, grid2, grid1;
    REF_DBL *solution3, *solution2, *solution1;
    REF_DBL *f3, *f2, *f1, *rate, *diff23, *diff12;
    REF_DBL h3, h2, h1;
    REF_INTERP ref_interp;
    REF_INT ldim, dim;
    REF_INT i, node;
    REIS(1, rate_pos,
         "required args: --rate "
         "coarse_grid.ext coarse_solution.solb "
         "medium_grid.ext medium_solution.solb "
         "fine_grid.ext fine_solution.solb");
    if (8 > argc) {
      printf(
          "required args: --rate "
          "coarse_grid.ext coarse_solution.solb "
          "medium_grid.ext medium_solution.solb "
          "fine_grid.ext fine_solution.solb\n");
      return REF_FAILURE;
    }

    if (ref_mpi_once(ref_mpi)) printf("coarse %s %s\n", argv[2], argv[3]);
    RSS(ref_part_by_extension(&grid3, ref_mpi, argv[2]),
        "part coarse_grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(grid3), &dim, &solution3, argv[3]),
        "unable to load coarse_solution in position 3");
    ldim = dim;

    if (ref_mpi_once(ref_mpi)) printf("medium %s %s\n", argv[4], argv[5]);
    RSS(ref_part_by_extension(&grid2, ref_mpi, argv[4]),
        "part medium_grid in position 4");
    RSS(ref_part_scalar(ref_grid_node(grid2), &dim, &solution2, argv[5]),
        "unable to load medium_solution in position 5");
    REIS(ldim, dim, "expected same ldim");

    if (ref_mpi_once(ref_mpi)) printf("fine %s %s\n", argv[6], argv[7]);
    RSS(ref_part_by_extension(&grid1, ref_mpi, argv[6]),
        "part fine_grid in position 6");
    RSS(ref_part_scalar(ref_grid_node(grid1), &dim, &solution1, argv[7]),
        "unable to load fine_solution in position 7");
    REIS(ldim, dim, "expected same ldim");

    f3 = solution3;
    ref_malloc(f2, ldim * ref_node_max(ref_grid_node(grid3)), REF_DBL);
    ref_malloc(f1, ldim * ref_node_max(ref_grid_node(grid3)), REF_DBL);
    ref_malloc(rate, ldim * ref_node_max(ref_grid_node(grid3)), REF_DBL);
    ref_malloc(diff23, ldim * ref_node_max(ref_grid_node(grid3)), REF_DBL);
    ref_malloc(diff12, ldim * ref_node_max(ref_grid_node(grid3)), REF_DBL);

    if (ref_mpi_once(ref_mpi)) printf("interp medium to coarse\n");
    RSS(ref_interp_create(&ref_interp, grid2, grid3), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    RSS(ref_interp_scalar(ref_interp, ldim, solution2, f2), "interp scalar");
    RSS(ref_interp_free(ref_interp), "interp free");

    if (ref_mpi_once(ref_mpi)) printf("interp fine to coarse\n");
    RSS(ref_interp_create(&ref_interp, grid1, grid3), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    RSS(ref_interp_scalar(ref_interp, ldim, solution1, f1), "interp scalar");
    RSS(ref_interp_free(ref_interp), "interp free");

    h3 = pow((REF_DBL)ref_node_n_global(ref_grid_node(grid3)), -1.0 / 3.0);
    h2 = pow((REF_DBL)ref_node_n_global(ref_grid_node(grid2)), -1.0 / 3.0);
    h1 = pow((REF_DBL)ref_node_n_global(ref_grid_node(grid1)), -1.0 / 3.0);
    h1 = h1 / h3;
    h2 = h2 / h3;
    h3 = 1.0;
    if (ref_mpi_once(ref_mpi)) printf("h3 %f h2 %f h1 %f\n", h3, h2, h1);

    each_ref_node_valid_node(ref_grid_node(grid3), node) {
      for (i = 0; i < ldim; i++) {
        diff23[i + ldim * node] = f2[i + ldim * node] - f3[i + ldim * node];
        diff12[i + ldim * node] = f1[i + ldim * node] - f2[i + ldim * node];
      }
    }

    each_ref_node_valid_node(ref_grid_node(grid3), node) {
      for (i = 0; i < ldim; i++) {
        RSS(ref_interp_convergence_rate(
                f3[i + ldim * node], h3, f2[i + ldim * node], h2,
                f1[i + ldim * node], h1, &(rate[i + ldim * node])),
            "rate");
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("gather ref_interp_rate.tec (rate)\n");
    RSS(ref_gather_scalar_by_extension(grid3, ldim, rate, NULL,
                                       "ref_interp_rate.tec"),
        "gather");
    if (ref_mpi_once(ref_mpi))
      printf("gather ref_interp_diff23.tec (medium-coarse)\n");
    RSS(ref_gather_scalar_by_extension(grid3, ldim, diff23, NULL,
                                       "ref_interp_diff23.tec"),
        "gather");
    if (ref_mpi_once(ref_mpi))
      printf("gather ref_interp_diff12.tec (fine-medium)\n");
    RSS(ref_gather_scalar_by_extension(grid3, ldim, diff12, NULL,
                                       "ref_interp_diff12.tec"),
        "gather");

    ref_free(diff23);
    ref_free(diff12);
    ref_free(rate);
    ref_free(f1);
    ref_free(f2);

    ref_free(solution1);
    ref_grid_free(grid1);
    ref_free(solution2);
    ref_grid_free(grid2);
    ref_free(solution3);
    ref_grid_free(grid3);

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != error_pos) {
    REF_GRID truth_grid, candidate_grid;
    REF_DBL *truth_scalar, *candidate_scalar, *interp_scalar;
    REF_INTERP ref_interp;
    REF_INT p;
    REF_DBL h, error;
    REF_INT ldim;
    REIS(1, error_pos,
         "required args: --error truth_mesh.ext truth_solution.solb "
         "canidate_mesh.ext canidate_solution.solb norm-order\n");
    if (7 > argc) {
      printf(
          "required args: --error truth_mesh.ext truth_solution.solb "
          "canidate_mesh.ext canidate_solution.solb norm-order\n");
      return REF_FAILURE;
    }

    RSS(ref_part_by_extension(&truth_grid, ref_mpi, argv[2]),
        "part truth grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(truth_grid), &ldim, &truth_scalar,
                        argv[3]),
        "unable to load scalar in position 3");
    REIS(1, ldim, "expected one truth scalar");

    RSS(ref_part_by_extension(&candidate_grid, ref_mpi, argv[4]),
        "part candidate grid in position 4");
    RSS(ref_part_scalar(ref_grid_node(candidate_grid), &ldim, &candidate_scalar,
                        argv[5]),
        "unable to load scalar in position 5");
    REIS(1, ldim, "expected one truth scalar");

    p = atoi(argv[6]);

    RSS(ref_interp_create(&ref_interp, candidate_grid, truth_grid),
        "make interp");
    RSS(ref_interp_locate(ref_interp), "map");

    ref_malloc(interp_scalar, ref_node_max(ref_grid_node(truth_grid)), REF_DBL);

    RSS(ref_interp_scalar(ref_interp, 1, candidate_scalar, interp_scalar),
        "interp scalar");
    RSS(ref_interp_integrate(truth_grid, interp_scalar, truth_scalar, p,
                             &error),
        "integrate error");
    h = pow((REF_DBL)ref_node_n_global(ref_grid_node(candidate_grid)),
            (-1.0 / 3.0));
    if (ref_mpi_once(ref_mpi)) {
      printf("%e %e\n", h, error);
    }

    ref_free(interp_scalar);
    RSS(ref_interp_free(ref_interp), "interp free");
    ref_free(candidate_scalar);
    RSS(ref_grid_free(candidate_grid), "free");
    ref_free(truth_scalar);
    RSS(ref_grid_free(truth_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != subset_pos) {
    REF_GRID old_grid, new_grid;
    REF_DBL *old_subset, *new_subset;
    REF_INTERP ref_interp;
    REF_INT ldim;
    REIS(1, subset_pos,
         "required args: --subset old_mesh.ext old_solution.solb "
         "new_mesh.ext new_solution.solb\n");
    if (6 > argc) {
      printf(
          "required args: --subset old_mesh.ext old_solution.solb "
          "new_mesh.ext new_solution.solb\n");
      return REF_FAILURE;
    }

    if (ref_mpi_once(ref_mpi)) printf("read/part old grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&old_grid, ref_mpi, argv[2]),
        "read/part old grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read old grid");
    if (ref_mpi_once(ref_mpi)) printf("read/part old subset %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(old_grid), &ldim, &old_subset, argv[3]),
        "read/part old scalar subset in position 3");
    ref_mpi_stopwatch_stop(ref_mpi, "read old subset");
    if (ref_mpi_once(ref_mpi)) printf("read/part new grid %s\n", argv[4]);
    RSS(ref_part_by_extension(&new_grid, ref_mpi, argv[4]),
        "read/part new grid in position 4");
    if (ref_mpi_once(ref_mpi)) {
      printf("%d leading dim from " REF_GLOB_FMT " old nodes to " REF_GLOB_FMT
             " new nodes\n",
             ldim, ref_node_n_global(ref_grid_node(old_grid)),
             ref_node_n_global(ref_grid_node(new_grid)));
    }
    ref_mpi_stopwatch_stop(ref_mpi, "read new grid");
    RSS(ref_interp_create(&ref_interp, old_grid, new_grid), "make interp");
    RSS(ref_interp_locate_subset(ref_interp), "map");
    ref_mpi_stopwatch_stop(ref_mpi, "locate");

    ref_malloc(new_subset, ldim * ref_node_max(ref_grid_node(new_grid)),
               REF_DBL);

    RSS(ref_interp_scalar(ref_interp, ldim, old_subset, new_subset),
        "interp scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "interp");

    if (ref_mpi_once(ref_mpi)) printf("write/gather new subset %s\n", argv[5]);
    RSS(ref_gather_scalar(new_grid, ldim, new_subset, argv[5]),
        "write/gather new subset");
    ref_mpi_stopwatch_stop(ref_mpi, "write new subset");

    ref_free(new_subset);
    RSS(ref_interp_free(ref_interp), "interp free");
    RSS(ref_grid_free(new_grid), "free");
    ref_free(old_subset);
    RSS(ref_grid_free(old_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != field_pos) {
    REF_GRID old_grid, new_grid;
    REF_DBL *old_field, *new_field;
    REF_INTERP ref_interp;
    REF_INT ldim;
    REIS(1, field_pos,
         "required args: --field old_mesh.ext old_solution.solb "
         "new_mesh.ext new_solution.solb\n");
    if (6 > argc) {
      printf(
          "required args: --field old_mesh.ext old_solution.solb "
          "new_mesh.ext new_solution.solb\n");
      return REF_FAILURE;
    }

    if (ref_mpi_once(ref_mpi)) printf("read/part old grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&old_grid, ref_mpi, argv[2]),
        "read/part old grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read old grid");
    if (ref_mpi_once(ref_mpi)) printf("read/part old field %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(old_grid), &ldim, &old_field, argv[3]),
        "read/part old scalar field in position 3");
    ref_mpi_stopwatch_stop(ref_mpi, "read old field");
    if (ref_mpi_once(ref_mpi)) printf("read/part new grid %s\n", argv[4]);
    RSS(ref_part_by_extension(&new_grid, ref_mpi, argv[4]),
        "read/part new grid in position 4");
    if (ref_mpi_once(ref_mpi)) {
      printf("%d leading dim from " REF_GLOB_FMT " old nodes to " REF_GLOB_FMT
             " new nodes\n",
             ldim, ref_node_n_global(ref_grid_node(old_grid)),
             ref_node_n_global(ref_grid_node(new_grid)));
    }
    ref_mpi_stopwatch_stop(ref_mpi, "read new grid");
    RSS(ref_interp_create(&ref_interp, old_grid, new_grid), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    ref_mpi_stopwatch_stop(ref_mpi, "locate");

    ref_malloc(new_field, ldim * ref_node_max(ref_grid_node(new_grid)),
               REF_DBL);

    RSS(ref_interp_scalar(ref_interp, ldim, old_field, new_field),
        "interp scalar");
    ref_mpi_stopwatch_stop(ref_mpi, "interp");

    if (ref_mpi_once(ref_mpi)) printf("write/gather new field %s\n", argv[5]);
    RSS(ref_gather_scalar(new_grid, ldim, new_field, argv[5]),
        "write/gather new field");
    ref_mpi_stopwatch_stop(ref_mpi, "write new field");

    ref_free(new_field);
    RSS(ref_interp_free(ref_interp), "interp free");
    RSS(ref_grid_free(new_grid), "free");
    ref_free(old_field);
    RSS(ref_grid_free(old_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != mach_pos) {
    REF_GRID ref_grid;
    REF_DBL *field, *mach;
    REF_INT ldim, node;
    REF_DBL gamma = 1.4;

    REIS(1, mach_pos,
         "required args: --mach grid.ext solution.solb mach.solb\n");
    if (5 > argc) {
      printf("required args: --mach grid.ext solution.solb mach.solb\n");
      return REF_FAILURE;
    }

    if (REF_EMPTY != gamma_pos && gamma_pos <= argc - 2) {
      gamma = atof(argv[gamma_pos + 1]);
    }
    if (ref_mpi_once(ref_mpi)) printf("gamma %f\n", gamma);

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "part grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "part grid");

    if (ref_mpi_once(ref_mpi)) printf("reading solution %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load field in position 3");
    ref_mpi_stopwatch_stop(ref_mpi, "part scalar");

    ref_malloc(mach, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, u, v, w, p, temp;
      rho = field[0 + ldim * node];
      u = field[1 + ldim * node];
      v = field[2 + ldim * node];
      w = field[3 + ldim * node];
      p = field[4 + ldim * node];
      temp = gamma * p / rho;
      mach[node] = sqrt((u * u + v * v + w * w) / temp);
    }
    ref_mpi_stopwatch_stop(ref_mpi, "mach");

    if (ref_mpi_once(ref_mpi)) printf("writing mach to %s\n", argv[4]);
    RSS(ref_gather_scalar(ref_grid, 1, mach, argv[4]), "export mach");
    ref_mpi_stopwatch_stop(ref_mpi, "gather scalar");

    ref_free(mach);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != cust_pos) {
    REF_GRID ref_grid;
    REF_DBL *field, *cust;
    REF_INT ldim, node;
    REIS(1, cust_pos,
         "required args: --cust grid.ext solution.solb cust.solb\n");
    if (5 > argc) {
      printf("required args: --cust grid.ext solution.solb cust.solb\n");
      return REF_FAILURE;
    }

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "part grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load field in position 3");

    ref_malloc(cust, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, u, v, w, p, temp;
      REF_DBL gamma = 1.4;
      rho = field[0 + ldim * node];
      u = field[1 + ldim * node];
      v = field[2 + ldim * node];
      w = field[3 + ldim * node];
      p = field[4 + ldim * node];
      temp = gamma * p / rho;
      /* scaled mach and density */
      cust[node] = sqrt((u * u + v * v + w * w) / temp) - 0.2 * rho;
    }

    RSS(ref_gather_scalar(ref_grid, 1, cust, argv[4]), "export cust");

    ref_free(cust);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != entropy_pos) {
    REF_GRID ref_grid;
    REF_DBL *field, *entropy;
    REF_INT ldim, node;
    REIS(1, entropy_pos,
         "required args: --entropy grid.ext solution.solb entropy.solb\n");
    if (5 > argc) {
      printf("required args: --entropy grid.ext solution.solb entropy.solb\n");
      return REF_FAILURE;
    }

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "part grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load field in position 3");

    ref_malloc(entropy, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, p;
      REF_DBL gamma = 1.4;
      rho = field[0 + ldim * node];
      p = field[4 + ldim * node];
      entropy[node] = log(p * gamma / pow(rho, gamma));
    }

    RSS(ref_gather_scalar(ref_grid, 1, entropy, argv[4]), "export entropy");

    ref_free(entropy);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != entropyadj_pos) {
    REF_GRID ref_grid;
    REF_DBL *field, *output;
    REF_INT ldim, odim, node;
    REIS(1, entropyadj_pos,
         "required args: --entropyadj grid.ext solution.solb entropy.solb\n");
    if (5 > argc) {
      printf(
          "required args: --entropyadj grid.ext solution.solb entropy.solb\n");
      return REF_FAILURE;
    }

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "part grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load field in position 3");

    odim = 20;
    ref_malloc(output, odim * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, u, v, w, p, s, e;
      REF_DBL gamma = 1.4;
      rho = field[0 + ldim * node];
      u = field[1 + ldim * node];
      v = field[2 + ldim * node];
      w = field[3 + ldim * node];
      p = field[4 + ldim * node];
      /* entropy adjoint. Is pressure the right nondimensionalization? */
      s = log(p / pow(rho, gamma));
      output[0 + odim * node] =
          (gamma - s) / (gamma - 1.0) - 0.5 * rho * (u * u + v * v + w * w) / p;
      output[1 + odim * node] = rho * u / p;
      output[2 + odim * node] = rho * v / p;
      output[3 + odim * node] = rho * w / p;
      output[4 + odim * node] = -u / p;

      e = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

      output[0 + 5 + odim * node] = rho * u;
      output[1 + 5 + odim * node] = rho * u * u + p;
      output[2 + 5 + odim * node] = rho * u * v;
      output[3 + 5 + odim * node] = rho * u * w;
      output[4 + 5 + odim * node] = u * (e + p);

      output[0 + 10 + odim * node] = rho * v;
      output[1 + 10 + odim * node] = rho * v * u;
      output[2 + 10 + odim * node] = rho * v * v + p;
      output[3 + 10 + odim * node] = rho * v * w;
      output[4 + 10 + odim * node] = v * (e + p);

      output[0 + 15 + odim * node] = rho * w;
      output[1 + 15 + odim * node] = rho * w * u;
      output[2 + 15 + odim * node] = rho * w * v;
      output[3 + 15 + odim * node] = rho * w * w + p;
      output[4 + 15 + odim * node] = w * (e + p);
    }

    RSS(ref_gather_scalar(ref_grid, odim, output, argv[4]), "export output");

    ref_free(output);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (REF_EMPTY != heat_pos) {
    REF_GRID ref_grid;
    REF_DBL *field, *temp, *grad, *heat, *hits;
    REF_INT ldim, node, faceid;
    REIS(1, heat_pos, "required args: --heat grid.ext solution.solb faceid\n");
    if (5 > argc) {
      printf("required args: --heat grid.ext solution.solb faceid\n");
      return REF_FAILURE;
    }

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "part grid in position 2");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field, argv[3]),
        "unable to load field in position 3");
    faceid = atoi(argv[4]);

    ref_malloc(temp, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL rho, p;
      REF_DBL gamma = 1.4;
      rho = field[0 + ldim * node];
      p = field[4 + ldim * node];
      temp[node] = gamma * p / rho;
    }

    ref_malloc(grad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc_init(heat, ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    ref_malloc_init(hits, ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    RSS(ref_recon_gradient(ref_grid, temp, grad, REF_RECON_KEXACT), "grad");
    {
      REF_NODE ref_node = ref_grid_node(ref_grid);
      REF_CELL ref_cell = ref_grid_tri(ref_grid);
      REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
      REF_DBL normal[3], dtdn;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (faceid == nodes[ref_cell_node_per(ref_cell)]) {
          RSS(ref_node_tri_normal(ref_node, nodes, normal), "tri norm");
          RSS(ref_math_normalize(normal), "zero normal tri");
          each_ref_cell_cell_node(ref_cell, cell_node) {
            node = nodes[cell_node];
            dtdn = ref_math_dot(&(grad[3 * node]), normal);
            heat[node] += dtdn;
            hits[node] += 1.0;
          }
        }
      }
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        if (hits[node] > 0.5) {
          heat[node] /= hits[node];
          printf("%.15e %.15e %.15e %.15e\n", ref_node_xyz(ref_node, 0, node),
                 ref_node_xyz(ref_node, 1, node),
                 ref_node_xyz(ref_node, 2, node), heat[node]);
        }
      }
    }

    ref_free(hits);
    ref_free(heat);
    ref_free(grad);
    ref_free(temp);
    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  {
    REF_INTERP ref_interp;
    RSS(ref_interp_setup(&ref_interp, ref_mpi), "setup");
    RSS(ref_interp_teardown(ref_interp), "teardown");
  }

  if (!ref_mpi_para(ref_mpi)) { /* locate between tet2 */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INTERP ref_interp;
    REF_INT node0, node1, new_node;
    REF_GLOB global;
    REF_DBL max_error;
    RSS(ref_fixture_tet2_grid(&ref_grid, ref_mpi), "2 tet fixture");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_grid_cache_background(ref_grid), "cache");
    ref_interp = ref_grid_interp(ref_grid);
    ref_interp_continuously(ref_grid_interp(ref_grid)) = REF_TRUE;
    node0 = 0;
    node1 = 4;
    RSS(ref_node_next_global(ref_node, &global), "next global");
    RSS(ref_node_add(ref_node, global, &new_node), "new node");

    RSS(ref_node_interpolate_edge(ref_node, node0, node1, 0.5, new_node),
        "interp new node");
    RSS(ref_interp_locate_between(ref_interp, node0, node1, new_node),
        "locate");
    RUS(REF_EMPTY, ref_interp_cell(ref_interp, new_node), "not located");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RGDS(1.0e-15, max_error, "large interp error");
    ref_grid_free(ref_grid);
  }

  { /* bricks */
    REF_GRID from, to;
    char file[] = "ref_interp_test.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&from, ref_mpi), "brick");
      RSS(ref_export_by_extension(from, file), "export");
      RSS(ref_grid_free(from), "free");
    }
    RSS(ref_part_by_extension(&from, ref_mpi, file), "import");
    RSS(ref_part_by_extension(&to, ref_mpi, file), "import");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(file), "test clean up");

    RSS(ref_interp_shift_cube_interior(ref_grid_node(to)), "shift");

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(26, ref_interp->n_walk, "walk count");
      REIS(30, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.121 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_interp_create(&ref_interp, to, from), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(26, ref_interp->n_walk, "walk count");
      REIS(30, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.121 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");
  }

  { /* bricks subset */
    REF_GRID from, to;
    char file[] = "ref_interp_test.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&from, ref_mpi), "brick");
      RSS(ref_export_by_extension(from, file), "export");
      RSS(ref_grid_free(from), "free");
    }
    RSS(ref_part_by_extension(&from, ref_mpi, file), "import");
    RSS(ref_part_by_extension(&to, ref_mpi, file), "import");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(file), "test clean up");

    RSS(ref_interp_shift_cube_interior(ref_grid_node(to)), "shift");

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    RSS(ref_interp_locate_subset(ref_interp), "map");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(30, ref_interp->n_walk, "walk count");
      REIS(34, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.121 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_interp_create(&ref_interp, to, from), "make interp");
    RSS(ref_interp_locate_subset(ref_interp), "map");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(30, ref_interp->n_walk, "walk count");
      REIS(34, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.121 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");
  }

  { /* odd split one brick */
    REF_GRID from, to;
    char even[] = "ref_interp_test_even.meshb";
    char odd[] = "ref_interp_test_odd.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if (ref_mpi_once(ref_mpi)) {
      REF_GRID ref_grid;

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, even), "export");
      RSS(ref_grid_free(ref_grid), "free");

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_split_edge_pattern(ref_grid, 1, 2), "split");
      RSS(ref_export_by_extension(ref_grid, odd), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&from, ref_mpi, even), "import");
    RSS(ref_part_by_extension(&to, ref_mpi, odd), "import");
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove(even), "test clean up");
      REIS(0, remove(odd), "test clean up");
    }

    RSS(ref_interp_shift_cube_interior(ref_grid_node(to)), "shift");

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(129, ref_interp->n_walk, "walk count");
      REIS(66, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.241 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(9.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_interp_create(&ref_interp, to, from), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(26, ref_interp->n_walk, "walk count");
      REIS(30, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.241 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");
  }

  { /* odd/even split bricks */
    REF_GRID from, to;
    char even[] = "ref_interp_test_even.meshb";
    char odd[] = "ref_interp_test_odd.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if (ref_mpi_once(ref_mpi)) {
      REF_GRID ref_grid;

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_split_edge_pattern(ref_grid, 0, 2), "split");
      RSS(ref_export_by_extension(ref_grid, even), "export");
      RSS(ref_grid_free(ref_grid), "free");

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_split_edge_pattern(ref_grid, 1, 2), "split");
      RSS(ref_export_by_extension(ref_grid, odd), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&from, ref_mpi, even), "import");
    RSS(ref_part_by_extension(&to, ref_mpi, odd), "import");
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove(even), "test clean up");
      REIS(0, remove(odd), "test clean up");
    }

    RSS(ref_interp_shift_cube_interior(ref_grid_node(to)), "shift");

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(129, ref_interp->n_walk, "walk count");
      REIS(66, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.241 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(9.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_interp_create(&ref_interp, to, from), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");
    REIS(8, ref_interp->n_geom, "geom missing");
    REIS(0, ref_interp->n_geom_fail, "geom fail");
    if (!ref_mpi_para(ref_mpi)) {
      REIS(121, ref_interp->n_walk, "walk count");
      REIS(75, ref_interp->n_tree, "tree count");
    }
    RSS(ref_interp_min_bary(ref_interp, &min_bary), "min bary");
    RAS(-0.241 < min_bary, "large extrapolation");
    RSS(ref_interp_max_error(ref_interp, &max_error), "err");
    RAS(7.0e-16 > max_error, "large interp error");
    RSS(ref_interp_free(ref_interp), "interp free");

    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");
  }

  { /* interp scalar for odd/even split bricks, with curved boundary */
    REF_GRID from, to;
    char even[] = "ref_interp_test_even.meshb";
    char odd[] = "ref_interp_test_odd.meshb";
    REF_INTERP ref_interp;
    REF_DBL *from_scalar, *to_scalar;
    REF_INT node, i;
    REF_DBL dist2;

    if (ref_mpi_once(ref_mpi)) {
      REF_GRID ref_grid;

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_split_edge_pattern(ref_grid, 0, 2), "split");
      RSS(ref_export_by_extension(ref_grid, even), "export");
      RSS(ref_grid_free(ref_grid), "free");

      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_split_edge_pattern(ref_grid, 1, 2), "split");
      RSS(ref_export_by_extension(ref_grid, odd), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&from, ref_mpi, even), "import");
    RSS(ref_interp_shift_cube_interior(ref_grid_node(from)), "shift");
    RSS(ref_part_by_extension(&to, ref_mpi, odd), "import");
    RSS(ref_interp_shift_cube_interior(ref_grid_node(to)), "shift");
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove(even), "test clean up");
      REIS(0, remove(odd), "test clean up");
    }

    ref_malloc(from_scalar, 3 * ref_node_max(ref_grid_node(from)), REF_DBL);
    ref_malloc_init(to_scalar, 3 * ref_node_max(ref_grid_node(to)), REF_DBL,
                    0.0);

    each_ref_node_valid_node(ref_grid_node(from), node) {
      for (i = 0; i < 3; i++) {
        from_scalar[i + 3 * node] = ref_node_xyz(ref_grid_node(from), i, node);
      }
    }

    RSS(ref_interp_create(&ref_interp, from, to), "make interp");
    RSS(ref_interp_locate(ref_interp), "map");

    RSS(ref_interp_scalar(ref_interp, 3, from_scalar, to_scalar), "interp");

    each_ref_node_valid_node(ref_grid_node(to), node) {
      if (!ref_node_owned(ref_grid_node(to), node)) continue;
      dist2 = 0.0;
      for (i = 0; i < 3; i++) {
        dist2 += pow(
            to_scalar[i + 3 * node] - ref_node_xyz(ref_grid_node(to), i, node),
            2);
      }
      RWDS(0.0, dist2, 2.0e-3, "interp scalar xyz not matching");
    }

    RSS(ref_interp_free(ref_interp), "free");
    ref_free(to_scalar);
    ref_free(from_scalar);
    RSS(ref_grid_free(to), "free");
    RSS(ref_grid_free(from), "free");
  }

  { /* integrate scalar */
    char grid[] = "ref_interp_test.meshb";
    REF_GRID ref_grid;
    REF_DBL *truth_scalar, *candidate_scalar;
    REF_INT p;
    REF_DBL error;

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, grid), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, grid), "import");
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove(grid), "test clean up");
    }

    ref_malloc_init(candidate_scalar, ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 3.0);
    ref_malloc_init(truth_scalar, ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 1.0);
    p = 2;
    RSS(ref_interp_integrate(ref_grid, candidate_scalar, truth_scalar, p,
                             &error),
        "int");
    RWDS(2.0, error, -1.0, "expected sqrt(2^2)");

    ref_free(truth_scalar);
    ref_free(candidate_scalar);
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* first order convergence rate uniform */
    REF_DBL f3 = 1.00, f2 = 0.50, f1 = 0.25;
    REF_DBL h3 = 1.00, h2 = 0.50, h1 = 0.25;
    REF_DBL rate;
    RSS(ref_interp_convergence_rate(f3, h3, f2, h2, f1, h1, &rate),
        "conv rate");
    RWDS(1.0, rate, 0.0001, "first order");
  }

  { /* first order convergence rate non-uniform */
    REF_DBL f3 = 1.00, f2 = 0.60, f1 = 0.35;
    REF_DBL h3 = 1.00, h2 = 0.60, h1 = 0.35;
    REF_DBL rate;
    RSS(ref_interp_convergence_rate(f3, h3, f2, h2, f1, h1, &rate),
        "conv rate");
    RWDS(1.0, rate, 0.0001, "first order");
  }

  { /* first order convergence rate non-uniform */
    REF_DBL f3 = 1.00, f2 = 0.89, f1 = 0.81;
    REF_DBL h3 = 1.00, h2 = 0.89, h1 = 0.81;
    REF_DBL rate;
    RSS(ref_interp_convergence_rate(f3, h3, f2, h2, f1, h1, &rate),
        "conv rate");
    RWDS(1.0, rate, 0.0001, "first order");
  }

  { /* second order convergence rate uniform */
    REF_DBL f3 = 1.00, f2 = 0.25, f1 = 0.0625;
    REF_DBL h3 = 1.00, h2 = 0.50, h1 = 0.25;
    REF_DBL rate;
    RSS(ref_interp_convergence_rate(f3, h3, f2, h2, f1, h1, &rate),
        "conv rate");
    RWDS(2.0, rate, 0.0001, "second order");
  }

  { /* second order convergence rate non-uniform */
    REF_DBL f3 = 1.00, f2 = 0.36, f1 = 0.25;
    REF_DBL h3 = 1.00, h2 = 0.60, h1 = 0.50;
    REF_DBL rate;
    RSS(ref_interp_convergence_rate(f3, h3, f2, h2, f1, h1, &rate),
        "conv rate");
    RWDS(2.0, rate, 0.0001, "second order");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
