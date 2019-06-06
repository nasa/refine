
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

#include "ref_recon.h"

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

#include "ref_cloud.h"

#include "ref_clump.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  if (argc == 3) {
    REF_GRID ref_grid;
    REF_DBL *function, *derivatives, *scalar, *grad;
    REF_INT ldim, node, i, dir;
    REF_RECON_RECONSTRUCTION reconstruction = REF_RECON_L2PROJECTION;
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[1]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[1]),
        "unable to load grid in position 1");

    if (ref_mpi_once(ref_mpi)) printf("reading function %s\n", argv[2]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &function, argv[2]),
        "unable to load function in position 2");
    ref_malloc(derivatives, 3 * ldim * ref_node_max(ref_grid_node(ref_grid)),
               REF_DBL);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(grad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    if (ref_mpi_once(ref_mpi)) printf("reconstruct\n");
    for (i = 0; i < ldim; i++) {
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        scalar[node] = function[i + ldim * node];
      }
      RSS(ref_recon_gradient(ref_grid, scalar, grad, reconstruction), "grad");
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        for (dir = 0; dir < 3; dir++) {
          derivatives[dir + 3 * i + 3 * ldim * node] =
              ABS(grad[dir + 3 * node]);
        }
      }
    }
    if (ref_mpi_once(ref_mpi)) printf("gather %s\n", "ref_recon_deriv.tec");
    RSS(ref_gather_scalar_by_extension(ref_grid, 3 * ldim, derivatives, NULL,
                                       "ref_recon_deriv.tec"),
        "export derivatives");

    ref_free(function);
    ref_free(derivatives);
    ref_free(scalar);
    ref_free(grad);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
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

    RSS(ref_recon_gradient(ref_grid, scalar, grad, REF_RECON_L2PROJECTION),
        "l2 grad");

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

    RSS(ref_recon_hessian(ref_grid, scalar, hessian, REF_RECON_L2PROJECTION),
        "l2 hess");

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

  { /* zeroth order extrapolate boundary averaging constant recon */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_DBL *recon;
    REF_BOOL *replace;
    REF_INT node;
    REF_INT i, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(recon, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(replace, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_INT);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      recon[0 + 6 * node] = 100.0;
      recon[1 + 6 * node] = 7.0;
      recon[2 + 6 * node] = 22.0;
      recon[3 + 6 * node] = 200.0;
      recon[4 + 6 * node] = 15.0;
      recon[5 + 6 * node] = 300.0;
    }
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      each_ref_cell_cell_node(ref_cell, cell_node) {
        for (i = 0; i < 6; i++) recon[i + 6 * nodes[cell_node]] = 0.0;
      }
    }
    RSS(ref_recon_mask_tri(ref_grid, replace, 6), "mask");
    RSS(ref_recon_extrapolate_zeroth(ref_grid, recon, replace, 6),
        "bound extrap");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(100.0, recon[0 + 6 * node], tol, "m11");
      RWDS(7.0, recon[1 + 6 * node], tol, "m12");
      RWDS(22.0, recon[2 + 6 * node], tol, "m13");
      RWDS(200.0, recon[3 + 6 * node], tol, "m22");
      RWDS(15.0, recon[4 + 6 * node], tol, "m23");
      RWDS(300.0, recon[5 + 6 * node], tol, "m33");
    }

    ref_free(replace);
    ref_free(recon);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* seq k-exact gradient for small variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *gradient;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(gradient, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL y = ref_node_xyz(ref_node, 1, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      scalar[node] = 0.5 + 0.01 * x + 0.02 * y + 0.06 * z;
    }
    RSS(ref_recon_gradient(ref_grid, scalar, gradient, REF_RECON_KEXACT),
        "k-exact hess");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.01, gradient[0 + 3 * node], tol, "g[0]");
      RWDS(0.02, gradient[1 + 3 * node], tol, "g[1]");
      RWDS(0.06, gradient[2 + 3 * node], tol, "g[2]");
    }

    ref_free(gradient);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* para file k-exact gradient for small variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *gradient;
    REF_DBL tol = -1.0;
    char file[] = "ref_recon_test.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, file), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "import");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(file), "test clean up");

    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(gradient, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL y = ref_node_xyz(ref_node, 1, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      scalar[node] = 0.5 + 0.01 * x + 0.02 * y + 0.06 * z;
    }
    RSS(ref_recon_gradient(ref_grid, scalar, gradient, REF_RECON_KEXACT),
        "k-exact hess");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.01, gradient[0 + 3 * node], tol, "g[0]");
      RWDS(0.02, gradient[1 + 3 * node], tol, "g[1]");
      RWDS(0.06, gradient[2 + 3 * node], tol, "g[2]");
    }

    ref_free(gradient);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* seq k-exact hessian for small variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *hessian;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(hessian, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL y = ref_node_xyz(ref_node, 1, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      scalar[node] = 0.5 + 0.01 * (0.5 * x * x) + 0.02 * x * y +
                     0.04 * (0.5 * y * y) + 0.06 * (0.5 * z * z);
    }
    RSS(ref_recon_hessian(ref_grid, scalar, hessian, REF_RECON_KEXACT),
        "k-exact hess");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.01, hessian[0 + 6 * node], tol, "m[0]");
      RWDS(0.02, hessian[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, hessian[2 + 6 * node], tol, "m[2]");
      RWDS(0.04, hessian[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, hessian[4 + 6 * node], tol, "m[4]");
      RWDS(0.06, hessian[5 + 6 * node], tol, "m[5]");
    }

    ref_free(hessian);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* para file k-exact hessian for small variation */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *hessian;
    REF_DBL tol = -1.0;
    char file[] = "ref_recon_test.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, file), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "import");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(file), "test clean up");

    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(hessian, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL y = ref_node_xyz(ref_node, 1, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      scalar[node] = 0.5 + 0.01 * (0.5 * x * x) + 0.02 * x * y +
                     0.04 * (0.5 * y * y) + 0.06 * (0.5 * z * z);
    }
    RSS(ref_recon_hessian(ref_grid, scalar, hessian, REF_RECON_KEXACT),
        "k-exact hess");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.01, hessian[0 + 6 * node], tol, "m[0]");
      RWDS(0.02, hessian[1 + 6 * node], tol, "m[1]");
      RWDS(0.00, hessian[2 + 6 * node], tol, "m[2]");
      RWDS(0.04, hessian[3 + 6 * node], tol, "m[3]");
      RWDS(0.00, hessian[4 + 6 * node], tol, "m[4]");
      RWDS(0.06, hessian[5 + 6 * node], tol, "m[5]");
    }

    ref_free(hessian);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* k-exact 2D */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *scalar, *hessian;
    REF_DBL tol = -1.0;

    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(hessian, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      scalar[node] =
          0.5 + 0.01 * (0.5 * x * x) + 0.02 * x * z + 0.06 * (0.5 * z * z);
    }
    RSS(ref_recon_hessian(ref_grid, scalar, hessian, REF_RECON_KEXACT),
        "k-exact hess");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(0.01, hessian[0 + 6 * node], tol, "m[0] xx");
      RWDS(0.00, hessian[1 + 6 * node], tol, "m[1] xy");
      RWDS(0.02, hessian[2 + 6 * node], tol, "m[2] xz");
      RWDS(0.00, hessian[3 + 6 * node], tol, "m[3] yy");
      RWDS(0.00, hessian[4 + 6 * node], tol, "m[4] yz");
      RWDS(0.06, hessian[5 + 6 * node], tol, "m[5] zz");
    }

    ref_free(hessian);
    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* imply recon right tet */
    REF_DBL tol = 1.0e-12;
    REF_GRID ref_grid;
    REF_DBL *hess;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "tet");

    ref_malloc_init(hess, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);

    RSS(ref_recon_roundoff_limit(hess, ref_grid), "imply");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RAS(tol < hess[0 + 6 * node], "m[0]");
      RAS(tol < hess[3 + 6 * node], "m[3]");
      RAS(tol < hess[5 + 6 * node], "m[5]");
    }

    ref_free(hess);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* ghost cloud one_layer cloud */
    REF_NODE ref_node;
    REF_INT local, ghost;
    REF_GLOB global;
    REF_CLOUD one_layer[] = {NULL, NULL};
    REF_DBL xyzs[4];

    RSS(ref_node_create(&ref_node, ref_mpi), "create");

    /* set up 0-local and 1-ghost with global=part */
    global = ref_mpi_rank(ref_mpi);
    RSS(ref_node_add(ref_node, global, &local), "add");
    ref_node_part(ref_node, local) = (REF_INT)global;

    if (ref_mpi_para(ref_mpi)) {
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global = 0;
      RSS(ref_node_add(ref_node, global, &ghost), "add");
      ref_node_part(ref_node, ghost) = (REF_INT)global;
    }

    /* set up local cloud with both nodes if para */
    each_ref_node_valid_node(ref_node, local) {
      RSS(ref_cloud_create(&(one_layer[local]), 4), "cloud storage");
    }
    local = 0;
    xyzs[0] = 10.0 * ref_node_part(ref_node, local);
    xyzs[1] = 20.0 * ref_node_part(ref_node, local);
    xyzs[2] = 30.0 * ref_node_part(ref_node, local);
    xyzs[3] = 40.0 * ref_node_part(ref_node, local);
    global = ref_node_global(ref_node, local);
    RSS(ref_cloud_store(one_layer[0], global, xyzs), "store cloud stencil");
    local = 1;
    if (ref_node_valid(ref_node, local)) {
      xyzs[0] = 10.0 * ref_node_part(ref_node, local);
      xyzs[1] = 20.0 * ref_node_part(ref_node, local);
      xyzs[2] = 30.0 * ref_node_part(ref_node, local);
      xyzs[3] = 40.0 * ref_node_part(ref_node, local);
      global = ref_node_global(ref_node, local);
      RSS(ref_cloud_store(one_layer[0], global, xyzs), "store cloud stencil");
    }
    RSS(ref_recon_ghost_cloud(one_layer, ref_node), "update ghosts");

    if (ref_mpi_para(ref_mpi)) {
      REIS(2, ref_cloud_n(one_layer[0]), "local");
      global = ref_mpi_rank(ref_mpi);
      RAS(ref_cloud_has_global(one_layer[0], global), "local");
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global = 0;
      RAS(ref_cloud_has_global(one_layer[0], global), "local");

      REIS(2, ref_cloud_n(one_layer[1]), "ghost");
      global = ref_mpi_rank(ref_mpi) + 1;
      if (global >= ref_mpi_n(ref_mpi)) global -= ref_mpi_n(ref_mpi);
      RAS(ref_cloud_has_global(one_layer[1], global), "local");
      global = ref_mpi_rank(ref_mpi) + 2;
      if (global >= ref_mpi_n(ref_mpi)) global -= ref_mpi_n(ref_mpi);
      RAS(ref_cloud_has_global(one_layer[1], global), "local");
    } else {
      REIS(1, ref_cloud_n(one_layer[0]), "no one");
    }

    each_ref_node_valid_node(ref_node, local) {
      ref_cloud_free(one_layer[local]); /* no-op for null */
    }
    RSS(ref_node_free(ref_node), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
