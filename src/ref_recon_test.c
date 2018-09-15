
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

#include "ref_geom.h"

#include "ref_clump.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

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

    RSS(ref_recon_l2_projection_grad(ref_grid, scalar, grad), "l2 grad");

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

    RSS(ref_recon_l2_projection_hessian(ref_grid, scalar, hessian), "l2 hess");

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

  { /* boundary averaging constant recon */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *recon;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(recon, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      recon[0 + 6 * node] = 100.0;
      recon[1 + 6 * node] = 7.0;
      recon[2 + 6 * node] = 22.0;
      recon[3 + 6 * node] = 200.0;
      recon[4 + 6 * node] = 15.0;
      recon[5 + 6 * node] = 300.0;
    }

    RSS(ref_recon_extrapolate_boundary(recon, ref_grid), "bound extrap");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(100.0, recon[0 + 6 * node], tol, "m11");
      RWDS(7.0, recon[1 + 6 * node], tol, "m12");
      RWDS(22.0, recon[2 + 6 * node], tol, "m13");
      RWDS(200.0, recon[3 + 6 * node], tol, "m22");
      RWDS(15.0, recon[4 + 6 * node], tol, "m23");
      RWDS(300.0, recon[5 + 6 * node], tol, "m33");
    }

    ref_free(recon);

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* multipass boundary averaging constant recon */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *recon;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");

    ref_malloc(recon, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      recon[0 + 6 * node] = 100.0;
      recon[1 + 6 * node] = 7.0;
      recon[2 + 6 * node] = 22.0;
      recon[3 + 6 * node] = 200.0;
      recon[4 + 6 * node] = 15.0;
      recon[5 + 6 * node] = 300.0;
    }

    RSS(ref_recon_extrapolate_boundary_multipass(recon, ref_grid),
        "bound extrap");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RWDS(100.0, recon[0 + 6 * node], tol, "m11");
      RWDS(7.0, recon[1 + 6 * node], tol, "m12");
      RWDS(22.0, recon[2 + 6 * node], tol, "m13");
      RWDS(200.0, recon[3 + 6 * node], tol, "m22");
      RWDS(15.0, recon[4 + 6 * node], tol, "m23");
      RWDS(300.0, recon[5 + 6 * node], tol, "m33");
    }

    ref_free(recon);

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* k-exact for small variation */
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
    RSS(ref_recon_kexact_hessian(ref_grid, scalar, hessian), "k-exact hess");
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
    RSS(ref_recon_kexact_hessian(ref_grid, scalar, hessian), "k-exact hess");
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

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
