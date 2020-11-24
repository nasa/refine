
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

#include "ref_phys.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_args.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_histogram.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_recon.h"

static REF_STATUS ref_phys_tri_grad_nodes(REF_NODE ref_node, REF_INT *nodes,
                                          REF_DBL *scalar, REF_DBL *gradient) {
  REF_DBL area2, dot, side_length;
  REF_DBL grad1[3], grad2[3], edge02[3], edge01[3], norm02[3], norm01[3];
  REF_INT i;
  gradient[0] = 0.0;
  gradient[1] = 0.0;
  gradient[2] = 0.0;

  RSS(ref_node_tri_area(ref_node, nodes, &area2), "area");
  area2 *= 2;

  for (i = 0; i < 3; i++)
    edge01[i] = ref_node_xyz(ref_node, i, nodes[1]) -
                ref_node_xyz(ref_node, i, nodes[0]);
  for (i = 0; i < 3; i++)
    edge02[i] = ref_node_xyz(ref_node, i, nodes[2]) -
                ref_node_xyz(ref_node, i, nodes[0]);

  for (i = 0; i < 3; i++) norm01[i] = edge01[i];
  for (i = 0; i < 3; i++) norm02[i] = edge02[i];
  RSS(ref_math_normalize(norm01), "normalize zero length n0 -> n1");
  RSS(ref_math_normalize(norm02), "normalize zero length n0 -> n2");

  dot = ref_math_dot(edge01, norm02);
  side_length = sqrt(ref_math_dot(edge02, edge02));
  for (i = 0; i < 3; i++) grad1[i] = edge01[i] - dot * norm02[i];
  RSS(ref_math_normalize(grad1), "normalize zero length grad1");
  for (i = 0; i < 3; i++) grad1[i] *= side_length;

  dot = ref_math_dot(edge02, norm01);
  side_length = sqrt(ref_math_dot(edge01, edge01));
  for (i = 0; i < 3; i++) grad2[i] = edge02[i] - dot * norm01[i];
  RSS(ref_math_normalize(grad2), "normalize zero length grad2");
  for (i = 0; i < 3; i++) grad2[i] *= side_length;

  for (i = 0; i < 3; i++)
    gradient[i] =
        (scalar[1] - scalar[0]) * grad1[i] + (scalar[2] - scalar[0]) * grad2[i];

  if (ref_math_divisible(gradient[0], area2) &&
      ref_math_divisible(gradient[1], area2) &&
      ref_math_divisible(gradient[2], area2)) {
    gradient[0] /= area2;
    gradient[1] /= area2;
    gradient[2] /= area2;
  } else {
    gradient[0] = 0.0;
    gradient[1] = 0.0;
    gradient[2] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

static REF_STATUS xy_primitive(REF_INT ldim, REF_DBL *volume, REF_INT node,
                               REF_DBL *primitive) {
  REF_DBL tempu;
  REF_INT i;
  for (i = 0; i < 5; i++) primitive[i] = volume[i + ldim * node];
  tempu = primitive[3];
  primitive[3] = primitive[2];
  primitive[2] = tempu;
  return REF_SUCCESS;
}

static REF_STATUS norm_check_square(REF_DBL x, REF_DBL y, REF_DBL *n) {
  REF_DBL tol = -1.0;

  /* printf("x %f y %f n %f %f\n", x, y, n[0], n[1]); */

  if (x > 0.999 && y < 0.9 && y > 0.1) {
    RWDS(1, n[0], tol, "x=1 nx");
    RWDS(0, n[1], tol, "x=1 ny");
  }

  if (x < 0.001 && y < 0.9 && y > 0.1) {
    RWDS(-1, n[0], tol, "x=0 nx");
    RWDS(0, n[1], tol, "x=0 ny");
  }

  if (y > 0.999 && x < 0.9 && x > 0.1) {
    RWDS(0, n[0], tol, "y=1 nx");
    RWDS(1, n[1], tol, "y=1 ny");
  }

  if (y < 0.001 && x < 0.9 && x > 0.1) {
    RWDS(0, n[0], tol, "y=0 nx");
    RWDS(-1, n[1], tol, "y=0 ny");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_phys_flipper(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri_cell = ref_grid_tri(ref_grid);
  REF_CELL edg_cell = ref_grid_edg(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT right, wrong;
  REF_INT node0, node1, ncell, cell_list[1];
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normal[3];

  right = 0;
  wrong = 0;
  each_ref_cell_valid_cell_with_nodes(tri_cell, cell, nodes) {
    RSS(ref_node_tri_normal(ref_node, nodes, normal), "norm inside of area");
    if (normal[2] < 0.0) {
      right++;
    } else {
      ref_cell_c2n(tri_cell, 0, cell) = nodes[1];
      ref_cell_c2n(tri_cell, 1, cell) = nodes[0];
      wrong++;
    }
  }
  printf("tri %d right %d wrong (per EGADS)\n", right, wrong);

  right = 0;
  wrong = 0;
  each_ref_cell_valid_cell_with_nodes(edg_cell, cell, nodes) {
    node0 = nodes[0];
    node1 = nodes[1];
    RSS(ref_cell_list_with2(tri_cell, node0, node1, 1, &ncell, cell_list),
        "found more than one tri with two nodes");
    REIS(ncell, 1, "boundry one tri with two nodes");
    RSS(ref_cell_nodes(tri_cell, cell_list[0], tri_nodes), "tri nodes");
    if ((node1 == tri_nodes[0] && node0 == tri_nodes[1]) ||
        (node1 == tri_nodes[1] && node0 == tri_nodes[2]) ||
        (node1 == tri_nodes[2] && node0 == tri_nodes[0])) {
      right++;
    } else {
      ref_cell_c2n(edg_cell, 0, cell) = node1;
      ref_cell_c2n(edg_cell, 1, cell) = node0;
      wrong++;
    }
  }
  printf("edg %d right %d wrong (per EGADS)\n", right, wrong);

  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_INT laminar_flux_pos = REF_EMPTY;
  REF_INT euler_flux_pos = REF_EMPTY;
  REF_INT convdiff_flux_pos = REF_EMPTY;
  REF_INT mask_pos = REF_EMPTY;
  REF_INT cont_res_pos = REF_EMPTY;
  REF_INT uplus_pos = REF_EMPTY;
  REF_INT turb1_pos = REF_EMPTY;
  REF_INT entropy_output_pos = REF_EMPTY;

  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--laminar-flux", &laminar_flux_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--euler-flux", &euler_flux_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--convdiff", &convdiff_flux_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--mask", &mask_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--cont-res", &cont_res_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--uplus", &uplus_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--turb1", &turb1_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--entropy-output", &entropy_output_pos),
      REF_NOT_FOUND, "arg search");

  if (laminar_flux_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL mach, re, temperature;
    REF_DBL *primitive_dual, *dual_flux;
    REF_RECON_RECONSTRUCTION recon = REF_RECON_L2PROJECTION;
    REF_DBL *prim, *grad, *onegrad;
    REF_INT ldim;
    REF_INT node, i, dir;
    REF_DBL direction[3], state[5], flux[5], gradient[15];
    REF_DBL turb = -1.0;

    REIS(1, laminar_flux_pos,
         "required args: --laminar-flux grid.meshb primitive_dual.solb Mach Re "
         "Temperature(Kelvin) dual_flux.solb");
    if (8 > argc) {
      printf(
          "required args: --laminar-flux grid.meshb primitive_dual.solb Mach "
          "Re Temperature(Kelvin) dual_flux.solb\n");
      return REF_FAILURE;
    }
    mach = atof(argv[4]);
    re = atof(argv[5]);
    temperature = atof(argv[6]);
    if (ref_mpi_once(ref_mpi))
      printf("Reference Mach %f Re %e temperature %f\n", mach, re, temperature);

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    ref_malloc(dual_flux, 20 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    if (ref_mpi_once(ref_mpi)) printf("reading primitive_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &primitive_dual,
                        argv[3]),
        "unable to load primitive_dual in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected 10 (rho,u,v,w,p,5*adj) or 12 (rho,u,v,w,p,turb,6*adj)");

    if (ref_mpi_once(ref_mpi)) printf("copy dual\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 5; i++) {
        dual_flux[i + 20 * node] = primitive_dual[ldim / 2 + i + ldim * node];
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("zero flux\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 15; i++) {
        dual_flux[5 + i + 20 * node] = 0.0;
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("Euler flux\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 5; i++) {
        state[i] = primitive_dual[i + ldim * node];
      }
      for (dir = 0; dir < 3; dir++) {
        direction[0] = 0;
        direction[1] = 0;
        direction[2] = 0;
        direction[dir] = 1;
        RSS(ref_phys_euler(state, direction, flux), "euler");
        for (i = 0; i < 5; i++) {
          dual_flux[i + 5 + 5 * dir + 20 * node] += flux[i];
        }
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("reconstruct gradient\n");
    ref_malloc(grad, 15 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(prim, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(onegrad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    for (i = 0; i < 5; i++) {
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        prim[node] = primitive_dual[i + ldim * node];
      }
      RSS(ref_recon_gradient(ref_grid, prim, onegrad, recon), "grad_lam");
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        for (dir = 0; dir < 3; dir++) {
          grad[dir + 3 * i + 15 * node] = onegrad[dir + 3 * node];
        }
      }
    }
    ref_free(onegrad);
    ref_free(prim);

    if (ref_mpi_once(ref_mpi)) printf("Laminar flux\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 5; i++) {
        state[i] = primitive_dual[i + ldim * node];
      }
      if (6 == ldim) {
        turb = primitive_dual[5 + ldim * node];
      } else {
        turb = -1.0; /* laminar */
      }
      for (i = 0; i < 15; i++) {
        gradient[i] = grad[i + 15 * node];
      }
      for (dir = 0; dir < 3; dir++) {
        direction[0] = 0;
        direction[1] = 0;
        direction[2] = 0;
        direction[dir] = 1;
        RSS(ref_phys_viscous(state, gradient, turb, mach, re, temperature,
                             direction, flux),
            "laminar");
        for (i = 0; i < 5; i++) {
          dual_flux[i + 5 + 5 * dir + 20 * node] -= flux[i];
        }
      }
    }
    ref_free(grad);
    ref_free(primitive_dual);

    if (ref_mpi_once(ref_mpi)) printf("writing dual_flux %s\n", argv[7]);
    RSS(ref_gather_scalar_by_extension(ref_grid, 20, dual_flux, NULL, argv[7]),
        "export dual_flux");

    ref_free(dual_flux);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (euler_flux_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *primitive_dual, *dual_flux;
    REF_INT ldim;

    REIS(1, euler_flux_pos,
         "required args: --euler-flux grid.meshb primitive_dual.solb "
         "dual_flux.solb");
    if (5 > argc) {
      printf(
          "required args: --euler-flux grid.meshb primitive_dual.solb "
          "dual_flux.solb\n");
      return REF_FAILURE;
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    ref_malloc(dual_flux, 20 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    if (ref_mpi_once(ref_mpi)) printf("reading primitive_dual %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &primitive_dual,
                        argv[3]),
        "unable to load primitive_dual in position 3");
    RAS(10 == ldim || 12 == ldim,
        "expected 10 (rho,u,v,w,p,5*adj) or 12 primitive_dual");

    if (ref_mpi_once(ref_mpi)) printf("form dual_flux\n");
    RSS(ref_phys_euler_dual_flux(ref_grid, ldim, primitive_dual, dual_flux),
        "euler dual_flux");

    if (ref_mpi_once(ref_mpi)) printf("writing dual_flux %s\n", argv[4]);
    RSS(ref_gather_scalar_by_extension(ref_grid, 20, dual_flux, NULL, argv[4]),
        "export dual_flux");

    ref_free(dual_flux);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (convdiff_flux_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *primitive, *dual, *dual_flux;
    REF_DBL *grad, flux[1], direction[3], diffusivity;
    REF_INT ldim;
    REF_INT node, i, dir;
    REF_RECON_RECONSTRUCTION recon = REF_RECON_L2PROJECTION;

    REIS(1, convdiff_flux_pos,
         "required args: --convdiff-flux grid.meshb primitive.solb "
         "dual.solb diffusivity "
         "dual_flux.solb");
    if (7 > argc) {
      printf(
          "required args: --convdiff-flux grid.meshb primitive.solb "
          "dual.solb diffusivity "
          "dual_flux.solb\n");
      return REF_FAILURE;
    }

    diffusivity = atof(argv[5]);
    if (ref_mpi_once(ref_mpi)) printf("diffusivity %f\n", diffusivity);

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    ref_malloc(dual_flux, 4 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    if (ref_mpi_once(ref_mpi)) printf("reading primitive %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &primitive, argv[3]),
        "unable to load primitive in position 3");
    REIS(1, ldim, "expected 1 (s) primitive");

    if (ref_mpi_once(ref_mpi)) printf("reading dual %s\n", argv[4]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &dual, argv[4]),
        "unable to load primitive in position 3");
    REIS(1, ldim, "expected 1 (lambda) dual");

    if (ref_mpi_once(ref_mpi)) printf("copy dual\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      i = 0;
      dual_flux[i + 4 * node] = dual[i + 1 * node];
    }
    if (ref_mpi_once(ref_mpi)) printf("zero flux\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 1; i < 4; i++) {
        dual_flux[i + 4 * node] = 0.0;
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("reconstruct gradient\n");
    ref_malloc(grad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_recon_gradient(ref_grid, primitive, grad, recon), "grad_lam");

    if (ref_mpi_once(ref_mpi)) printf("add convdiff flux\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (dir = 0; dir < 3; dir++) {
        direction[0] = 0;
        direction[1] = 0;
        direction[2] = 0;
        direction[dir] = 1;
        RSS(ref_phys_convdiff(&(primitive[node]), &(grad[3 * node]),
                              diffusivity, direction, flux),
            "convdiff");
        dual_flux[1 + dir + 4 * node] += flux[0];
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("writing dual_flux %s\n", argv[6]);
    RSS(ref_gather_scalar_by_extension(ref_grid, 4, dual_flux, NULL, argv[6]),
        "export dual_flux");

    ref_free(grad);
    ref_free(dual_flux);
    ref_free(dual);
    ref_free(primitive);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (mask_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DICT ref_dict;
    REF_BOOL *replace;
    REF_DBL *primitive_dual;
    REF_INT ldim;

    REIS(1, mask_pos,
         "required args: --mask grid.meshb grid.mapbc primitive_dual.solb "
         "strong_replacement.solb");
    if (6 > argc) {
      printf(
          "required args: --mask grid.meshb grid.mapbc primitive_dual.solb "
          "strong_replacement.solb\n");
      return REF_FAILURE;
    }

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");

    if (ref_mpi_once(ref_mpi)) printf("reading bc map %s\n", argv[3]);
    RSS(ref_dict_create(&ref_dict), "create");
    RSS(ref_phys_read_mapbc(ref_dict, argv[3]),
        "unable to mapbc in position 3");

    if (ref_mpi_once(ref_mpi)) printf("reading primitive_dual %s\n", argv[4]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &primitive_dual,
                        argv[4]),
        "unable to load primitive_dual in position 3");
    RAS(10 == ldim || 12 == ldim || 15 == ldim || 18 == ldim,
        "expected rho,u,v,w,p,5*adj,[5*dfdq] or "
        "rho,u,v,w,p,turb,6*adj,[6*dfdq]");

    ref_malloc(replace, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_BOOL);
    RSS(ref_phys_mask_strong_bcs(ref_grid, ref_dict, replace, ldim), "mask");
    RSS(ref_recon_extrapolate_zeroth(ref_grid, primitive_dual, replace, ldim),
        "extrapolate zeroth order");
    ref_free(replace);

    if (ref_mpi_once(ref_mpi))
      printf("writing strong_replacement %s\n", argv[5]);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, primitive_dual, NULL,
                                       argv[5]),
        "export primitive_dual");

    ref_free(primitive_dual);
    RSS(ref_dict_free(ref_dict), "free");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (cont_res_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;
    REF_DBL *dual_flux, *system, *flux, *res, *weight, *packed_weight;
    REF_INT ldim, n;
    REF_INT equ, dir, node;
    REF_DBL convergence_rate, exponent, total, l2res;
    REF_DBL min_weight, max_weight, median, minmax;
    REF_INT nsystem, nequ;

    REIS(1, cont_res_pos,
         "required args: --cont-res grid.meshb dual_flux.solb convergence_rate "
         "weight.solb");
    if (6 > argc) {
      printf(
          "required args: --cont-res grid.meshb dual_flux.solb "
          "convergence_rate exponent weight.solb\n");
      return REF_FAILURE;
    }

    convergence_rate = atof(argv[4]);
    exponent = 0.25;
    if (convergence_rate > 0.0) exponent = -1.0 / convergence_rate;
    if (ref_mpi_once(ref_mpi))
      printf("convergence rate %f exponent %f\n", convergence_rate, exponent);

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);

    if (ref_mpi_once(ref_mpi)) printf("reading dual flux %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &dual_flux, argv[3]),
        "unable to load scalar in position 3");
    RAS(4 == ldim || 20 == ldim,
        "expected 4 (adj,xflux,yflux,zflux) "
        "or 20 (5*adj,5*xflux,5*yflux,5*zflux)");
    nequ = ldim / 4;
    nsystem = 1 + nequ + nequ;
    ref_malloc_init(weight, ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc_init(flux, 3 * nequ * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 0.0);
    ref_malloc_init(res, nequ * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc_init(system, nsystem * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 0.0);

    if (ref_mpi_once(ref_mpi)) printf("compute residual\n");
    each_ref_node_valid_node(ref_node, node) {
      for (dir = 0; dir < 3; dir++) {
        for (equ = 0; equ < nequ; equ++) {
          flux[equ + dir * nequ + 3 * nequ * node] =
              dual_flux[equ + dir * nequ + nequ + ldim * node];
        }
      }
    }
    RSS(ref_phys_cc_fv_embed(ref_grid, nequ, flux, res), "res");
    each_ref_node_valid_node(ref_node, node) {
      for (equ = 0; equ < nequ; equ++) {
        system[equ + nsystem * node] = res[equ + nequ * node];
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("copy dual to system\n");
    each_ref_node_valid_node(ref_node, node) {
      for (equ = 0; equ < nequ; equ++) {
        system[equ + nequ + nsystem * node] = dual_flux[equ + ldim * node];
      }
    }

    if (ref_mpi_once(ref_mpi)) printf("compute weight\n");
    total = 0;
    l2res = 0;
    min_weight = 1.0e300;
    max_weight = -1.0e300;
    each_ref_node_valid_node(ref_node, node) {
      for (equ = 0; equ < nequ; equ++) {
        weight[node] +=
            ABS(dual_flux[equ + ldim * node] * res[equ + nequ * node]);
        if (ref_node_owned(ref_node, node))
          l2res += res[equ + nequ * node] * res[equ + nequ * node];
      }
      /* approximate boundary with double weight */
      if (!ref_cell_node_empty(ref_cell, node)) weight[node] *= 2.0;
    }

    ref_malloc_init(packed_weight, ref_node_max(ref_node), REF_DBL, 0.0);
    n = 0;
    each_ref_node_valid_node(ref_node, node) {
      packed_weight[n] = weight[node];
      n++;
    }
    REIS(n, ref_node_n(ref_node), "node miscount");
    RSS(ref_search_selection(ref_mpi, n, packed_weight,
                             ref_node_n_global(ref_node) / 2, &median),
        "parallel median selection");
    each_ref_node_valid_node(ref_node, node) {
      weight[node] /= median;
      /* weight expected in length scale */
      if (weight[node] > 0.0) weight[node] = pow(weight[node], exponent);
      if (ref_node_owned(ref_node, node)) {
        total += weight[node];
        min_weight = MIN(min_weight, weight[node]);
        max_weight = MAX(max_weight, weight[node]);
      }
      system[nsystem - 1 + nsystem * node] = weight[node];
    }

    RSS(ref_mpi_allsum(ref_mpi, &total, 1, REF_DBL_TYPE), "sum total");
    total /= (REF_DBL)ref_node_n_global(ref_node);

    minmax = min_weight;
    RSS(ref_mpi_min(ref_mpi, &minmax, &min_weight, REF_DBL_TYPE), "mpi min");
    RSS(ref_mpi_bcast(ref_mpi, &min_weight, 1, REF_DBL_TYPE), "mbast");

    minmax = max_weight;
    RSS(ref_mpi_max(ref_mpi, &minmax, &max_weight, REF_DBL_TYPE), "mpi max");
    RSS(ref_mpi_bcast(ref_mpi, &max_weight, 1, REF_DBL_TYPE), "mbast");

    RSS(ref_mpi_allsum(ref_mpi, &l2res, 1, REF_DBL_TYPE), "sum l2res");
    l2res /= (REF_DBL)ref_node_n_global(ref_node);
    l2res = sqrt(l2res);
    if (ref_mpi_once(ref_mpi)) {
      printf("median %e\n", median);
      printf("L2 res %e\n", l2res);
      printf("L1 h scale weight %e\n", total);
      printf("min max %e %e\n", min_weight, max_weight);
    }
    if (ref_mpi_once(ref_mpi)) printf("writing weight %s\n", argv[5]);
    RSS(ref_gather_scalar_by_extension(ref_grid, 1, weight, NULL, argv[5]),
        "export weight");
    if (ref_mpi_once(ref_mpi))
      printf("writing res,dual,weight ref_phys_system.tec\n");
    RSS(ref_gather_scalar_by_extension(ref_grid, nsystem, system, NULL,
                                       "ref_phys_system.tec"),
        "export primitive_dual");
    if (ref_mpi_once(ref_mpi)) printf("writing histogram.tec\n");
    RSS(ref_histogram_node_tec(ref_grid, weight), "export histogram");

    ref_free(weight);
    ref_free(flux);
    ref_free(res);
    ref_free(system);
    ref_free(dual_flux);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (uplus_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *scalar, yplus, uplus;
    REF_INT ldim, node;
    ref_mpi_stopwatch_start(ref_mpi);
    REIS(1, uplus_pos,
         "required args: --uplus grid.meshb dist.solb [yplus=1] uplus.solb");
    if (6 > argc) {
      printf(
          "required args: --uplus grid.meshb dist.solb [yplus=1] uplus.solb");
      return REF_FAILURE;
    }
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    if (ref_mpi_once(ref_mpi)) printf("reading distance %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, argv[3]),
        "unable to load distance in position 3");
    REIS(1, ldim, "expected one distance");
    ref_mpi_stopwatch_stop(ref_mpi, "read distance");

    yplus = atof(argv[4]);
    if (ref_mpi_once(ref_mpi)) printf("yplus=1 of %f\n", yplus);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RAS(ref_math_divisible(scalar[node], yplus),
          "distance not divisible by yplus=1")
      RSS(ref_phys_spalding_uplus(scalar[node] / yplus, &uplus), "uplus");
      scalar[node] = uplus;
    }
    ref_mpi_stopwatch_stop(ref_mpi, "compute uplus");

    if (ref_mpi_once(ref_mpi)) printf("writing uplus %s\n", argv[5]);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, NULL, argv[5]),
        "export uplus");
    ref_mpi_stopwatch_stop(ref_mpi, "write uplus");

    ref_free(scalar);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (turb1_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *distance, *field, turb1;
    REF_INT ldim, node;
    ref_mpi_stopwatch_start(ref_mpi);
    REIS(1, turb1_pos,
         "required args: --turb1 grid.meshb dist.solb volume.solb");
    if (5 > argc) {
      printf("required args: --turb1 grid.meshb dist.solb volume.solb");
      return REF_FAILURE;
    }
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    if (ref_mpi_once(ref_mpi)) printf("reading distance %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &distance, argv[3]),
        "unable to load distance in position 3");
    REIS(1, ldim, "expected one distance");
    ref_mpi_stopwatch_stop(ref_mpi, "read distance");

    ldim = 6;
    ref_malloc(field, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RSS(ref_phys_sa_surrogate(distance[node], &turb1), "sa soln");
      field[0 + ldim * node] = 1.0;
      field[1 + ldim * node] = 0.5;
      field[2 + ldim * node] = 0.0;
      field[3 + ldim * node] = 0.0;
      field[4 + ldim * node] = 1.0 / 1.4;
      field[5 + ldim * node] = turb1;
    }
    ref_mpi_stopwatch_stop(ref_mpi, "compute turb1 field");

    if (ref_mpi_once(ref_mpi)) printf("writing field %s\n", argv[4]);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, field, NULL, argv[4]),
        "export turb1");
    ref_mpi_stopwatch_stop(ref_mpi, "write field");

    ref_free(field);
    ref_free(distance);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  /* src/pde/NS/OutputEuler2D.h src/pde/NS/OutputNavierStokes2D.h */
  if (entropy_output_pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_DBL *volume;
    REF_INT ldim;
    REF_DBL primitive[5], primitive0[5], primitive1[5], primitive2[5];
    REF_DBL dual[5], gradient[15], tri_grad[3], scalar[3];
    REF_DBL flux[3], laminar_flux[5];
    REF_CELL ref_cell;
    REF_NODE ref_node;
    REF_INT i, cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL inviscid_total, viscous_boundary, viscous_interior;
    REF_DBL area, dx[3], normal[3];
    REF_INT part;
    REF_DBL *grad, *prim, *onegrad;
    REF_INT dir, node;
    REF_RECON_RECONSTRUCTION recon = REF_RECON_L2PROJECTION;
    REF_DBL mach = 0.5, re = 100.0, temperature = 288.15, turb = -1.0;

    ref_mpi_stopwatch_start(ref_mpi);
    REIS(1, entropy_output_pos,
         "required args: --entropy-output grid.meshb volume.solb");
    if (4 > argc) {
      printf("required args: --entropy-output grid.meshb volume.solb");
      return REF_FAILURE;
    }
    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_mpi_stopwatch_stop(ref_mpi, "read grid");

    RSS(ref_phys_flipper(ref_grid), "flip it");

    if (ref_mpi_once(ref_mpi)) printf("reading volume %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &volume, argv[3]),
        "unable to load volume in position 3");
    RAS(ldim >= 5, "expected a ldim of at least 5");
    ref_mpi_stopwatch_stop(ref_mpi, "read volume");

    inviscid_total = 0.0;
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_edg(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) != part) continue;
      if (1 == nodes[2] || 4 == nodes[2]) continue;
      RSS(xy_primitive(ldim, volume, nodes[0], primitive0), "prim 0");
      RSS(xy_primitive(ldim, volume, nodes[1], primitive1), "prim 1");
      for (i = 0; i < 5; i++)
        primitive[i] = 0.5 * (primitive0[i] + primitive1[i]);
      RSS(ref_phys_entropy_flux(primitive, flux), "flux1");
      for (i = 0; i < 3; i++)
        dx[i] = ref_node_xyz(ref_node, i, nodes[1]) -
                ref_node_xyz(ref_node, i, nodes[0]);
      normal[0] = dx[1];
      normal[1] = -dx[0];
      normal[2] = 0.0;
      area = sqrt(ref_math_dot(dx, dx));
      RSS(ref_math_normalize(normal), "norm");
      RSS(norm_check_square(ref_node_xyz(ref_node, 0, nodes[0]),
                            ref_node_xyz(ref_node, 1, nodes[0]), normal),
          "square norm");
      inviscid_total += area * ref_math_dot(normal, flux);
    }
    RSS(ref_mpi_allsum(ref_mpi, &inviscid_total, 1, REF_DBL_TYPE), "mpi sum");
    if (ref_mpi_once(ref_mpi))
      printf("inviscid total   = %9.6f\n", inviscid_total);

    if (ref_mpi_once(ref_mpi)) printf("reconstruct gradient\n");
    ref_malloc(grad, 15 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(prim, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    ref_malloc(onegrad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    for (i = 0; i < 5; i++) {
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        RSS(xy_primitive(ldim, volume, node, primitive), "prim 0");
        prim[node] = primitive[i];
      }
      RSS(ref_recon_gradient(ref_grid, prim, onegrad, recon), "grad_lam");
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        for (dir = 0; dir < 3; dir++) {
          grad[dir + 3 * i + 15 * node] = onegrad[dir + 3 * node];
        }
      }
    }
    ref_free(onegrad);
    ref_free(prim);

    viscous_boundary = 0.0;
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_edg(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) != part) continue;
      if (1 == nodes[2] || 4 == nodes[2]) continue;
      RSS(xy_primitive(ldim, volume, nodes[0], primitive0), "prim 0");
      RSS(xy_primitive(ldim, volume, nodes[1], primitive1), "prim 1");
      for (i = 0; i < 5; i++)
        primitive[i] = 0.5 * (primitive0[i] + primitive1[i]);
      /* dual is same as v, symmetric */
      RSS(ref_phys_entropy_adjoint(primitive, dual), "flux1");

      for (i = 0; i < 3; i++)
        dx[i] = ref_node_xyz(ref_node, i, nodes[1]) -
                ref_node_xyz(ref_node, i, nodes[0]);
      normal[0] = dx[1];
      normal[1] = -dx[0];
      normal[2] = 0.0;
      area = sqrt(ref_math_dot(dx, dx));
      RSS(ref_math_normalize(normal), "norm");
      for (i = 0; i < 15; i++)
        gradient[i] = 0.5 * (grad[i + 15 * nodes[0]] + grad[i + 15 * nodes[1]]);
      RSS(ref_phys_viscous(primitive, gradient, turb, mach, re, temperature,
                           normal, laminar_flux),
          "laminar");
      for (i = 0; i < 5; i++)
        viscous_boundary -= area * dual[i] * laminar_flux[i];
    }
    RSS(ref_mpi_allsum(ref_mpi, &viscous_boundary, 1, REF_DBL_TYPE), "mpi sum");
    if (ref_mpi_once(ref_mpi))
      printf("viscous boundary = %9.6f\n", viscous_boundary);

    ref_free(grad);

    viscous_interior = 0.0;
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) != part) continue;
      RSS(xy_primitive(ldim, volume, nodes[0], primitive0), "prim 0");
      RSS(xy_primitive(ldim, volume, nodes[1], primitive1), "prim 1");
      RSS(xy_primitive(ldim, volume, nodes[2], primitive2), "prim 1");
      for (i = 0; i < 5; i++)
        primitive[i] =
            (1.0 / 3.0) * (primitive0[i] + primitive1[i] + primitive2[i]);

      RSS(ref_node_tri_area(ref_node, nodes, &area), "area");

      for (dir = 0; dir < 3; dir++) {
        normal[0] = 0;
        normal[1] = 0;
        normal[2] = 0;
        normal[dir] = 1;

        for (i = 0; i < 5; i++) {
          scalar[0] = primitive0[i];
          scalar[1] = primitive1[i];
          scalar[2] = primitive2[i];
          RSS(ref_phys_tri_grad_nodes(ref_node, nodes, scalar, tri_grad), "tg");
          gradient[0 + 3 * i] = tri_grad[0];
          gradient[1 + 3 * i] = tri_grad[1];
          gradient[2 + 3 * i] = tri_grad[2];
        }
        RSS(ref_phys_viscous(primitive, gradient, turb, mach, re, temperature,
                             normal, laminar_flux),
            "laminar");
        for (i = 0; i < 5; i++) {
          RSS(ref_phys_entropy_adjoint(primitive0, dual), "flux1");
          scalar[0] = dual[i];
          RSS(ref_phys_entropy_adjoint(primitive1, dual), "flux1");
          scalar[1] = dual[i];
          RSS(ref_phys_entropy_adjoint(primitive2, dual), "flux1");
          scalar[2] = dual[i];
          RSS(ref_phys_tri_grad_nodes(ref_node, nodes, scalar, tri_grad), "tg");
          gradient[0 + 3 * i] = tri_grad[0];
          gradient[1 + 3 * i] = tri_grad[1];
          gradient[2 + 3 * i] = tri_grad[2];
        }
        for (i = 0; i < 5; i++)
          viscous_interior += area * gradient[dir + 3 * i] * laminar_flux[i];
      }
    }
    RSS(ref_mpi_allsum(ref_mpi, &viscous_interior, 1, REF_DBL_TYPE), "mpi sum");
    if (ref_mpi_once(ref_mpi))
      printf("viscous interior = %9.6f\n", viscous_interior);

    if (ref_mpi_once(ref_mpi))
      printf("total            = %9.6f\n",
             (inviscid_total + viscous_boundary + viscous_interior));
    if (ref_mpi_once(ref_mpi))
      printf("%d %e %e %e %e %e # conv\n", (REF_INT)ref_node_n_global(ref_node),
             pow((REF_DBL)ref_node_n_global(ref_node), (-1.0 / 2.0)),
             inviscid_total, viscous_boundary, viscous_interior,
             (inviscid_total + viscous_boundary + viscous_interior));

    ref_free(volume);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* converts primitive and conserved */
    REF_DBL state[5];
    REF_DBL primitive[5];
    REF_DBL conserved[5];
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    RSS(ref_phys_make_conserved(state, conserved), "cons");
    RSS(ref_phys_make_primitive(conserved, primitive), "prim");
    RWDS(state[0], primitive[0], -1, "rho");
    RWDS(state[1], primitive[1], -1, "u");
    RWDS(state[2], primitive[2], -1, "v");
    RWDS(state[3], primitive[3], -1, "w");
    RWDS(state[4], primitive[4], -1, "p");
  }

  { /* entropy adjoint */
    REF_DBL primitive[5];
    REF_DBL dual[5];
    primitive[0] = 1.0;
    primitive[1] = 0.5;
    primitive[2] = 0.1;
    primitive[3] = 0.2;
    primitive[4] = 0.8 / 1.4;
    RSS(ref_phys_entropy_adjoint(primitive, dual), "entropy adj");
    RWDS(4.636539469838557, dual[0], -1, "cont");
    RWDS(0.875, dual[1], -1, "x-mom");
    RWDS(0.175, dual[2], -1, "y-mom");
    RWDS(0.35, dual[3], -1, "z-mom");
    RWDS(-1.75, dual[4], -1, "energy");
  }

  { /* x-Euler flux */
    REF_DBL state[5], direction[3];
    REF_DBL flux[5];
    state[0] = 1.0;
    state[1] = 0.2;
    state[2] = 0.0;
    state[3] = 0.0;
    state[4] = 1.0 / 1.4;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    RSS(ref_phys_euler(state, direction, flux), "euler");
    RWDS(0.2, flux[0], -1, "mass flux");
    RWDS(0.04 + 1.0 / 1.4, flux[1], -1, "x mo flux");
    RWDS(0.0, flux[2], -1, "y mo flux");
    RWDS(0.0, flux[3], -1, "z mo flux");
    RWDS(0.504, flux[4], -1, "energy flux");
  }

  { /* dEuler/drho */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], pert[5], cons_offset[5], prim_offset[5], direction[3];
    REF_DBL jac[5], fp[5], fm[5];
    REF_DBL step = 1.0e-8;
    REF_DBL tol = 1.e-8;
    REF_INT i;
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    direction[0] = 1.0 / 9.0;
    direction[1] = 4.0 / 9.0;
    direction[2] = 8.0 / 9.0;

    pert[0] = step;
    pert[1] = 0.0;
    pert[2] = 0.0;
    pert[3] = 0.0;
    pert[4] = 0.0;

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] += pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fp), "euler");

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] -= pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fm), "euler");

    for (i = 0; i < 5; i++) jac[i] = (fp[i] - fm[i]) / (2.0 * step);

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");

    i = 0;
    RWDS(jac[0], dflux_dcons[0 + i * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + i * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + i * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + i * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + i * 5], tol, "energy flux");
  }

  { /* dEuler/dmo-x */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], pert[5], cons_offset[5], prim_offset[5], direction[3];
    REF_DBL jac[5], fp[5], fm[5];
    REF_DBL step = 1.0e-8;
    REF_DBL tol = 1.e-8;
    REF_INT i;
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    direction[0] = 1.0 / 9.0;
    direction[1] = 4.0 / 9.0;
    direction[2] = 8.0 / 9.0;

    pert[0] = 0.0;
    pert[1] = step;
    pert[2] = 0.0;
    pert[3] = 0.0;
    pert[4] = 0.0;

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] += pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fp), "euler");

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] -= pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fm), "euler");

    for (i = 0; i < 5; i++) jac[i] = (fp[i] - fm[i]) / (2.0 * step);

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");

    i = 1;
    RWDS(jac[0], dflux_dcons[0 + i * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + i * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + i * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + i * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + i * 5], tol, "energy flux");
  }

  { /* dEuler/dmo-y */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], pert[5], cons_offset[5], prim_offset[5], direction[3];
    REF_DBL jac[5], fp[5], fm[5];
    REF_DBL step = 1.0e-8;
    REF_DBL tol = 1.e-8;
    REF_INT i;
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    direction[0] = 1.0 / 9.0;
    direction[1] = 4.0 / 9.0;
    direction[2] = 8.0 / 9.0;

    pert[0] = 0.0;
    pert[1] = 0.0;
    pert[2] = step;
    pert[3] = 0.0;
    pert[4] = 0.0;

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] += pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fp), "euler");

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] -= pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fm), "euler");

    for (i = 0; i < 5; i++) jac[i] = (fp[i] - fm[i]) / (2.0 * step);

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");

    i = 2;
    RWDS(jac[0], dflux_dcons[0 + i * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + i * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + i * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + i * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + i * 5], tol, "energy flux");
  }

  { /* dEuler/dmo-z */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], pert[5], cons_offset[5], prim_offset[5], direction[3];
    REF_DBL jac[5], fp[5], fm[5];
    REF_DBL step = 1.0e-8;
    REF_DBL tol = 1.e-8;
    REF_INT i;
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    direction[0] = 1.0 / 9.0;
    direction[1] = 4.0 / 9.0;
    direction[2] = 8.0 / 9.0;

    pert[0] = 0.0;
    pert[1] = 0.0;
    pert[2] = 0.0;
    pert[3] = step;
    pert[4] = 0.0;

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] += pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fp), "euler");

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] -= pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fm), "euler");

    for (i = 0; i < 5; i++) jac[i] = (fp[i] - fm[i]) / (2.0 * step);

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");

    i = 3;
    RWDS(jac[0], dflux_dcons[0 + i * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + i * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + i * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + i * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + i * 5], tol, "energy flux");
  }

  { /* dEuler/dE */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], pert[5], cons_offset[5], prim_offset[5], direction[3];
    REF_DBL jac[5], fp[5], fm[5];
    REF_DBL step = 1.0e-8;
    REF_DBL tol = 1.e-8;
    REF_INT i;
    state[0] = 0.8;
    state[1] = 0.2;
    state[2] = 0.03;
    state[3] = 0.3;
    state[4] = 0.9 / 1.4;
    direction[0] = 1.0 / 9.0;
    direction[1] = 4.0 / 9.0;
    direction[2] = 8.0 / 9.0;

    pert[0] = 0.0;
    pert[1] = 0.0;
    pert[2] = 0.0;
    pert[3] = 0.0;
    pert[4] = step;

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] += pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fp), "euler");

    RSS(ref_phys_make_conserved(state, cons_offset), "cons");
    for (i = 0; i < 5; i++) cons_offset[i] -= pert[i];
    RSS(ref_phys_make_primitive(cons_offset, prim_offset), "cons");
    RSS(ref_phys_euler(prim_offset, direction, fm), "euler");

    for (i = 0; i < 5; i++) jac[i] = (fp[i] - fm[i]) / (2.0 * step);

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");

    i = 4;
    RWDS(jac[0], dflux_dcons[0 + i * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + i * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + i * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + i * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + i * 5], tol, "energy flux");
  }

  { /* dEuler/dE */
    REF_DBL dflux_dcons[25];
    REF_DBL state[5], jac[5], direction[3];
    REF_DBL tol = -1;

    state[0] = 0.9000090795737405;
    state[1] = 0.1;
    state[2] = 0.15;
    state[3] = 0.0;
    state[4] = 0.5102040816326531;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;

    RSS(ref_phys_euler_jac(state, direction, dflux_dcons), "jac");
    /*#include "ref_matrix.h"
      ref_matrix_show_ab(5, 5, dflux_dcons);*/

    jac[0] = 0;
    jac[1] = -0.0035;
    jac[2] = -0.015;
    jac[3] = 0;
    jac[4] = -0.199385696763196;

    RWDS(jac[0], dflux_dcons[0 + 0 * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + 0 * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + 0 * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + 0 * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + 0 * 5], tol, "energy flux");

    jac[0] = 1;
    jac[1] = 0.16;
    jac[2] = 0.15;
    jac[3] = 0;
    jac[4] = 1.99635696763196;

    RWDS(jac[0], dflux_dcons[0 + 1 * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + 1 * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + 1 * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + 1 * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + 1 * 5], tol, "energy flux");

    jac[0] = 0;
    jac[1] = -0.06;
    jac[2] = 0.1;
    jac[3] = 0;
    jac[4] = -0.006;

    RWDS(jac[0], dflux_dcons[0 + 2 * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + 2 * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + 2 * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + 2 * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + 2 * 5], tol, "energy flux");

    jac[0] = 0;
    jac[1] = 0.4;
    jac[2] = 0.0;
    jac[3] = 0;
    jac[4] = 0.14;

    RWDS(jac[0], dflux_dcons[0 + 4 * 5], tol, "mass flux");
    RWDS(jac[1], dflux_dcons[1 + 4 * 5], tol, "mox flux");
    RWDS(jac[2], dflux_dcons[2 + 4 * 5], tol, "moy flux");
    RWDS(jac[3], dflux_dcons[3 + 4 * 5], tol, "moz flux");
    RWDS(jac[4], dflux_dcons[4 + 4 * 5], tol, "energy flux");
  }

  { /* Couette laminar flux */
    REF_DBL state[5], gradient[15], direction[3];
    REF_DBL flux[5];
    REF_DBL turb = -1.0;
    REF_DBL mach = 0.1, re = 40.0, temp = 273.11;
    REF_DBL dudy = 1.0, mu = 1.0;
    REF_DBL thermal_conductivity = mu / ((1.4 - 1.0) * 0.72);
    REF_DBL dpdx = 1.0 / 1.4, dtdx = 1.0;
    REF_INT i;
    for (i = 0; i < 15; i++) gradient[i] = 0.0;
    gradient[1 + 3 * 1] = dudy;
    gradient[0 + 3 * 4] = dpdx;

    state[0] = 1.0;
    state[1] = 0.1;
    state[2] = 0.0;
    state[3] = 0.0;
    state[4] = 1.0 / 1.4;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    RSS(ref_phys_viscous(state, gradient, turb, mach, re, temp, direction,
                         flux),
        "laminar");
    RWDS(0.0, flux[0], -1, "mass flux");
    RWDS(0.0, flux[1], -1, "x mo flux");
    RWDS(mach / re * mu * dudy * direction[0], flux[2], -1, "y mo flux");
    RWDS(0.0, flux[3], -1, "z mo flux");
    RWDS(
        mach / re *
            (mu * dudy * direction[0] * state[2] + thermal_conductivity * dtdx),
        flux[4], -1, "energy flux");
  }

  { /* bulk visc laminar flux */
    REF_DBL state[5], gradient[15], direction[3];
    REF_DBL flux[5];
    REF_DBL turb = -1.0;
    REF_DBL mach = 0.1, re = 40.0, temp = 273.11;
    REF_DBL dvdy = 1.0, mu = 1.0;
    REF_DBL thermal_conductivity = mu / ((1.4 - 1.0) * 0.72), dtdy = 1.0,
            dpdy = 1.0 / 1.4;
    REF_INT i;
    for (i = 0; i < 15; i++) gradient[i] = 0.0;
    gradient[1 + 3 * 2] = dvdy;
    gradient[1 + 3 * 4] = dpdy;

    state[0] = 1.0;
    state[1] = 0.1;
    state[2] = 0.0;
    state[3] = 0.0;
    state[4] = 1.0 / 1.4;
    direction[0] = 0.0;
    direction[1] = 1.0;
    direction[2] = 0.0;
    RSS(ref_phys_viscous(state, gradient, turb, mach, re, temp, direction,
                         flux),
        "laminar");
    RWDS(0.0, flux[0], -1, "mass flux");
    RWDS(0.0, flux[1], -1, "x mo flux");
    RWDS(mach / re * mu * (4.0 / 3.0) * dvdy * direction[1], flux[2], -1,
         "y mo flux");
    RWDS(0.0, flux[3], -1, "z mo flux");
    RWDS(mach / re *
             (mu * (4.0 / 3.0) * dvdy * state[2] + thermal_conductivity * dtdy),
         flux[4], -1, "energy flux");
  }
  { /* convert Spalart-Allmaras turbulence working variable to eddy viscosity
     */
    REF_DBL turb, rho, nu, mut_sa;
    REF_DBL chi, fv1;
    REF_DBL cv1 = 7.1;
    turb = -1.0;
    rho = 1.0;
    nu = 1.0;
    RSS(ref_phys_mut_sa(turb, rho, nu, &mut_sa), "eddy viscosity from SA turb");
    RWDS(0.0, mut_sa, -1, "negative SA turb");

    turb = 1.0;
    rho = 1.0;
    nu = 1.0;
    chi = turb / nu;
    fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));
    RSS(ref_phys_mut_sa(turb, rho, nu, &mut_sa), "eddy viscosity from SA turb");
    RWDS(rho * turb * fv1, mut_sa, -1, "negative SA turb");

    turb = 0.0;
    rho = 1.0;
    nu = 1.0;
    RSS(ref_phys_mut_sa(turb, rho, nu, &mut_sa), "eddy viscosity from SA turb");
    RWDS(0, mut_sa, -1, "negative SA turb");

    turb = 1.e100;
    rho = 1.0;
    nu = 1.e-100;
    REIS(REF_DIV_ZERO, ref_phys_mut_sa(turb, rho, nu, &mut_sa),
         "eddy viscosity from SA turb");
  }

  if (ref_mpi_once(ref_mpi)) {
    char file[] = "ref_phys_test.mapbc";
    FILE *f;
    REF_DICT ref_dict;
    REF_INT id, type;

    f = fopen(file, "w");
    fprintf(f, "3\n");
    fprintf(f, "1 1000\n");
    fprintf(f, "2    2000\n");
    fprintf(f, "3  3000 nobody description stuff\n");
    fprintf(f, "4 4000 ignored\n");
    fclose(f);

    RSS(ref_dict_create(&ref_dict), "create");

    RSS(ref_phys_read_mapbc(ref_dict, file), "read mapbc");

    REIS(3, ref_dict_n(ref_dict), "lines");
    id = 1;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(1000, type, "type");
    id = 2;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(2000, type, "type");
    id = 3;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(3000, type, "type");

    RSS(ref_dict_free(ref_dict), "free");
    REIS(0, remove(file), "test clean up");
  }

  if (ref_mpi_once(ref_mpi)) {
    char file[] = "ref_phys_test.mapbc";
    FILE *f;
    REF_DICT ref_dict;
    REF_INT id, type;

    f = fopen(file, "w");
    fprintf(f, "2\n");
    fprintf(f, "1 4000 solid_wall_top\n");
    fprintf(f, "2 4000 solid_wall_bottom\n");
    fclose(f);

    RSS(ref_dict_create(&ref_dict), "create");

    RSS(ref_phys_read_mapbc(ref_dict, file), "read mapbc");

    REIS(2, ref_dict_n(ref_dict), "lines");
    id = 1;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(4000, type, "type");
    id = 2;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(4000, type, "type");

    RSS(ref_dict_free(ref_dict), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* parse string of 1 bc tag into ref_dict */
    char tags[] = "7";
    REF_INT id, type;
    REF_DICT ref_dict;
    RSS(ref_dict_create(&ref_dict), "create");

    RSS(ref_phys_parse_tags(ref_dict, tags), "read mapbc");

    REIS(1, ref_dict_n(ref_dict), "number of tags");
    id = 7;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(4000, type, "type");

    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* parse string of 2 bc tag into ref_dict */
    char tags[] = "7,15";
    REF_INT id, type;
    REF_DICT ref_dict;
    RSS(ref_dict_create(&ref_dict), "create");

    RSS(ref_phys_parse_tags(ref_dict, tags), "read mapbc");

    REIS(2, ref_dict_n(ref_dict), "number of tags");
    id = 7;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(4000, type, "type");
    id = 15;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(4000, type, "type");

    RSS(ref_dict_free(ref_dict), "free");
  }

  { /* brick zeroth */
    REF_GRID ref_grid;
    FILE *f;
    REF_INT i, node, ldim;
    REF_DBL *field;
    REF_DICT ref_dict;
    REF_BOOL *replace;

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, "ref_phys_test.meshb"), "export");
      f = fopen("ref_phys_test.mapbc", "w");
      fprintf(f, "6\n");
      fprintf(f, "1 5000\n");
      fprintf(f, "2 5000\n");
      fprintf(f, "3 5000\n");
      fprintf(f, "4 5000\n");
      fprintf(f, "5 4000\n"); /* z = 0 */
      fprintf(f, "6 5000\n");
      fclose(f);
      ref_malloc(field, 10 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
      each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
        if (ref_node_xyz(ref_grid_node(ref_grid), 2, node) < 0.01) {
          for (i = 0; i < 10; i++) field[i + 10 * node] = -1.0;
        } else {
          for (i = 0; i < 10; i++) field[i + 10 * node] = (REF_DBL)i;
        }
      }
      RSS(ref_gather_scalar_by_extension(ref_grid, 10, field, NULL,
                                         "ref_phys_test.solb"),
          "gather");
      ref_free(field);
      RSS(ref_grid_free(ref_grid), "free");
    }

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, "ref_phys_test.meshb"),
        "import");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field,
                        "ref_phys_test.solb"),
        "part field");

    RSS(ref_dict_create(&ref_dict), "create");
    RSS(ref_phys_read_mapbc(ref_dict, "ref_phys_test.mapbc"),
        "unable to mapbc");
    ref_malloc(replace, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_BOOL);
    RSS(ref_phys_mask_strong_bcs(ref_grid, ref_dict, replace, ldim), "mask");
    RSS(ref_dict_free(ref_dict), "free");
    RSS(ref_recon_extrapolate_zeroth(ref_grid, field, replace, ldim),
        "extrapolate zeroth order");
    ref_free(replace);

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 6; i < 10; i++) {
        RWDS((REF_DBL)i, field[i + 10 * node], -1.0, "not replaced");
      }
    }
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove("ref_phys_test.meshb"), "meshb clean up");
      REIS(0, remove("ref_phys_test.mapbc"), "mapbc clean up");
      REIS(0, remove("ref_phys_test.solb"), "solb clean up");
    }

    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");
  }

  {
    REF_DBL uplus, yplus;
    uplus = -0.00001;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(-0.00001, yplus, -1, "yplus");

    uplus = 0;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(0, yplus, -1, "yplus");

    uplus = 0.00001;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(uplus, yplus, -1, "yplus");

    uplus = 100;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(0.1108 * exp(0.4 * uplus), yplus, -1, "yplus");

    uplus = 1.0e4;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(0.1108 * exp(0.4 * uplus), yplus, -1, "yplus");

    uplus = 1.0e8;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(0.1108 * exp(0.4 * uplus), yplus, -1, "yplus");

    uplus = 1.0e20;
    RSS(ref_phys_spalding_yplus(uplus, &yplus), "yplus");
    RWDS(0.1108 * exp(0.4 * uplus), yplus, -1, "yplus");
  }

  {
    REF_DBL uplus, dyplus_duplus, yplus0, yplus1, duplus;
    uplus = 0;
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS(1.0, dyplus_duplus, -1, "yplus");
    uplus = 1.0e-11;
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS(1.0, dyplus_duplus, -1, "yplus");
    uplus = 100.0;
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS(0.1108 * 0.4 * exp(0.4 * uplus), dyplus_duplus, -1, "yplus");

    /* finite-difference */
    uplus = 5;
    duplus = 1.0e-6;
    RSS(ref_phys_spalding_yplus(uplus - duplus, &yplus0), "yplus");
    RSS(ref_phys_spalding_yplus(uplus + duplus, &yplus1), "yplus");
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS((yplus1 - yplus0) / (2.0 * duplus), dyplus_duplus, duplus, "yplus");

    /* finite-difference */
    uplus = 11;
    duplus = 1.0e-6;
    RSS(ref_phys_spalding_yplus(uplus - duplus, &yplus0), "yplus");
    RSS(ref_phys_spalding_yplus(uplus + duplus, &yplus1), "yplus");
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS((yplus1 - yplus0) / (2.0 * duplus), dyplus_duplus, duplus, "yplus");

    /* finite-difference */
    uplus = 50;
    duplus = 1.0e-4;
    RSS(ref_phys_spalding_yplus(uplus - duplus, &yplus0), "yplus");
    RSS(ref_phys_spalding_yplus(uplus + duplus, &yplus1), "yplus");
    RSS(ref_phys_spalding_dyplus_duplus(uplus, &dyplus_duplus), "yplus");
    RWDS((yplus1 - yplus0) / (2.0 * duplus), dyplus_duplus, 0.01, "yplus");
  }

  { /* eval spalding forward and back */
    REF_DBL yplus, uplus, y;

    yplus = -6.0e-6;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = -0.00001;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RWDS(-0.00001, uplus, -1, "uplus");

    yplus = 0;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RWDS(0, uplus, -1, "uplus");

    yplus = 1.0e-10;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 1;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 10;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 15;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 30;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 100;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 1.0e4;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 1.0e10;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");

    yplus = 1.0e20;
    RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
    RSS(ref_phys_spalding_yplus(uplus, &y), "y");
    RWDS(y, yplus, -1, "uplus");
  }

  { /* sa surrogate */
    REF_DBL wall_distance, nu_tilde;
    wall_distance = -1.0;
    RSS(ref_phys_sa_surrogate(wall_distance, &nu_tilde), "sa soln");
    RWDS(0, nu_tilde, -1, "nu_tilde neg");
    wall_distance = 0.0;
    RSS(ref_phys_sa_surrogate(wall_distance, &nu_tilde), "sa soln");
    RWDS(0, nu_tilde, -1, "nu_tilde 0");
    wall_distance = 0.05;
    RSS(ref_phys_sa_surrogate(wall_distance, &nu_tilde), "sa soln");
    RWDS(1000.0, nu_tilde, -1, "nu_tilde half BL");
    wall_distance = 1e6;
    RSS(ref_phys_sa_surrogate(wall_distance, &nu_tilde), "sa soln");
    RWDS(0, nu_tilde, -1, "nu_tilde far");
  }

  { /* mid tri signed dist */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *field, *distance;
    REF_DBL offset = 0.5;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "tri");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(field, ref_node_max(ref_node), REF_DBL);
    ref_malloc(distance, ref_node_max(ref_node), REF_DBL);
    each_ref_node_valid_node(ref_node, node) {
      field[node] = ref_node_xyz(ref_node, 0, node) - offset;
    }
    RSS(ref_phys_signed_distance(ref_grid, field, distance), "dist");
    each_ref_node_valid_node(ref_node, node) {
      RWDS(ref_node_xyz(ref_node, 0, node) - offset, distance[node], -1,
           "dist");
    }
    ref_free(distance);
    ref_free(field);
    ref_grid_free(ref_grid);
  }

  { /* shifted tri signed dist */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *field, *distance;
    REF_DBL offset = 0.2;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "tri");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(field, ref_node_max(ref_node), REF_DBL);
    ref_malloc(distance, ref_node_max(ref_node), REF_DBL);
    each_ref_node_valid_node(ref_node, node) {
      field[node] = ref_node_xyz(ref_node, 0, node) - offset;
    }
    RSS(ref_phys_signed_distance(ref_grid, field, distance), "dist");
    each_ref_node_valid_node(ref_node, node) {
      RWDS(ref_node_xyz(ref_node, 0, node) - offset, distance[node], -1,
           "dist");
    }
    ref_free(distance);
    ref_free(field);
    ref_grid_free(ref_grid);
  }

  { /* two tri signed dist */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_DBL *field, *distance;
    REF_DBL offset = 0.5;

    RSS(ref_fixture_tri2_grid(&ref_grid, ref_mpi), "tri");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(field, ref_node_max(ref_node), REF_DBL);
    ref_malloc(distance, ref_node_max(ref_node), REF_DBL);
    each_ref_node_valid_node(ref_node, node) {
      field[node] = ref_node_xyz(ref_node, 0, node) - offset;
    }
    RSS(ref_phys_signed_distance(ref_grid, field, distance), "dist");
    each_ref_node_valid_node(ref_node, node) {
      RWDS(ref_node_xyz(ref_node, 0, node) - offset, distance[node], -1,
           "dist");
    }
    ref_free(distance);
    ref_free(field);
    ref_grid_free(ref_grid);
  }

  { /* tri wall dist */
    REF_GRID ref_grid;
    REF_DICT ref_dict;
    REF_DBL *distance;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "tri");
    RSS(ref_dict_create(&ref_dict), "dict");
    ref_malloc_init(distance, ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    -1.0);
    RSS(ref_dict_store(ref_dict, 10, 4000), "store");
    RSS(ref_phys_wall_distance(ref_grid, ref_dict, distance), "store");
    if (!ref_mpi_para(ref_mpi)) {
      RWDS(0.0, distance[0], -1, "n0");
      RWDS(0.0, distance[1], -1, "n1");
      RWDS(1.0, distance[2], -1, "n2");
    }
    ref_free(distance);
    ref_dict_free(ref_dict);
    ref_grid_free(ref_grid);
  }

  { /* tet wall dist */
    REF_GRID ref_grid;
    REF_DICT ref_dict;
    REF_DBL *distance;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "tri");
    RSS(ref_dict_create(&ref_dict), "dict");
    ref_malloc_init(distance, ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    -1.0);
    RSS(ref_dict_store(ref_dict, 10, 4000), "store");
    RSS(ref_phys_wall_distance(ref_grid, ref_dict, distance), "store");
    if (!ref_mpi_para(ref_mpi)) {
      RWDS(0.0, distance[0], -1, "n0");
      RWDS(0.0, distance[1], -1, "n1");
      RWDS(0.0, distance[2], -1, "n2");
      RWDS(1.0, distance[3], -1, "n3");
    }
    ref_free(distance);
    ref_dict_free(ref_dict);
    ref_grid_free(ref_grid);
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
