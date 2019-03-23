
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

#include "ref_phys.h"

#include "ref_args.h"
#include "ref_mpi.h"

#include "ref_grid.h"

#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_malloc.h"
#include "ref_part.h"

#include "ref_recon.h"

static REF_STATUS ref_phys_mask_strong_bcs(REF_GRID ref_grid, REF_DICT ref_dict,
                                           REF_BOOL *replace, REF_INT ldim) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], cell_node;
  REF_INT first, last, i, node, bc;

  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < ldim; i++) {
      replace[i + ldim * node] = REF_FALSE;
    }
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
        "dict bc");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      node = nodes[cell_node];
      if (4000 == bc) {
        first = ldim / 2 + 1; /* first momentum */
        last = ldim / 2 + 4;  /* energy */
        for (i = first; i <= last; i++) replace[i + ldim * node] = REF_TRUE;
        if (12 == ldim) replace[11 + ldim * node] = REF_TRUE; /* turb */
      }
    }
  }

  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_INT laminar_flux_pos = REF_EMPTY;
  REF_INT euler_flux_pos = REF_EMPTY;
  REF_INT mask_pos = REF_EMPTY;
  REF_INT cont_res_pos = REF_EMPTY;

  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--laminar-flux", &laminar_flux_pos),
      REF_NOT_FOUND, "arg search");
  RXS(ref_args_find(argc, argv, "--euler-flux", &euler_flux_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--mask", &mask_pos), REF_NOT_FOUND,
      "arg search");
  RXS(ref_args_find(argc, argv, "--cont-res", &cont_res_pos), REF_NOT_FOUND,
      "arg search");

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
          dual_flux[i + 5 + 5 * dir + 20 * node] += flux[i];
        }
      }
    }
    ref_free(grad);
    ref_free(primitive_dual);

    if (ref_mpi_once(ref_mpi)) printf("writing dual_flux %s\n", argv[7]);
    RSS(ref_gather_scalar(ref_grid, 20, dual_flux, argv[7]),
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
    REF_INT node, i, dir;
    REF_DBL direction[3], state[5], flux[5];

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
    REIS(10, ldim, "expected 10 (rho,u,v,w,p,5*adj) primitive_dual");

    if (ref_mpi_once(ref_mpi)) printf("copy dual\n");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 0; i < 5; i++) {
        dual_flux[i + 20 * node] = primitive_dual[5 + i + 10 * node];
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
        state[i] = primitive_dual[i + 10 * node];
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

    if (ref_mpi_once(ref_mpi)) printf("writing dual_flux %s\n", argv[4]);
    RSS(ref_gather_scalar(ref_grid, 20, dual_flux, argv[4]),
        "export dual_flux");

    ref_free(dual_flux);

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
    RAS(10 == ldim || 12 == ldim,
        "expected 10 (rho,u,v,w,p,5*adj) or 12 (rho,u,v,w,p,turb,6*adj)");

    ref_malloc(replace, ldim * ref_node_max(ref_grid_node(ref_grid)), REF_BOOL);
    RSS(ref_phys_mask_strong_bcs(ref_grid, ref_dict, replace, ldim), "mask");
    RSS(ref_recon_extrapolate_zeroth(ref_grid, primitive_dual, replace, ldim),
        "extrapolate zeroth order");
    ref_free(replace);

    if (ref_mpi_once(ref_mpi))
      printf("writing strong_replacement %s\n", argv[5]);
    RSS(ref_gather_scalar(ref_grid, ldim, primitive_dual, argv[5]),
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
    REF_DBL *dual_flux, *system, *flux, *weight;
    REF_INT ldim;
    REF_INT equ, dir, node, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL cell_vol, flux_grad[3];
    REF_DBL convergence_rate, exponent;
    REF_INT nsystem = 11;

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
    if (convergence_rate > 0.0) exponent = 1.0 / convergence_rate;
    if (ref_mpi_once(ref_mpi))
      printf("convergence rate %f exponent %f\n", convergence_rate, exponent);

    if (ref_mpi_once(ref_mpi)) printf("reading grid %s\n", argv[2]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]),
        "unable to load target grid in position 2");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tet(ref_grid);

    if (ref_mpi_once(ref_mpi)) printf("reading dual flux %s\n", argv[3]);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &dual_flux, argv[3]),
        "unable to load scalar in position 3");
    REIS(20, ldim, "expected 20 (5*adj,5*xflux,5*yflux,5*zflux) scalar");

    ref_malloc_init(weight, ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                    0.0);
    ref_malloc_init(flux, ref_node_max(ref_grid_node(ref_grid)), REF_DBL, 0.0);
    ref_malloc_init(system, nsystem * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 0.0);

    for (dir = 0; dir < 3; dir++) {
      for (equ = 0; equ < 5; equ++) {
        each_ref_node_valid_node(ref_node, node) {
          flux[node] = dual_flux[equ + dir * 5 + 5 + ldim * node];
        }
        each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
          RSS(ref_node_tet_vol(ref_node, nodes, &cell_vol), "vol");
          RSS(ref_node_tet_grad(ref_node, nodes, flux, flux_grad), "grad");
          each_ref_cell_cell_node(ref_cell, cell_node) {
            system[equ + nsystem * nodes[cell_node]] +=
                0.25 * flux_grad[dir] * cell_vol;
          }
        }
      }
    }

    each_ref_node_valid_node(ref_node, node) {
      for (equ = 0; equ < 5; equ++) {
        system[equ + 5 + nsystem * node] = dual_flux[equ + ldim * node];
      }
    }
    each_ref_node_valid_node(ref_node, node) {
      for (equ = 0; equ < 5; equ++) {
        system[nsystem - 1 + nsystem * node] +=
            ABS(dual_flux[equ + ldim * node] * system[equ + 6 * node]);
        weight[node] +=
            ABS(dual_flux[equ + ldim * node] * system[equ + 6 * node]);
      }
      if (weight[node] > 0.0) weight[node] = pow(weight[node], exponent);
    }
    if (ref_mpi_once(ref_mpi)) printf("writing system.tec\n");
    RSS(ref_gather_scalar_by_extension(ref_grid, nsystem, system, NULL,
                                       "system.tec"),
        "export primitive_dual");
    if (ref_mpi_once(ref_mpi)) printf("writing weight %s\n", argv[5]);
    RSS(ref_gather_scalar(ref_grid, 1, weight, argv[5]), "export weight");

    ref_free(weight);
    ref_free(flux);
    ref_free(system);
    ref_free(dual_flux);

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
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

  { /* Couette laminar flux */
    REF_DBL state[5], gradient[15], direction[3];
    REF_DBL flux[5];
    REF_DBL turb = -1.0;
    REF_DBL mach = 0.1, re = 10.0, temp = 273.0;
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
    REF_DBL mach = 0.1, re = 10.0, temp = 273.0;
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

  if (!ref_mpi_para(ref_mpi)) {
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

  if (!ref_mpi_para(ref_mpi)) {
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

  { /* brick zeroth */
    REF_GRID ref_grid;
    FILE *f;
    REF_INT i, node, ldim;
    REF_DBL *field;

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
    REIS(
        0,
        system("./ref_phys_test --mask ref_phys_test.meshb ref_phys_test.mapbc "
               "ref_phys_test.solb ref_phys_test_replace.solb > /dev/null"),
        "mask");

    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &field,
                        "ref_phys_test_replace.solb"),
        "part field");

    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      for (i = 6; i < 10; i++) {
        RWDS((REF_DBL)i, field[i + 10 * node], -1.0, "not repalced");
      }
    }
    if (ref_mpi_once(ref_mpi)) {
      REIS(0, remove("ref_phys_test.meshb"), "meshb clean up");
      REIS(0, remove("ref_phys_test.mapbc"), "mapbc clean up");
      REIS(0, remove("ref_phys_test.solb"), "solb clean up");
      REIS(0, remove("ref_phys_test_replace.solb"), "solb clean up");
    }

    ref_free(field);
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
