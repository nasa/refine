
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

#include "ref_phys.h"

#include "ref_malloc.h"
#include "ref_math.h"

#include "ref_recon.h"

REF_STATUS ref_phys_euler(REF_DBL *state, REF_DBL *direction, REF_DBL *flux) {
  REF_DBL rho, u, v, w, p, e, speed;
  REF_DBL gamma = 1.4;

  rho = state[0];
  u = state[1];
  v = state[2];
  w = state[3];
  p = state[4];

  e = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);

  speed = u * direction[0] + v * direction[1] + w * direction[2];

  flux[0] = rho * speed;
  flux[1] = rho * speed * u + p * direction[0];
  flux[2] = rho * speed * v + p * direction[1];
  flux[3] = rho * speed * w + p * direction[2];
  flux[4] = speed * (e + p);

  return REF_SUCCESS;
}
/* I do like CFD, vol 2, page 77, (3.6.8) */
REF_STATUS ref_phys_euler_jac(REF_DBL *state, REF_DBL *direction,
                              REF_DBL *dflux_dcons) {
  REF_DBL rho, u, v, w, p, q2, e, qn;
  REF_DBL gamma = 1.4;
  REF_DBL K, H;
  REF_DBL nx, ny, nz;
  nx = direction[0];
  ny = direction[1];
  nz = direction[2];

  K = gamma - 1;
  rho = state[0];
  u = state[1];
  v = state[2];
  w = state[3];
  p = state[4];

  q2 = (u * u + v * v + w * w);
  e = p / K + 0.5 * rho * q2;
  H = e + p;

  qn = u * nx + v * ny + w * nz;

  dflux_dcons[0 + 0 * 5] = 0.0;
  dflux_dcons[1 + 0 * 5] = 0.5 * K * q2 * nx - u * qn;
  dflux_dcons[2 + 0 * 5] = 0.5 * K * q2 * ny - v * qn;
  dflux_dcons[3 + 0 * 5] = 0.5 * K * q2 * nz - w * qn;
  dflux_dcons[4 + 0 * 5] = (0.5 * K * q2 - H) * qn;

  dflux_dcons[0 + 1 * 5] = nx;
  dflux_dcons[1 + 1 * 5] = u * nx - K * u * nx + qn;
  dflux_dcons[2 + 1 * 5] = v * nx - K * u * ny;
  dflux_dcons[3 + 1 * 5] = w * nx - K * u * nz;
  dflux_dcons[4 + 1 * 5] = H * nx - K * u * qn;

  dflux_dcons[0 + 2 * 5] = nx;
  dflux_dcons[1 + 2 * 5] = u * ny - K * v * nx;
  dflux_dcons[2 + 2 * 5] = v * ny - K * v * ny + qn;
  dflux_dcons[3 + 2 * 5] = w * ny - K * v * nz;
  dflux_dcons[4 + 2 * 5] = H * ny - K * v * qn;

  dflux_dcons[0 + 3 * 5] = nz;
  dflux_dcons[1 + 3 * 5] = u * nz - K * w * nx;
  dflux_dcons[2 + 3 * 5] = v * nz - K * w * ny;
  dflux_dcons[3 + 3 * 5] = w * nz - K * w * nz + qn;
  dflux_dcons[4 + 3 * 5] = H * nz - K * w * qn;

  dflux_dcons[0 + 4 * 5] = 0;
  dflux_dcons[1 + 4 * 5] = K * nx;
  dflux_dcons[2 + 4 * 5] = K * ny;
  dflux_dcons[3 + 4 * 5] = K * nz;
  dflux_dcons[4 + 4 * 5] = gamma * qn;

  return REF_SUCCESS;
}

REF_STATUS ref_phys_viscous(REF_DBL *state, REF_DBL *grad, REF_DBL turb,
                            REF_DBL mach, REF_DBL re, REF_DBL reference_temp,
                            REF_DBL *dir, REF_DBL *flux) {
  REF_DBL rho, u, v, w, p, mu, mu_t, t;
  REF_DBL gamma = 1.4;
  REF_DBL sutherland_constant = 110.56;
  REF_DBL sutherland_temp;
  REF_DBL pr = 0.72;
  REF_DBL turbulent_pr = 0.90;
  REF_DBL tau[3][3], qdot[3];
  REF_INT i, j, k;
  REF_DBL thermal_conductivity, dtdx;

  rho = state[0];
  u = state[1];
  v = state[2];
  w = state[3];
  p = state[4];
  t = gamma * p / rho;

  sutherland_temp = sutherland_constant / reference_temp;
  mu = (1.0 + sutherland_temp) / (t + sutherland_temp) * t * sqrt(t);

  RSS(ref_phys_mut_sa(turb, rho, mu / rho, &mu_t), "eddy viscosity");
  thermal_conductivity =
      -(mu / (pr * (gamma - 1.0)) + mu_t / (turbulent_pr * (gamma - 1.0)));
  mu += mu_t;

  mu *= mach / re;
  thermal_conductivity *= mach / re;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tau[i][j] = mu * (grad[j + 3 * (1 + i)] + grad[i + 3 * (1 + j)]);
    }
    for (k = 0; k < 3; k++) {
      tau[i][i] += (-2.0 / 3.0) * mu * grad[k + 3 * (1 + k)];
    }
  }

  for (i = 0; i < 3; i++) {
    /* t = gamma * p / rho quotient rule */
    dtdx = gamma * (grad[i + 3 * 4] * rho - p * grad[i + 3 * 0]) / rho / rho;
    qdot[i] = thermal_conductivity * dtdx;
  }

  flux[0] = 0.0;
  flux[1] = dir[0] * tau[0][0] + dir[1] * tau[0][1] + dir[2] * tau[0][2];
  flux[2] = dir[0] * tau[1][0] + dir[1] * tau[1][1] + dir[2] * tau[1][2];
  flux[3] = dir[0] * tau[2][0] + dir[1] * tau[2][1] + dir[2] * tau[2][2];
  flux[4] = dir[0] * (u * tau[0][0] + v * tau[0][1] + w * tau[0][2] - qdot[0]) +
            dir[1] * (u * tau[1][0] + v * tau[1][1] + w * tau[1][2] - qdot[1]) +
            dir[2] * (u * tau[2][0] + v * tau[2][1] + w * tau[2][2] - qdot[2]);

  return REF_SUCCESS;
}

REF_STATUS ref_phys_mut_sa(REF_DBL turb, REF_DBL rho, REF_DBL nu,
                           REF_DBL *mut_sa) {
  if (turb > 0.0 && rho > 0.0 && nu > 0.0) {
    REF_DBL chi, chi3, fv1;
    REF_DBL cv1 = 7.1;
    if (!ref_math_divisible(turb, nu)) return REF_DIV_ZERO;
    chi = turb / nu;
    chi3 = pow(chi, 3);
    fv1 = chi3 / (chi3 + pow(cv1, 3));
    *mut_sa = rho * turb * fv1;
  } else {
    *mut_sa = 0.0;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_phys_convdiff(REF_DBL *state, REF_DBL *grad, REF_DBL diffusivity,
                             REF_DBL *dir, REF_DBL *flux) {
  REF_DBL velocity[3];

  velocity[0] = 1.0;
  velocity[1] = 1.0;
  velocity[2] = 1.0;

  flux[0] = ref_math_dot(dir, velocity) * state[0];

  flux[0] -= diffusivity * ref_math_dot(dir, grad);

  return REF_SUCCESS;
}

REF_STATUS ref_phys_read_mapbc(REF_DICT ref_dict, const char *mapbc_filename) {
  FILE *file;
  REF_INT i, n, id, type;
  char buffer[1024];
  file = fopen(mapbc_filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", mapbc_filename);
  RNS(file, "unable to open file");
  RES(1, fscanf(file, "%d", &n), "number of lines");
  for (i = 0; i < n; i++) {
    RES(1, fscanf(file, "%d", &id), "read id");
    RES(1, fscanf(file, "%d", &type), "read type");
    fgets(buffer, sizeof(buffer), file);
    RSS(ref_dict_store(ref_dict, id, type), "store");
  }
  fclose(file);
  return REF_SUCCESS;
}

/* continuous galerkin spike
   for (equ = 0; equ < nequ; equ++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_node_tet_vol(ref_node, nodes, &cell_vol), "vol");
      each_ref_cell_cell_node(ref_cell, cell_node) {
        for (i = 0; i < 3; i++)
          tri_nodes[i] = ref_cell_f2n(ref_cell, i, cell_node, cell);
        RSS(ref_node_tri_normal(ref_node, tri_nodes, normal), "vol");
        for (dir = 0; dir < 3; dir++) {
          normal[dir] /= cell_vol;
        }
        for (i = 0; i < 4; i++) {
          for (dir = 0; dir < 3; dir++) {
            system[equ + nsystem * nodes[cell_node]] +=
                0.25 *
                dual_flux[equ + dir * nequ + nequ + ldim * nodes[i]] *
                normal[dir] * cell_vol;
          }
        }
      }
    }
  }
*/

REF_STATUS ref_phys_cc_fv_res(REF_GRID ref_grid, REF_INT nequ, REF_DBL *flux,
                              REF_DBL *res) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT equ, dir, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL cell_vol, flux_grad[3];
  REF_DBL tet_flux[4], *xyzs[4];

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    each_ref_cell_cell_node(ref_cell, cell_node) {
      xyzs[cell_node] = ref_node_xyz_ptr(ref_node, nodes[cell_node]);
    }
    RSS(ref_node_xyz_vol(xyzs, &cell_vol), "vol");
    for (dir = 0; dir < 3; dir++) {
      for (equ = 0; equ < nequ; equ++) {
        each_ref_cell_cell_node(ref_cell, cell_node) {
          tet_flux[cell_node] =
              flux[equ + dir * nequ + 3 * nequ * nodes[cell_node]];
        }
        RSS(ref_node_xyz_grad(xyzs, tet_flux, flux_grad), "grad");
        each_ref_cell_cell_node(ref_cell, cell_node) {
          res[equ + nequ * nodes[cell_node]] +=
              0.25 * flux_grad[dir] * cell_vol;
        }
      }
    }
  }

  RSS(ref_node_ghost_dbl(ref_node, res, nequ), "ghost res");

  return REF_SUCCESS;
}

REF_STATUS ref_phys_cc_fv_embed(REF_GRID ref_grid, REF_INT nequ, REF_DBL *flux,
                                REF_DBL *res) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT i, equ, dir, cell, cell_node, cell_edge, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL cell_vol, flux_grad[3];
  REF_DBL tet_flux[4], *xyzs[4];
  REF_DBL macro_flux[10], macro_xyz[10][3];
  REF_INT m2n[8][4] = {{0, 4, 5, 6}, {1, 8, 7, 4}, {2, 7, 9, 6}, {3, 6, 9, 8},
                       {4, 6, 9, 5}, {7, 8, 9, 4}, {7, 9, 5, 4}, {8, 6, 9, 4}};
  REF_INT n0, n1, macro;
  REF_RECON_RECONSTRUCTION recon = REF_RECON_L2PROJECTION;
  REF_INT node;
  REF_DBL sdot0, sdot1, *direqu, *fluxgrad;
  REF_BOOL high_order = REF_TRUE;
  /* macro node [edge]
                                  3------9[5]--------2
                                 / \              . /
                                /   \          .   /
                               /     \      .     /
                              /       \  .       /
                             /        .\        /
                          6[2]   5[1]  8[4]   7[3]
                           /    .        \    /
                          /  .            \  /
                         /.                \/
                        0------4[0]--------1
  */

  ref_malloc(direqu, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(fluxgrad, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  for (dir = 0; dir < 3; dir++) {
    for (equ = 0; equ < nequ; equ++) {
      each_ref_node_valid_node(ref_node, node) {
        direqu[node] = flux[equ + dir * nequ + 3 * nequ * node];
      }
      RSS(ref_recon_gradient(ref_grid, direqu, fluxgrad, recon), "grad");
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        each_ref_cell_cell_node(ref_cell, cell_node) {
          for (i = 0; i < 3; i++) {
            macro_xyz[cell_node][i] =
                ref_node_xyz(ref_node, i, nodes[cell_node]);
          }
        }
        each_ref_cell_cell_edge(ref_cell, cell_edge) {
          n0 = nodes[ref_cell_e2n_gen(ref_cell, 0, cell_edge)];
          n1 = nodes[ref_cell_e2n_gen(ref_cell, 1, cell_edge)];
          for (i = 0; i < 3; i++) {
            macro_xyz[4 + cell_edge][i] = 0.5 * (ref_node_xyz(ref_node, i, n0) +
                                                 ref_node_xyz(ref_node, i, n1));
          }
        }

        each_ref_cell_cell_node(ref_cell, cell_node) {
          i = equ + dir * nequ + 3 * nequ * nodes[cell_node];
          macro_flux[cell_node] = flux[i];
        }
        each_ref_cell_cell_edge(ref_cell, cell_edge) {
          n0 = nodes[ref_cell_e2n_gen(ref_cell, 0, cell_edge)];
          n1 = nodes[ref_cell_e2n_gen(ref_cell, 1, cell_edge)];
          macro_flux[4 + cell_edge] =
              0.5 * (flux[equ + dir * nequ + 3 * nequ * n0] +
                     flux[equ + dir * nequ + 3 * nequ * n1]);
          sdot0 = 0.0;
          sdot1 = 0.0;
          for (i = 0; i < 3; i++) {
            sdot0 += (ref_node_xyz(ref_node, i, n1) -
                      ref_node_xyz(ref_node, i, n0)) *
                     fluxgrad[i + 3 * n0];
            sdot1 += (ref_node_xyz(ref_node, i, n1) -
                      ref_node_xyz(ref_node, i, n0)) *
                     fluxgrad[i + 3 * n1];
          }
          if (high_order) macro_flux[4 + cell_edge] += 0.125 * (sdot0 - sdot1);
        }
        for (macro = 0; macro < 8; macro++) {
          each_ref_cell_cell_node(ref_cell, cell_node) {
            xyzs[cell_node] = macro_xyz[m2n[macro][cell_node]];
            tet_flux[cell_node] = macro_flux[m2n[macro][cell_node]];
          }
          RSS(ref_node_xyz_vol(xyzs, &cell_vol), "vol");
          RSS(ref_node_xyz_grad(xyzs, tet_flux, flux_grad), "grad");
          each_ref_cell_cell_node(ref_cell, cell_node) {
            res[equ + nequ * nodes[cell_node]] +=
                0.25 * flux_grad[dir] * cell_vol;
          }
        }
      }
    }
  }

  ref_free(fluxgrad);
  ref_free(direqu);

  RSS(ref_node_ghost_dbl(ref_node, res, nequ), "ghost res");

  return REF_SUCCESS;
}
