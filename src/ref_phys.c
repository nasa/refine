
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
#include <string.h>

#include "ref_egads.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_part.h"
#include "ref_recon.h"
#include "ref_sort.h"

REF_STATUS ref_phys_make_primitive(REF_DBL *conserved, REF_DBL *primitive) {
  REF_DBL rho, u, v, w, p, e;
  REF_DBL gamma = 1.4;
  rho = conserved[0];
  u = conserved[1] / conserved[0];
  v = conserved[2] / conserved[0];
  w = conserved[3] / conserved[0];
  e = conserved[4];
  p = (gamma - 1.0) * (e - 0.5 * rho * (u * u + v * v + w * w));

  primitive[0] = rho;
  primitive[1] = u;
  primitive[2] = v;
  primitive[3] = w;
  primitive[4] = p;

  return REF_SUCCESS;
}
REF_STATUS ref_phys_make_conserved(REF_DBL *primitive, REF_DBL *conserved) {
  REF_DBL rho, u, v, w, p, e;
  REF_DBL gamma = 1.4;
  rho = primitive[0];
  u = primitive[1];
  v = primitive[2];
  w = primitive[3];
  p = primitive[4];
  e = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
  conserved[0] = rho;
  conserved[1] = rho * u;
  conserved[2] = rho * v;
  conserved[3] = rho * w;
  conserved[4] = e;

  return REF_SUCCESS;
}

REF_STATUS ref_phys_entropy_adjoint(REF_DBL *primitive, REF_DBL *dual) {
  REF_DBL rho, u, v, w, p, s, vel2;
  REF_DBL gamma = 1.4;
  rho = primitive[0];
  u = primitive[1];
  v = primitive[2];
  w = primitive[3];
  p = primitive[4];
  /* entropy adjoint, Equ. (11), AIAA 2009-3790 */
  s = log(p / pow(rho, gamma));
  vel2 = u * u + v * v + w * w;
  dual[0] = (gamma - s) / (gamma - 1.0) - 0.5 * rho * vel2 / p;
  dual[1] = rho * u / p;
  dual[2] = rho * v / p;
  dual[3] = rho * w / p;
  dual[4] = -rho / p;

  return REF_SUCCESS;
}

REF_STATUS ref_phys_entropy_flux(REF_DBL *primitive, REF_DBL *flux) {
  REF_DBL rho, u, v, w, p, s, U;
  REF_DBL gamma = 1.4;

  rho = primitive[0];
  u = primitive[1];
  v = primitive[2];
  w = primitive[3];
  p = primitive[4];

  s = log(p / pow(rho, gamma));
  U = -rho * s / (gamma - 1.0);

  flux[0] = U * u;
  flux[1] = U * v;
  flux[2] = U * w;

  return REF_SUCCESS;
}

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

  rho = state[0];
  u = state[1];
  v = state[2];
  w = state[3];
  p = state[4];

  K = gamma - 1;
  q2 = (u * u + v * v + w * w);
  e = p / K + 0.5 * rho * q2;
  H = (e + p) / rho;

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

  dflux_dcons[0 + 2 * 5] = ny;
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

  mu = 1.0;
  if (reference_temp > 0.0) {
    sutherland_temp = sutherland_constant / reference_temp;
    mu = (1.0 + sutherland_temp) / (t + sutherland_temp) * t * sqrt(t);
  }

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

REF_STATUS ref_phys_euler_dual_flux(REF_GRID ref_grid, REF_INT ldim,
                                    REF_DBL *primitive_dual,
                                    REF_DBL *dual_flux) {
  REF_INT i, node, dir, nvar;
  REF_DBL direction[3], state[5], flux[5];

  RAS(0 == ldim % 2, "expect even ldim");
  RAS(ldim >= 10, "expect ldim >= 10");
  nvar = ldim / 2;
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < 5; i++) {
      dual_flux[i + 20 * node] = primitive_dual[nvar + i + ldim * node];
    }
  }

  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < 15; i++) {
      dual_flux[5 + i + 20 * node] = 0.0;
    }
  }

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
  return REF_SUCCESS;
}

REF_STATUS ref_phys_mask_strong_bcs(REF_GRID ref_grid, REF_DICT ref_dict,
                                    REF_BOOL *replace, REF_INT ldim) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], cell_node;
  REF_INT first, last, i, node, bc;
  REF_INT nequ;

  nequ = 0;
  if (0 == ldim % 5) nequ = 5;
  if (0 == ldim % 6) nequ = 6;

  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < ldim; i++) {
      replace[i + ldim * node] = REF_FALSE;
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  if (ref_grid_twod(ref_grid)) ref_cell = ref_grid_edg(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    bc = REF_EMPTY;
    RXS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
        REF_NOT_FOUND, "bc");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      node = nodes[cell_node];
      if (4000 == bc) {
        first = nequ + 1; /* first momentum */
        last = nequ + 4;  /* energy */
        for (i = first; i <= last; i++) replace[i + ldim * node] = REF_TRUE;
        if (6 == nequ) replace[nequ + 5 + ldim * node] = REF_TRUE; /* turb */
      }
    }
  }

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

REF_STATUS ref_phys_read_mapbc_token(REF_DICT ref_dict,
                                     const char *mapbc_filename,
                                     const char *token) {
  FILE *file;
  REF_INT i, n, id, type;
  char name[1024];
  file = fopen(mapbc_filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", mapbc_filename);
  RNS(file, "unable to open file");
  RES(1, fscanf(file, "%d", &n), "number of lines");
  for (i = 0; i < n; i++) {
    RES(1, fscanf(file, "%d", &id), "read id");
    RES(1, fscanf(file, "%d", &type), "read type");
    REIS(32, fgetc(file), "expected space");
    fgets(name, sizeof(name), file);
    /* printf(">%s<>%s<\n",name,token); */
    if (0 == strncmp(name, token, strlen(token))) {
      RSS(ref_dict_store(ref_dict, id, type), "store");
    }
  }
  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_phys_parse_tags(REF_DICT ref_dict, const char *tags) {
  REF_INT id, type;
  char *token, *copy;
  RNS(tags, "tags string is NULL");
  copy = strdup(tags);
  RNS(copy, "tags copy is NULL");
  type = 4000;
  token = strtok(copy, ",");
  while (token != NULL) {
    id = atoi(token);
    RSS(ref_dict_store(ref_dict, id, type), "store");
    token = strtok(NULL, ",");
  }
  free(copy);

  return REF_SUCCESS;
}

REF_STATUS ref_phys_av_tag_attributes(REF_DICT ref_dict, REF_GEOM ref_geom) {
  const char *patch_type;
  REF_INT id, type;
  REF_STATUS ref_status;
  each_ref_geom_face_id(ref_geom, id) {
    ref_status = ref_egads_get_attribute(ref_geom, REF_GEOM_FACE, id,
                                         "av:patch_type", &patch_type);
    if (REF_SUCCESS == ref_status) {
      if (0 == strcmp("noslipwall", patch_type)) {
        type = 4000;
        RSS(ref_dict_store(ref_dict, id, type), "store");
      }
    }
  }
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

REF_STATUS ref_phys_spalding_yplus(REF_DBL uplus, REF_DBL *yplus) {
  *yplus = uplus + 0.1108 * (exp(0.4 * uplus) - 1.0 - 0.4 * uplus);
  return REF_SUCCESS;
}

REF_STATUS ref_phys_spalding_dyplus_duplus(REF_DBL uplus,
                                           REF_DBL *dyplus_duplus) {
  *dyplus_duplus = 1.0 + 0.1108 * (exp(0.4 * uplus) * 0.4 - 0.4);
  return REF_SUCCESS;
}

/*
   y = 0.1108 * exp(0.4 * u)
   y / 0.1108 = exp(0.4 * u)
   log(y / 0.1108) = 0.4 * u
   log(y / 0.1108) / 0.4 = u
*/
REF_STATUS ref_phys_spalding_uplus(REF_DBL yplus, REF_DBL *uplus) {
  REF_DBL u, y, error, dyplus_duplus;
  REF_BOOL keep_going;
  REF_INT iters;
  u = yplus;
  if (u > 12) u = log(yplus / 0.1108) / 0.4;

  iters = 0;
  keep_going = REF_TRUE;
  while (keep_going) {
    RSS(ref_phys_spalding_yplus(u, &y), "yplus");
    error = y - yplus;

    RSS(ref_phys_spalding_dyplus_duplus(u, &dyplus_duplus), "dyplus");

    u = u - error / dyplus_duplus;

    if (ref_math_divisible(error, y) && ABS(y) > 1.0e-3) {
      keep_going = (ABS(error / y) > 1.0e-12);
    } else {
      keep_going = (ABS(error) > 1.0e-15);
    }

    iters++;
    RAB(iters < 100, "iteration count exceeded",
        { printf(" y %e u %e err %e dydu %e\n", y, u, error, dyplus_duplus); });
  }

  *uplus = u;

  return REF_SUCCESS;
}

REF_STATUS ref_phys_minspac(REF_DBL reynolds_number, REF_DBL *yplus1) {
  REF_DBL cf;
  *yplus1 = -1.0;
  /* laminar: 1.0 / ( sqrt(reynolds_number) *20.0 ); */
  RAS(0.06 * reynolds_number > 0, "expected positive reynolds_number");
  RAS(ref_math_divisible(0.445, pow(log(0.06 * reynolds_number), 2)),
      "can not invert for cf");
  cf = 0.455 / pow(log(0.06 * reynolds_number), 2);
  RAS(ref_math_divisible(2.0, cf), "can not invert cf");
  RAS(ref_math_divisible(sqrt(2.0 / cf), reynolds_number),
      "can not invert reynolds_number for yplus1");

  *yplus1 = sqrt(2.0 / cf) / reynolds_number;
  return REF_SUCCESS;
}

REF_STATUS ref_phys_sa_surrogate(REF_DBL wall_distance, REF_DBL *nu_tilde) {
  REF_DBL max_nu_tilde = 1000;
  REF_DBL freestream_nu_tilde = 3;
  REF_DBL boundary_layer_height = 0.1;
  REF_DBL x, height;

  *nu_tilde = 0.0;
  if (wall_distance > 0.0) {
    height = wall_distance / boundary_layer_height;
    x = 4.0 * (height - height * height);
    *nu_tilde = max_nu_tilde * x;
    if (height > 0.5) *nu_tilde = MAX(*nu_tilde, freestream_nu_tilde);
  }
  return REF_SUCCESS;
}

/* use? Fast sweeping methods for eikonal equations on triangular meshes */
REF_STATUS ref_phys_signed_distance(REF_GRID ref_grid, REF_DBL *field,
                                    REF_DBL *distance) {
  REF_DBL node_min, node_max;
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL gradient[3], slope, dist, len;
  REF_INT cell_node, cell, node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *set;
  REF_BOOL update;
  REF_INT passes;
  ref_malloc_init(set, ref_node_max(ref_node), REF_INT, 0);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    node_min = MIN(MIN(field[nodes[0]], field[nodes[1]]), field[nodes[2]]);
    node_max = MAX(MAX(field[nodes[0]], field[nodes[1]]), field[nodes[2]]);
    if (node_min <= 0.0 && 0.0 <= node_max) {
      RSS(ref_node_tri_grad_nodes(ref_node, nodes, field, gradient), "grad");
      each_ref_cell_cell_node(ref_cell, cell_node) {
        slope = sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                     gradient[2] * gradient[2]);
        if (ref_math_divisible(field[nodes[cell_node]], slope)) {
          distance[nodes[cell_node]] = field[nodes[cell_node]] / slope;
        } else {
          distance[nodes[cell_node]] = 0.0;
        }
        set[nodes[cell_node]] = 1;
      }
    }
  }

  update = REF_TRUE;
  passes = 0;
  while (update) {
    update = REF_FALSE;
    passes++;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (2 == set[nodes[0]] + set[nodes[1]] + set[nodes[2]]) {
        each_ref_cell_cell_node(ref_cell, cell_node) {
          if (0 == set[nodes[cell_node]]) {
            dist = 0;
            node = cell_node + 1;
            if (2 < node) node -= 3;
            len = sqrt(pow(ref_node_xyz(ref_node, 0, nodes[node]) -
                               ref_node_xyz(ref_node, 0, nodes[cell_node]),
                           2) +
                       pow(ref_node_xyz(ref_node, 1, nodes[node]) -
                               ref_node_xyz(ref_node, 1, nodes[cell_node]),
                           2) +
                       pow(ref_node_xyz(ref_node, 2, nodes[node]) -
                               ref_node_xyz(ref_node, 2, nodes[cell_node]),
                           2));
            dist = ABS(distance[nodes[node]]) + len;
            node = cell_node + 2;
            if (2 < node) node -= 3;
            len = sqrt(pow(ref_node_xyz(ref_node, 0, nodes[node]) -
                               ref_node_xyz(ref_node, 0, nodes[cell_node]),
                           2) +
                       pow(ref_node_xyz(ref_node, 1, nodes[node]) -
                               ref_node_xyz(ref_node, 1, nodes[cell_node]),
                           2) +
                       pow(ref_node_xyz(ref_node, 2, nodes[node]) -
                               ref_node_xyz(ref_node, 2, nodes[cell_node]),
                           2));

            dist = MIN(dist, ABS(distance[nodes[node]]) + len);
            if (distance[nodes[node]] >= 0.0) {
              distance[nodes[cell_node]] = dist;
            } else {
              distance[nodes[cell_node]] = -dist;
            }
            set[nodes[cell_node]] = 1;
            update = REF_TRUE;
          }
        }
      }
    }
  }

  ref_free(set);
  return REF_SUCCESS;
}

REF_BOOL ref_phys_wall_distance_bc(REF_INT bc) {
  return (bc == 4000 ||  /* viscous_solid */
          bc == -4000 || /* viscous_solid_trs */
          bc == 4075 ||  /* viscous_wall_rough */
          bc == 4100 ||  /* viscous_wall_function */
          bc == 4110 ||  /* viscous_weak_wall */
          bc == -4110 || /* viscous_weak_trs */
          bc == -4100 || /* viscous_wf_trs */
          bc == 6200 ||  /* block_interface */
          bc == 6210);   /* filter */
}

static REF_STATUS ref_phys_local_wall(REF_GRID ref_grid, REF_DICT ref_dict,
                                      REF_INT *node_per_ptr,
                                      REF_INT *local_ncell_ptr,
                                      REF_DBL **local_xyz_ptr) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT local_ncell;
  REF_DBL *local_xyz;
  REF_INT node_per;
  REF_CELL ref_cell;
  REF_INT i, node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT bc;

  if (ref_grid_twod(ref_grid)) {
    ref_cell = ref_grid_edg(ref_grid);
  } else {
    ref_cell = ref_grid_tri(ref_grid);
  }
  node_per = ref_cell_node_per(ref_cell);
  local_ncell = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    bc = REF_EMPTY;
    RXS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
        REF_NOT_FOUND, "bc");
    if (ref_phys_wall_distance_bc(bc)) {
      local_ncell++;
    }
  }
  if (!ref_grid_twod(ref_grid)) { /* adds quads as two tri */
    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      bc = REF_EMPTY;
      RXS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
          REF_NOT_FOUND, "bc");
      if (ref_phys_wall_distance_bc(bc)) {
        local_ncell += 2;
      }
    }
    ref_cell = ref_grid_tri(ref_grid);
  }
  ref_malloc(local_xyz, 3 * node_per * local_ncell, REF_DBL);
  local_ncell = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    bc = REF_EMPTY;
    RXS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
        REF_NOT_FOUND, "bc");
    if (ref_phys_wall_distance_bc(bc)) {
      for (node = 0; node < node_per; node++) {
        for (i = 0; i < 3; i++) {
          local_xyz[i + 3 * node + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[node]);
        }
      }
      local_ncell++;
    }
  }
  if (!ref_grid_twod(ref_grid)) { /* adds quads as two tri */
    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      bc = REF_EMPTY;
      RXS(ref_dict_value(ref_dict, nodes[ref_cell_id_index(ref_cell)], &bc),
          REF_NOT_FOUND, "bc");
      if (ref_phys_wall_distance_bc(bc)) {
        for (i = 0; i < 3; i++) {
          local_xyz[i + 3 * 0 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[0]);
          local_xyz[i + 3 * 1 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[1]);
          local_xyz[i + 3 * 2 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[2]);
        }
        local_ncell += 1;
        for (i = 0; i < 3; i++) {
          local_xyz[i + 3 * 0 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[0]);
          local_xyz[i + 3 * 1 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[2]);
          local_xyz[i + 3 * 2 + 3 * node_per * local_ncell] =
              ref_node_xyz(ref_node, i, nodes[3]);
        }
        local_ncell += 1;
      }
    }
    ref_cell = ref_grid_tri(ref_grid);
  }

  *node_per_ptr = node_per;
  *local_ncell_ptr = local_ncell;
  *local_xyz_ptr = local_xyz;

  return REF_SUCCESS;
}

static REF_STATUS ref_phys_bcast_parts(REF_MPI ref_mpi, REF_INT *part_complete,
                                       REF_INT *part_ncell, REF_INT max_ncell,
                                       REF_INT node_per, REF_DBL *local_xyz,
                                       REF_INT *ncell, REF_DBL **xyz) {
  REF_INT total;
  REF_INT first_part, last_part;
  REF_INT offset, part, i;
  RAS(*part_complete < ref_mpi_n(ref_mpi), "already complete");
  first_part = *part_complete;
  last_part = *part_complete;
  (*ncell) = part_ncell[first_part];

  for (part = first_part + 1; part < ref_mpi_n(ref_mpi); part++) {
    if ((*ncell) + part_ncell[part] > max_ncell) break;
    last_part = part;
    (*ncell) += part_ncell[part];
  }
  (*part_complete) = last_part + 1;

  total = 3 * node_per * (*ncell);
  ref_malloc(*xyz, total, REF_DBL);

  offset = 0;
  for (part = first_part; part <= last_part; part++) {
    total = 3 * node_per * part_ncell[part];
    if (part == ref_mpi_rank(ref_mpi)) {
      if (part_ncell[part] > 0) {
        RSS(ref_mpi_bcast_from_rank(ref_mpi, local_xyz, total, REF_DBL_TYPE,
                                    part),
            "bcast ncell");
        for (i = 0; i < total; i++) (*xyz)[i + offset] = local_xyz[i];
      }
    } else {
      if (part_ncell[part] > 0) {
        RSS(ref_mpi_bcast_from_rank(ref_mpi, &((*xyz)[offset]), total,
                                    REF_DBL_TYPE, part),
            "bcast ncell");
      }
    }
    offset += total;
  }
  total = 3 * node_per * (*ncell);
  REIS(total, offset, "miscount");

  return REF_SUCCESS;
}

REF_STATUS ref_phys_wall_distance_static(REF_GRID ref_grid, REF_DICT ref_dict,
                                         REF_DBL *distance) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT ncell, local_ncell, *part_ncell, part_complete, max_ncell;
  REF_DBL *local_xyz, *xyz;
  REF_INT local_node_per, node_per;
  REF_INT node, cell;
  REF_INT *permutation, i;
  REF_SEARCH ref_search;
  REF_DBL center[3], radius;
  REF_DBL scale = 1.0 + 1.0e-8;
  REF_BOOL timing = REF_FALSE;

  if (timing) ref_mpi_stopwatch_start(ref_mpi);

  if (ref_grid_twod(ref_grid)) {
    node_per = 2;
  } else {
    node_per = 3;
  }

  each_ref_node_valid_node(ref_node, node) { distance[node] = REF_DBL_MAX; }

  RSS(ref_phys_local_wall(ref_grid, ref_dict, &local_node_per, &local_ncell,
                          &local_xyz),
      "local wall");
  REIS(node_per, local_node_per, "node_per miss match");

  ref_malloc(part_ncell, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &local_ncell, part_ncell, REF_INT_TYPE),
      "allgather part ncell");
  if (timing && ref_mpi_once(ref_mpi)) {
    for (i = 0; i < ref_mpi_n(ref_mpi); i++) {
      printf("%d part %d ncell\n", i, part_ncell[i]);
    }
  }
  if (timing) ref_mpi_stopwatch_stop(ref_mpi, "wall dist init");

  max_ncell = 1000000; /* (8 * 3 * node_per) bytes per element */

  part_complete = 0;
  while (part_complete < ref_mpi_n(ref_mpi)) {
    xyz = NULL;
    ncell = 0;
    RSS(ref_phys_bcast_parts(ref_mpi, &part_complete, part_ncell, max_ncell,
                             node_per, local_xyz, &ncell, &xyz),
        "bcast part");

    if (timing) ref_mpi_stopwatch_stop(ref_mpi, "form-scatter");

    RSS(ref_search_create(&ref_search, ncell), "make search");
    ref_malloc(permutation, ncell, REF_INT);
    RSS(ref_sort_shuffle(ncell, permutation), "shuffle");
    for (i = 0; i < ncell; i++) {
      cell = permutation[i];
      RSS(ref_node_bounding_sphere_xyz(&(xyz[3 * node_per * cell]), node_per,
                                       center, &radius),
          "bound");
      RSS(ref_search_insert(ref_search, cell, center, scale * radius), "ins");
    }
    ref_free(permutation);
    if (timing) ref_mpi_stopwatch_stop(ref_mpi, "create-insert");
    if (timing && ref_mpi_once(ref_mpi)) {
      REF_INT depth;
      RSS(ref_search_depth(ref_search, &depth), "depth");
      printf("ncell %d depth %d\n", ncell, depth);
    }
    if (timing) ref_mpi_stopwatch_stop(ref_mpi, "depth");

    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_search_nearest_element(
              ref_search, node_per, xyz,
              ref_node_xyz_ptr(ref_grid_node(ref_grid), node),
              &(distance[node])),
          "candidates");
    }
    RSS(ref_search_free(ref_search), "free");

    ref_free(xyz);
    if (timing) ref_mpi_stopwatch_stop(ref_mpi, "min(dist)");
  }

  ref_free(part_ncell);
  ref_free(local_xyz);
  return REF_SUCCESS;
}

REF_STATUS ref_phys_wall_distance(REF_GRID ref_grid, REF_DICT ref_dict,
                                  REF_DBL *distance) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nowned, nbalance, nstationary, *stationary;
  REF_INT ncell, local_ncell, *part_ncell, part_complete, max_ncell;
  REF_INT nnode, part;
  REF_INT a_total, b_total;
  REF_INT *a_size, *b_size;
  REF_INT *a_next;
  REF_DBL *a_xyz, *b_xyz;
  REF_DBL *a_dist, *b_dist;
  REF_DBL *local_xyz, *xyz;
  REF_INT local_node_per, node_per;
  REF_INT node, cell;
  REF_INT *permutation, i;
  REF_SEARCH ref_search;
  REF_DBL center[3], radius;
  REF_DBL scale = 1.0 + 1.0e-8;

  if (ref_grid_twod(ref_grid)) {
    node_per = 2;
  } else {
    node_per = 3;
  }

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  nowned = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      nowned++;
    }
  }
  /* limit the nodes that are balanced to prevent all2all overload */
  nbalance = MIN(nowned, REF_INT_MAX / ref_mpi_n(ref_mpi));

  /* count the balanced nodes based on destination, count stationary */
  nnode = 0;
  nstationary = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      if (nnode < nbalance) {
        part = ref_part_implicit(nbalance, ref_mpi_n(ref_mpi), nnode);
        if (ref_mpi_rank(ref_mpi) == part) {
          nstationary++;
        } else {
          a_size[part]++;
        }
        nnode++;
      }
    }
  }
  ref_malloc_init(stationary, nstationary, REF_INT, REF_EMPTY);

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) { a_total += a_size[part]; }
  ref_malloc(a_xyz, 3 * a_total, REF_DBL);
  ref_malloc(a_dist, a_total, REF_DBL);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) { b_total += b_size[part]; }
  ref_malloc(b_xyz, 3 * b_total, REF_DBL);
  ref_malloc(b_dist, b_total, REF_DBL);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  /* fill the balanced nodes based on destination, fill stationary */
  nnode = 0;
  nstationary = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      if (nnode < nbalance) {
        part = ref_part_implicit(nbalance, ref_mpi_n(ref_mpi), nnode);
        if (ref_mpi_rank(ref_mpi) == part) {
          stationary[nstationary] = node;
          nstationary++;
        } else {
          for (i = 0; i < 3; i++) {
            a_xyz[i + 3 * a_next[part]] = ref_node_xyz(ref_node, i, node);
          }
          a_next[part]++;
        }
        nnode++;
      }
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_xyz, a_size, b_xyz, b_size, 3, REF_DBL_TYPE),
      "alltoallv xyz");

  for (node = 0; node < b_total; node++) {
    b_dist[node] = REF_DBL_MAX;
  }
  each_ref_node_valid_node(ref_node, node) { distance[node] = REF_DBL_MAX; }

  RSS(ref_phys_local_wall(ref_grid, ref_dict, &local_node_per, &local_ncell,
                          &local_xyz),
      "local wall");
  REIS(node_per, local_node_per, "node_per miss match");

  ref_malloc(part_ncell, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &local_ncell, part_ncell, REF_INT_TYPE),
      "allgather part ncell");

  max_ncell = 1000000; /* (8 * 3 * node_per) bytes per element */

  part_complete = 0;
  while (part_complete < ref_mpi_n(ref_mpi)) {
    xyz = NULL;
    ncell = 0;
    RSS(ref_phys_bcast_parts(ref_mpi, &part_complete, part_ncell, max_ncell,
                             node_per, local_xyz, &ncell, &xyz),
        "bcast part");

    RSS(ref_search_create(&ref_search, ncell), "make search");
    ref_malloc(permutation, ncell, REF_INT);
    /* permutation randomizes tree insertion to improve balance/reduce depth */
    RSS(ref_sort_shuffle(ncell, permutation), "shuffle");
    for (i = 0; i < ncell; i++) {
      cell = permutation[i];
      RSS(ref_node_bounding_sphere_xyz(&(xyz[3 * node_per * cell]), node_per,
                                       center, &radius),
          "bound");
      RSS(ref_search_insert(ref_search, cell, center, scale * radius), "ins");
    }
    ref_free(permutation);

    for (i = 0; i < b_total; i++) {
      RSS(ref_search_nearest_element(ref_search, node_per, xyz, &(b_xyz[3 * i]),
                                     &(b_dist[i])),
          "balanced candidates");
    }
    for (i = 0; i < nstationary; i++) {
      node = stationary[i];
      RSS(ref_search_nearest_element(
              ref_search, node_per, xyz,
              ref_node_xyz_ptr(ref_grid_node(ref_grid), node),
              &(distance[node])),
          "stationary candidates");
    }
    RSS(ref_search_free(ref_search), "free");

    ref_free(xyz);
  }

  ref_free(part_ncell);
  ref_free(local_xyz);

  RSS(ref_mpi_alltoallv(ref_mpi, b_dist, b_size, a_dist, a_size, 1,
                        REF_DBL_TYPE),
      "alltoallv dist");

  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }
  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      if (nnode < nbalance) {
        part = ref_part_implicit(nbalance, ref_mpi_n(ref_mpi), nnode);
        if (ref_mpi_rank(ref_mpi) == part) {
        } else {
          distance[node] = a_dist[a_next[part]];
          a_next[part]++;
        }
        nnode++;
      }
    }
  }

  RSS(ref_node_ghost_dbl(ref_node, distance, 1), "ghost distance");

  ref_free(b_dist);
  ref_free(b_xyz);
  ref_free(a_dist);
  ref_free(a_xyz);
  ref_free(a_next);
  ref_free(stationary);
  ref_free(b_size);
  ref_free(a_size);
  return REF_SUCCESS;
}

REF_STATUS ref_phys_ddes_blend(REF_DBL mach, REF_DBL reynolds_number,
                               REF_DBL sqrt_vel_grad_dot_grad, REF_DBL distance,
                               REF_DBL nu, REF_DBL *fd) {
  REF_DBL rd;
  REF_DBL kappa = 0.41;
  REF_DBL num, denom;
  num = mach * nu / reynolds_number;
  denom = sqrt_vel_grad_dot_grad * kappa * kappa * distance * distance;
  if (ref_math_divisible(num, denom)) {
    rd = num / denom;
    *fd = 1.0 - tanh(pow(8.0 * rd, 3));
  } else {
    *fd = 0.0;
  };
  return REF_SUCCESS;
}
