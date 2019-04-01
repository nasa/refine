
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

REF_STATUS ref_phys_viscous(REF_DBL *state, REF_DBL *grad, REF_DBL turb,
                            REF_DBL mach, REF_DBL re, REF_DBL reference_temp,
                            REF_DBL *dir, REF_DBL *flux) {
  REF_DBL rho, u, v, w, p, mu, mu_t, t;
  REF_DBL gamma = 1.4;
  REF_DBL sutherland_constant = 198.6;
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
  mu = mach / re * mu;

  RSS(ref_phys_mut_sa(turb, rho, mu / rho, &mu_t), "eddy viscosity");
  thermal_conductivity =
      -(mu / (pr * (gamma - 1.0)) + mu_t / (turbulent_pr * (gamma - 1.0)));
  mu += mu_t;

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
  REF_INT equ, dir, node, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL cell_vol, flux_grad[3];
  REF_DBL *equ_flux;
  ref_malloc(equ_flux, ref_node_max(ref_node), REF_DBL);

  for (dir = 0; dir < 3; dir++) {
    for (equ = 0; equ < nequ; equ++) {
      each_ref_node_valid_node(ref_node, node) {
        equ_flux[node] = flux[equ + dir * nequ + 3 * nequ * node];
      }
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        RSS(ref_node_tet_vol(ref_node, nodes, &cell_vol), "vol");
        RSS(ref_node_tet_grad_nodes(ref_node, nodes, equ_flux, flux_grad),
            "grad");
        each_ref_cell_cell_node(ref_cell, cell_node) {
          res[equ + nequ * nodes[cell_node]] +=
              0.25 * flux_grad[dir] * cell_vol;
        }
      }
    }
  }
  ref_free(equ_flux);
  return REF_SUCCESS;
}
