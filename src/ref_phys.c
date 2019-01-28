
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

REF_STATUS ref_phys_laminar(REF_DBL *state, REF_DBL *grad, REF_DBL mach,
                            REF_DBL re, REF_DBL reference_temp, REF_DBL *dir,
                            REF_DBL *flux) {
  REF_DBL rho, u, v, w, p, mu, t;
  REF_DBL gamma = 1.4;
  REF_DBL sutherland_constant = 198.6;
  REF_DBL sutherland_temp;
  REF_DBL pr = 0.72;
  REF_DBL tau[3][3], qdot[3];
  REF_INT i, j, k;

  rho = state[0];
  u = state[1];
  v = state[2];
  w = state[3];
  p = state[4];
  t = gamma * p / rho;

  sutherland_temp = sutherland_constant / reference_temp;
  mu = (1.0 + sutherland_temp) / (t + sutherland_temp) * t * sqrt(t);
  mu = mach / re * mu;

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
    qdot[i] = -mu / (pr * (gamma - 1.0)) * gamma *
              (grad[i + 3 * 4] * rho - p * grad[i + 3 * 0]) / rho / rho;
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
