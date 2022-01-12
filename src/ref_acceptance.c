
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

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_phys.h"
#include "ref_sort.h"

static REF_STATUS flip_twod_yz(REF_NODE ref_node, REF_INT ldim,
                               REF_DBL *field) {
  REF_INT node;
  REF_DBL temp;
  REF_INT nequ;

  nequ = 0;
  if (ldim > 5 && 0 == ldim % 5) nequ = 5;
  if (ldim > 6 && 0 == ldim % 6) nequ = 6;

  each_ref_node_valid_node(ref_node, node) {
    if (ldim >= 5) {
      temp = field[2 + ldim * node];
      field[2 + ldim * node] = field[3 + ldim * node];
      field[3 + ldim * node] = temp;
    }
    if (nequ > 0) {
      temp = field[nequ + 2 + ldim * node];
      field[nequ + 2 + ldim * node] = field[nequ + 3 + ldim * node];
      field[nequ + 3 + ldim * node] = temp;
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_primal_trig(REF_DBL x, REF_DBL y,
                                             REF_DBL *primitive) {
  REF_DBL gamma = 1.4;
  REF_DBL x0, y0;
  REF_DBL r;

  x0 = 0.0;
  y0 = 0.4;
  r = sqrt(pow(x - x0, 2) + pow(y - y0, 2));
  primitive[0] = 0.1 * tanh(50.0 * (r - 0.5)) + 1.0;
  primitive[1] = 0.1 * sin(2.0 * ref_math_pi * x) + 0.2;
  primitive[2] = 0.2 * sin(2.0 * ref_math_pi * y) + 0.3;
  primitive[3] = 0.0;
  primitive[4] = 0.1 * tanh(50.0 * (x - y)) + 1.0 / gamma;
  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_primal_coax(REF_DBL x, REF_DBL y,
                                             REF_DBL *primitive) {
  REF_DBL r2 = 1.384, r1 = 1.0;
  REF_DBL r, t;
  REF_DBL o2 = 0.4, o1 = 0.1;
  REF_DBL a, b;
  REF_DBL u, p;
  a = (o2 * r2 * r2 - o1 * r1 * r1) / (r2 * r2 - r1 * r1);
  b = (o1 - o2) * r2 * r2 * r1 * r1 / (r2 * r2 - r1 * r1);

  r = sqrt(x * x + y * y);
  t = atan2(y, x);
  u = a * r + b / r;
  p = a * a * r * r / 2.0 + 2.0 * a * b * log(r) - b * b / (2 * r * r);
  primitive[0] = 1.0;
  primitive[1] = sin(t) * u;
  primitive[2] = -cos(t) * u;
  primitive[3] = 0;
  primitive[4] = p / 1.4;
  return REF_SUCCESS;
}

static REF_STATUS ringleb_v(REF_DBL x, REF_DBL y, REF_DBL *v) {
  REF_INT i;
  REF_DBL vp, vc, b, rho, l, denom;
  vp = 0.5;
  vc = 0.0;
  for (i = 0; i < 1000; i++) {
    RAS(1.0 - 0.2 * vp * vp > 0.0, "about to nan");
    b = sqrt(1.0 - 0.2 * vp * vp);
    rho = pow(b, 5);
    l = 1.0 / b + 1.0 / 3.0 / pow(b, 3) + 1.0 / 5.0 / pow(b, 5) -
        0.5 * log((1.0 + b) / (1.0 - b));
    denom = 2.0 * rho * sqrt(pow(x - 0.5 * l, 2) + y * y);
    vc = sqrt(1.0 / denom);
    RAS(isfinite(vc), "nan");
    if (ABS(vp - vc) < 1.0e-15) {
      *v = vc;
      return REF_SUCCESS;
    }
    vp = MIN(2.0, vc);
  }
  printf("x %f y %f vp %f vc %f err %e\n", x, y, vp, vc, ABS(vp - vc));

  return REF_FAILURE;
}

static REF_STATUS ref_acceptance_primal_ringleb(REF_DBL x, REF_DBL y,
                                                REF_DBL *primitive) {
  REF_DBL v, b, rho, p, l, psi, t;
  RSS(ringleb_v(x, y, &v), "fixed point v");
  b = sqrt(1.0 - 0.2 * v * v);
  rho = pow(b, 5);
  p = pow(b, 7);
  l = 1.0 / b + 1.0 / 3.0 / pow(b, 3) + 1.0 / 5.0 / pow(b, 5) -
      0.5 * log((1.0 + b) / (1.0 - b));
  psi = sqrt(0.5 / v / v - rho * (x - 0.5 * l));
  if (psi * v < (1.0 - 1.e14)) {
    t = asin(psi * v);
  } else {
    t = 0.5 * ref_math_pi;
  }
  RAS(isfinite(t), "theta is not finite");
  primitive[0] = rho;
  primitive[1] = -v * cos(t);
  primitive[2] = -v * sin(t);
  primitive[3] = 0;
  primitive[4] = p;

  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_primal_fp_sa(REF_DBL y, REF_DBL *primitive) {
  REF_DBL rho = 1.0;
  REF_DBL v = 0.0;
  REF_DBL w = 0.0;
  REF_DBL p = 1.0 / 1.4;
  REF_DBL u;
  REF_DBL y_bl_edge = 0.02;
  REF_DBL yplus_bl_edge = 10000.0;
  REF_DBL uplus_bl_edge = 27.0;
  REF_DBL yplus, uplus;
  REF_DBL mach = 0.2;
  REF_DBL nu_max = 200;
  REF_DBL nu_min = 3;
  REF_DBL nu;
  yplus = y / y_bl_edge * yplus_bl_edge;
  RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
  uplus = MIN(uplus, uplus_bl_edge);
  u = uplus / uplus_bl_edge * mach;

  nu = nu_max * sin(ref_math_pi * MIN(1.0, y / y_bl_edge));
  if (y > 0.5 * y_bl_edge) nu = MAX(nu, nu_min);

  primitive[0] = rho;
  primitive[1] = u;
  primitive[2] = v;
  primitive[3] = w;
  primitive[4] = p;
  primitive[5] = nu;

  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_primal_vortex(REF_DBL x, REF_DBL y,
                                               REF_DBL *primitive) {
  REF_DBL gamma = 1.4;
  REF_DBL rho, pressure, u, v, w, mach;
  REF_DBL ri = 1.0, mi = 2.25;
  REF_DBL rhoi = 1.0, pi = 1.0 / gamma, ai = 1.0;
  REF_DBL base, a, r, t;
  /* real 2D mesh */
  r = sqrt(x * x + y * y);
  t = atan2(y, x);
  base = 1 - pow(ri / r, 2);
  base = 1 + 0.5 * (gamma + 1.0) * mi * mi * base;
  rho = rhoi * pow(base, 1.0 / (gamma - 1.0));
  pressure = pi * pow(base, gamma / (gamma - 1.0));
  a = sqrt(gamma * pressure / rho);
  mach = ai * mi * ri / (a * r);
  /* fun3d 2D convention */
  u = sin(t) * mach * a;
  v = -cos(t) * mach * a;
  w = 0.0;
  primitive[0] = rho;
  primitive[1] = u;
  primitive[2] = v;
  primitive[3] = w;
  primitive[4] = pressure;

  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_u_lisbon(REF_DBL x, REF_DBL y,
                                          REF_DBL *scalar) {
  REF_DBL sigma, nu, u, v;
  sigma = 4.0;
  RAS(ref_math_divisible(y, x), "lisbon not defined for x=0");
  nu = sigma * y / x;
  u = erf(nu);
  v = 1.0 / (sigma * sqrt(ref_math_pi)) * (1.0 - exp(-nu * nu));
  *scalar = sqrt(u * u + v * v);

  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_u(REF_NODE ref_node, const char *function_name,
                                   REF_DBL *scalar) {
  REF_INT node;
  REF_DBL x, y, z;
  each_ref_node_valid_node(ref_node, node) {
    x = ref_node_xyz(ref_node, 0, node);
    y = ref_node_xyz(ref_node, 1, node);
    z = ref_node_xyz(ref_node, 2, node);
    if (strcmp(function_name, "u5") == 0) {
      scalar[node] = 2.0 * pow(x, 2) + 2.0 * pow(y, 2) + 2.0 * pow(z, 2);
    } else if (strcmp(function_name, "supcyl") == 0) {
      REF_DBL r;
      r = sqrt(x * x + y * y) - 0.5;
      scalar[node] = 1.0 - exp(-50.0 * r);
      if ((x - (-1.0)) < 0.5 * y * y) scalar[node] = 2.0;
    } else if (strcmp(function_name, "combo1") == 0) {
      REF_DBL xy, x3;
      xy = (2.0 * x - 1.0) * (2.0 * y - 1.0);
      x3 = pow((2.0 * x - 1.0), 3);
      scalar[node] = tanh(15.0 * xy - x3);
    } else if (strcmp(function_name, "combo2") == 0) {
      REF_DBL xy;
      xy = (2.0 * x - 1.0) * (2.0 * y - 1.0);
      scalar[node] = tanh(50.0 * (xy - 0.5));
    } else if (strcmp(function_name, "sin50xy") == 0) {
      REF_DBL xy;
      xy = (2.0 * x - 1.0) * (2.0 * y - 1.0);
      if (ABS(xy) <= (ref_math_pi / 50.0)) {
        scalar[node] = sin(50.0 * xy);
      } else if (ABS(xy) <= 2.0 * (ref_math_pi / 50.0)) {
        scalar[node] = sin(50.0 * xy);
      } else {
        scalar[node] = 0.01 * sin(50.0 * xy);
      }
    } else if (strcmp(function_name, "sinfun3") == 0) {
      REF_DBL xyz;
      xyz = (x - 0.4) * (y - 0.4) * (z - 0.4); /* sphere2 */
      if (xyz <= (-1.0 * ref_math_pi / 50.0)) {
        scalar[node] = 0.1 * sin(50. * xyz);
      } else if (xyz <= (2.0 * ref_math_pi / 50.0)) {
        scalar[node] = sin(50.0 * xyz);
      } else {
        scalar[node] = 0.1 * sin(50.0 * xyz);
      }
    } else if (strcmp(function_name, "sinatan3") == 0) {
      REF_DBL eps = 0.1;
      REF_DBL xz;
      xz = x * z;
      scalar[node] =
          0.1 * sin(50.0 * xz) + atan(eps / (sin(5.0 * y) - 2.0 * xz));
    } else if (strcmp(function_name, "tanh3") == 0) {
      scalar[node] = tanh(pow(x + 1.3, 20.0) * pow(y - 0.3, 9.0) * z);
    } else if (strcmp(function_name, "bl3") == 0) {
      REF_DBL f, g, h;
      REF_DBL a = 1.0;
      REF_DBL b = 1.0;
      REF_DBL c = 1.0;
      REF_DBL nu = 0.01;
      REF_DBL scale = -1;
      REF_DBL offset = 1;
      f = (1.0 - exp(-a * (1.0 - x) / nu)) / (1.0 - exp(-a / nu));
      g = (1.0 - exp(-b * (1.0 - y) / nu)) / (1.0 - exp(-b / nu));
      h = (1.0 - exp(-c * (1.0 - z) / nu)) / (1.0 - exp(-c / nu));
      scalar[node] = scale * f * g * h + offset;
    } else if (strcmp(function_name, "uplus") == 0) {
      REF_DBL uplus, yplus, scale;
      scale = 1.0e4;
      yplus = scale * y;
      RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
      scalar[node] = uplus + 5.0e-4 * x * x;
    } else if (strcmp(function_name, "circgap") == 0) {
      REF_DBL uplus0, uplus1, uplus, combined;
      REF_INT step;
      REF_DBL yplus0, yplus1, gap, r, radius, center, scale, hop;
      radius = 1.0;
      gap = 0.1;
      center = -gap - radius;
      scale = 1.0e4;
      hop = 25.0;
      yplus0 = scale * ABS(y);
      RSS(ref_phys_spalding_uplus(yplus0, &uplus0), "uplus");
      r = sqrt(x * x + (y - center) * (y - center));
      yplus1 = scale * (r - radius);
      RSS(ref_phys_spalding_uplus(yplus1, &uplus1), "uplus");
      uplus = MIN(uplus0, uplus1);
      step = MAX(0, (REF_INT)uplus - 25);
      combined = uplus + hop * (REF_DBL)step;
      scalar[node] = combined;
    } else if (strcmp(function_name, "sphuplus") == 0) {
      REF_DBL x0 = 1, y0 = 1, z0 = 1;
      REF_DBL uplus;
      REF_DBL yplus, r, radius, scale;
      radius = 0.5;
      scale = 1.0e4;
      r = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));
      yplus = scale * (r - radius);
      RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
      scalar[node] = uplus;
    } else if (strcmp(function_name, "cyluplus") == 0) {
      REF_DBL uplus;
      REF_DBL yplus, r, radius, scale;
      radius = 0.5;
      scale = 1.0e4;
      r = sqrt(x * x + y * y);
      yplus = scale * (r - radius);
      RSS(ref_phys_spalding_uplus(yplus, &uplus), "uplus");
      scalar[node] = uplus;
    } else if (strcmp(function_name, "dist") == 0) {
      REF_DBL r, radius;
      radius = 0.5;
      r = sqrt(x * x + y * y);
      scalar[node] = (r - radius);
    } else if (strcmp(function_name, "parabola") == 0) {
      scalar[node] = tanh(50 * (x - y * y));
    } else if (strcmp(function_name, "one") == 0) {
      scalar[node] = 1.0;
    } else if (strcmp(function_name, "half") == 0) {
      scalar[node] = 1.0;
      if (y < 0.5) scalar[node] = 0.5;
    } else if (strcmp(function_name, "lisbon") == 0) {
      REF_DBL vel;
      RSS(ref_acceptance_u_lisbon(x, y, &vel), "lisbon");
      scalar[node] = vel;
    } else if (strcmp(function_name, "trig") == 0) {
      REF_DBL rho, pressure, u, v, w, mach;
      REF_DBL primitive[5];
      x = ref_node_xyz(ref_node, 0, node);
      y = ref_node_xyz(ref_node, 1, node);
      RSS(ref_acceptance_primal_trig(ref_node_xyz(ref_node, 0, node),
                                     ref_node_xyz(ref_node, 1, node),
                                     primitive),
          "coax");
      rho = primitive[0];
      u = primitive[1];
      v = primitive[2];
      w = primitive[3];
      pressure = primitive[4];
      mach = sqrt(u * u + v * v + w * w) / sqrt(1.4 * pressure / rho);
      scalar[node] = mach;
    } else {
      printf("%s: %d: %s %s\n", __FILE__, __LINE__, "unknown user function",
             function_name);
      return REF_NOT_FOUND;
    }
  }
  return REF_SUCCESS;
}
static REF_STATUS ref_acceptance_g(REF_NODE ref_node, const char *function_name,
                                   REF_DBL *scalar) {
  REF_INT node;
  REF_DBL x, y;
  each_ref_node_valid_node(ref_node, node) {
    x = ref_node_xyz(ref_node, 0, node);
    y = ref_node_xyz(ref_node, 1, node);
    if (strcmp(function_name, "parabola") == 0) {
      /* tanh(50 * (x - y * y)) */
      scalar[0 + 3 * node] = 50.0 / pow(cosh(50 * (x - y * y)), 2);
      scalar[1 + 3 * node] = -100.0 * y / pow(cosh(50 * (x - y * y)), 2);
      scalar[2 + 3 * node] = 0.0;
    } else {
      printf("%s: %d: %s %s\n", __FILE__, __LINE__, "unknown user function",
             function_name);
      return REF_NOT_FOUND;
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_q(REF_NODE ref_node, const char *function_name,
                                   REF_INT *ldim, REF_DBL **scalar) {
  REF_INT node;

  *ldim = 5;
  if (strcmp(function_name, "fp-sa") == 0) {
    *ldim = 6;
  }
  if (strcmp(function_name, "sst") == 0) {
    *ldim = 7;
  }
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  each_ref_node_valid_node(ref_node, node) {
    if (strcmp(function_name, "trig") == 0) {
      REF_INT i;
      REF_DBL primitive[5];
      RSS(ref_acceptance_primal_trig(ref_node_xyz(ref_node, 0, node),
                                     ref_node_xyz(ref_node, 1, node),
                                     primitive),
          "trig");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
    } else if (strcmp(function_name, "vortex") == 0) {
      REF_INT i;
      REF_DBL primitive[5];
      RSS(ref_acceptance_primal_vortex(ref_node_xyz(ref_node, 0, node),
                                       ref_node_xyz(ref_node, 1, node),
                                       primitive),
          "vortex");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
    } else if (strcmp(function_name, "coax") == 0) {
      REF_INT i;
      REF_DBL primitive[5];
      RSS(ref_acceptance_primal_coax(ref_node_xyz(ref_node, 0, node),
                                     ref_node_xyz(ref_node, 1, node),
                                     primitive),
          "coax");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
    } else if (strcmp(function_name, "ringleb") == 0) {
      REF_INT i;
      REF_DBL primitive[5];
      RSS(ref_acceptance_primal_ringleb(ref_node_xyz(ref_node, 0, node),
                                        ref_node_xyz(ref_node, 1, node),
                                        primitive),
          "ringleb");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
    } else if (strcmp(function_name, "fp-sa") == 0) {
      REF_INT i;
      REF_DBL primitive[6];
      RSS(ref_acceptance_primal_fp_sa(ref_node_xyz(ref_node, 1, node),
                                      primitive),
          "ringleb");
      for (i = 0; i < 6; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
    } else if (strcmp(function_name, "sst") == 0) {
      REF_INT i;
      for (i = 0; i < 7; i++) (*scalar)[i + (*ldim) * node] = (REF_DBL)(i + 1);
    } else {
      printf("%s: %d: %s %s\n", __FILE__, __LINE__, "unknown user function",
             function_name);
      return REF_NOT_FOUND;
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_acceptance_pd(REF_NODE ref_node,
                                    const char *function_name, REF_INT *ldim,
                                    REF_DBL **scalar) {
  REF_INT node, i;
  *ldim = 10;
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  each_ref_node_valid_node(ref_node, node) {
    if (strcmp(function_name, "trig") == 0) {
      REF_DBL primitive[5], dual[5];
      RSS(ref_acceptance_primal_trig(ref_node_xyz(ref_node, 0, node),
                                     ref_node_xyz(ref_node, 1, node),
                                     primitive),
          "trig");
      RSS(ref_phys_entropy_adjoint(primitive, dual), "entropy adj");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
      for (i = 0; i < 5; i++) (*scalar)[i + 5 + (*ldim) * node] = dual[i];
    } else if (strcmp(function_name, "vortex") == 0) {
      REF_DBL primitive[5], dual[5];
      RSS(ref_acceptance_primal_vortex(ref_node_xyz(ref_node, 0, node),
                                       ref_node_xyz(ref_node, 1, node),
                                       primitive),
          "vortex");
      RSS(ref_phys_entropy_adjoint(primitive, dual), "entropy adj");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
      for (i = 0; i < 5; i++) (*scalar)[i + 5 + (*ldim) * node] = dual[i];
    } else if (strcmp(function_name, "coax") == 0) {
      REF_DBL primitive[5], dual[5];
      RSS(ref_acceptance_primal_coax(ref_node_xyz(ref_node, 0, node),
                                     ref_node_xyz(ref_node, 1, node),
                                     primitive),
          "coax");
      RSS(ref_phys_entropy_adjoint(primitive, dual), "entropy adj");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
      for (i = 0; i < 5; i++) (*scalar)[i + 5 + (*ldim) * node] = dual[i];
    } else if (strcmp(function_name, "ringleb") == 0) {
      REF_DBL primitive[5], dual[5];
      RSS(ref_acceptance_primal_ringleb(ref_node_xyz(ref_node, 0, node),
                                        ref_node_xyz(ref_node, 1, node),
                                        primitive),
          "ringleb");
      RSS(ref_phys_entropy_adjoint(primitive, dual), "entropy adj");
      for (i = 0; i < 5; i++) (*scalar)[i + (*ldim) * node] = primitive[i];
      for (i = 0; i < 5; i++) (*scalar)[i + 5 + (*ldim) * node] = dual[i];
    } else {
      printf("%s: %d: %s %s\n", __FILE__, __LINE__, "unknown user function",
             function_name);
      return REF_NOT_FOUND;
    }
  }
  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid = NULL;
  REF_NODE ref_node;
  REF_INT masabl_pos;
  REF_INT ugawg_pos;
  REF_INT xyz_pos;
  REF_INT twod_pos;
  REF_INT pos;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "-ugawg", &ugawg_pos), REF_NOT_FOUND, "arg");

  if (REF_EMPTY != ugawg_pos) {
    REF_BOOL metric_recognized = REF_FALSE;
    if (5 != argc) {
      printf("usage:\n");
      printf(
          "  %s -ugawg [linear,polar-1,polar-2] input.grid_format "
          "output.metric\n",
          argv[0]);
      RSS(ref_mpi_free(ref_mpi), "mpi free");
      RSS(ref_mpi_stop(), "stop");
      return (1);
    }
    printf("%s reading\n", argv[3]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");

    printf("%s field type\n", argv[2]);
    if (strcmp(argv[2], "linear") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg linear metric\n");
      RSS(ref_metric_olympic_node(ref_grid_node(ref_grid), 0.001), "lin");
    }
    if (strcmp(argv[2], "polar-1") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg polar-1 metric\n");
      RSS(ref_metric_ugawg_node(ref_grid_node(ref_grid), 1), "p1");
    }
    if (strcmp(argv[2], "polar-2") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg polar-2 metric\n");
      RSS(ref_metric_ugawg_node(ref_grid_node(ref_grid), 2), "p2");
    }
    if (strcmp(argv[2], "ring") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg ring metric\n");
      RSS(ref_metric_ring_node(ref_grid_node(ref_grid)), "ring");
    }
    if (strcmp(argv[2], "side") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg side metric\n");
      RSS(ref_metric_side_node(ref_grid_node(ref_grid)), "side");
    }
    if (strcmp(argv[2], "circle") == 0) {
      metric_recognized = REF_TRUE;
      printf(" -ugawg circle metric\n");
      RSS(ref_metric_circle_node(ref_grid_node(ref_grid)), "side");
    }
    RAS(metric_recognized, "did not recognize metric field name");
    printf("%s metric exported\n", argv[4]);
    RSS(ref_gather_metric(ref_grid, argv[4]), "in");

    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-twod", &twod_pos), REF_NOT_FOUND, "arg");

  if (REF_EMPTY != twod_pos) {
    if (5 != argc || 1 != twod_pos) {
      printf("usage:\n");
      printf("  %s -twod version input.grid_format output.metric\n", argv[0]);
      RSS(ref_mpi_free(ref_mpi), "mpi free");
      RSS(ref_mpi_stop(), "stop");
      return (1);
    }
    printf("%s reading\n", argv[3]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");

    printf("%s version\n", argv[2]);
    RSS(ref_metric_twod_analytic_node(ref_grid_node(ref_grid), argv[2]), "lin");
    printf("%s metric exported\n", argv[4]);
    RSS(ref_gather_metric(ref_grid, argv[4]), "in");
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-masabl", &masabl_pos), REF_NOT_FOUND, "arg");

  if (REF_EMPTY != masabl_pos) {
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "in");

    RSS(ref_metric_masabl_node(ref_grid_node(ref_grid)), "masabl");
    if (ref_grid_twod(ref_grid))
      RSS(ref_metric_twod_node(ref_grid_node(ref_grid)), "2d");

    RSS(ref_gather_metric(ref_grid, argv[2]), "in");
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-xyz", &xyz_pos), REF_NOT_FOUND, "arg");
  if (REF_EMPTY != xyz_pos) {
    REF_DBL *xyz;
    REF_INT node;
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "in");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(xyz, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

    each_ref_node_valid_node(ref_node, node) {
      REF_DBL x = ref_node_xyz(ref_node, 0, node);
      REF_DBL y = ref_node_xyz(ref_node, 1, node);
      REF_DBL z = ref_node_xyz(ref_node, 2, node);
      xyz[0 + 3 * node] = x + 0.10 * sin(2.0 * y * ref_math_pi);
      xyz[1 + 3 * node] = y - 0.10 * sin(2.0 * x * ref_math_pi);
      xyz[2 + 3 * node] = z + 0.10 * sin(z * ref_math_pi);
    }

    RSS(ref_gather_scalar_by_extension(ref_grid, 3, xyz, NULL, argv[3]), "in");

    each_ref_node_valid_node(ref_node, node) {
      ref_node_xyz(ref_node, 0, node) = xyz[0 + 3 * node];
      ref_node_xyz(ref_node, 1, node) = xyz[1 + 3 * node];
      ref_node_xyz(ref_node, 2, node) = xyz[2 + 3 * node];
    }

    RSS(ref_export_by_extension(ref_grid, argv[4]), "bent");

    return 0;
  }

  RXS(ref_args_find(argc, argv, "-u", &pos), REF_NOT_FOUND, "arg");
  if (REF_EMPTY != pos) {
    REF_DBL *scalar;
    REF_INT name_pos = 2;
    REIS(1, pos, "required args: -u id mesh.ext scalar.solb\n");
    REIS(5, argc, "required args: -u id mesh.ext scalar.solb\n");

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_acceptance_u(ref_grid_node(ref_grid), argv[name_pos], scalar),
        "fill u");
    RSS(ref_gather_scalar_by_extension(ref_grid, 1, scalar, NULL, argv[4]),
        "in");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-g", &pos), REF_NOT_FOUND, "arg");
  if (REF_EMPTY != pos) {
    REF_DBL *scalar;
    REF_INT name_pos = 2;
    REIS(1, pos, "required args: -g id mesh.ext scalar.solb\n");
    REIS(5, argc, "required args: -g id mesh.ext scalar.solb\n");

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(scalar, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
    RSS(ref_acceptance_g(ref_grid_node(ref_grid), argv[name_pos], scalar),
        "fill g");
    RSS(ref_gather_scalar_by_extension(ref_grid, 3, scalar, NULL, argv[4]),
        "in");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-q", &pos), REF_NOT_FOUND, "arg");
  if (REF_EMPTY != pos) {
    REF_DBL *scalar;
    REF_INT ldim;
    REF_INT name_pos = 2;
    REIS(1, pos, "required args: -q id mesh.ext scalar.solb\n");
    REIS(5, argc, "required args: -q id mesh.ext scalar.solb\n");

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_acceptance_q(ref_grid_node(ref_grid), argv[name_pos], &ldim,
                         &scalar),
        "fill u");
    if (ref_grid_twod(ref_grid)) {
      if (ref_mpi_once(ref_mpi)) printf("flip field v-w for twod\n");
      RSS(flip_twod_yz(ref_node, ldim, scalar), "flip");
    }
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, NULL, argv[4]),
        "in");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "-pd", &pos), REF_NOT_FOUND, "arg");

  if (REF_EMPTY != pos) {
    REF_DBL *scalar;
    REF_INT ldim;
    REF_INT name_pos = 2;
    REIS(1, pos, "required args: -pq id mesh.ext scalar.solb\n");
    REIS(5, argc, "required args: -pq id mesh.ext scalar.solb\n");

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[3]), "in");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_acceptance_pd(ref_grid_node(ref_grid), argv[name_pos], &ldim,
                          &scalar),
        "fill u");
    if (ref_grid_twod(ref_grid)) {
      if (ref_mpi_once(ref_mpi)) printf("flip field v-w for twod\n");
      RSS(flip_twod_yz(ref_node, ldim, scalar), "flip");
    }
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, NULL, argv[4]),
        "in");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (3 == argc) {
    switch (atoi(argv[1])) {
      case 1:
        RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
        break;
      case 2:
        RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi, 4), "brick");
        break;
      default:
        THROW("case not recognized");
    }
    RSS(ref_export_by_extension(ref_grid, argv[2]), "out");

    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (4 == argc) {
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "in");

    RSS(ref_metric_olympic_node(ref_grid_node(ref_grid), atof(argv[3])), "oly");
    if (ref_grid_twod(ref_grid))
      RSS(ref_metric_twod_node(ref_grid_node(ref_grid)), "2d");

    RSS(ref_gather_metric(ref_grid, argv[2]), "in");

    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (5 == argc) {
    REF_DBL x, z, h, c, r, k;
    REF_INT node;

    c = atof(argv[3]);
    k = atof(argv[4]);

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "in");
    ref_node = ref_grid_node(ref_grid);

    each_ref_node_valid_node(ref_node, node) {
      x = ref_node_xyz(ref_node, 0, node);
      z = ref_node_xyz(ref_node, 2, node);
      r = sqrt(x * x + z * z);
      r = MAX(r, 0.32 * c * c);
      h = c * pow(r, k);

      RSS(ref_node_metric_form(ref_node, node, 1.0 / (h * h), 0, 0, 1.0, 0,
                               1.0 / (h * h)),
          "set linear");
    }

    if (ref_grid_twod(ref_grid))
      RSS(ref_metric_twod_node(ref_grid_node(ref_grid)), "2d");

    RSS(ref_gather_metric(ref_grid, argv[2]), "in");

    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (5 < argc) {
    REF_DBL h0, r;
    REF_DBL z, hx, hz;
    REF_DBL x0, c, k, x, h;
    REF_DBL m0[6], m1[6], m2[6];

    REF_INT node;

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "in");
    ref_node = ref_grid_node(ref_grid);

    pos = 3;
    while (pos < argc) {
      if (strcmp(argv[pos], "-e") == 0) {
        printf("-p\n");
        pos++;
        hx = atof(argv[pos]);
        printf(" hx = %e\n", hx);
        pos++;
        h0 = atof(argv[pos]);
        printf(" h0 = %e\n", h0);
        pos++;
        r = atof(argv[pos]);
        printf(" r = %e\n", r);
        pos++;
        printf(" hz = h0 * (1+r) ** (z/h0)\n");

        each_ref_node_valid_node(ref_node, node) {
          z = ref_node_xyz(ref_node, 2, node);
          hz = h0 * pow(1.0 + r, z / h0);

          RSS(ref_node_metric_form(ref_node, node, 1.0 / (hx * hx), 0, 0, 1.0,
                                   0, 1.0 / (hz * hz)),
              "set linear");
        }

      } else if (strcmp(argv[pos], "-r") == 0) {
        printf("-r\n");
        pos++;
        x0 = atof(argv[pos]);
        printf(" x0 = %e\n", x0);
        pos++;
        c = atof(argv[pos]);
        printf(" c = %e\n", c);
        pos++;
        k = atof(argv[pos]);
        printf(" k = %e\n", k);
        pos++;
        printf(" at x0, h = c * r ** k\n");

        each_ref_node_valid_node(ref_node, node) {
          x = ref_node_xyz(ref_node, 0, node);
          z = ref_node_xyz(ref_node, 2, node);
          r = sqrt((x - x0) * (x - x0) + z * z);
          r = MAX(r, 0.32 * c * c);
          h = c * pow(r, k);

          RSS(ref_node_metric_get(ref_node, node, m0), "get intersect");

          m1[0] = 1.0 / (h * h);
          m1[1] = 0.0;
          m1[2] = 0.0;
          m1[3] = 1.0;
          m1[4] = 0.0;
          m1[5] = 1.0 / (h * h);

          RSS(ref_matrix_intersect(m0, m1, m2), "intersect");
          RSS(ref_node_metric_set(ref_node, node, m2), "set intersect");
        }

      } else if (strcmp(argv[pos], "-h") == 0) {
        printf(" usage\n");
        RSS(ref_grid_free(ref_grid), "grid free");
        RSS(ref_mpi_free(ref_mpi), "mpi free");
        RSS(ref_mpi_stop(), "stop");
        return (0);
      } else {
        fprintf(stderr, "Argument \"%s\" Ignored\n", argv[pos]);
        pos++;
      }
    }

    if (ref_grid_twod(ref_grid))
      RSS(ref_metric_twod_node(ref_grid_node(ref_grid)), "2d");

    RSS(ref_gather_metric(ref_grid, argv[2]), "in");

    RSS(ref_grid_free(ref_grid), "grid free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RSS(ref_grid_free(ref_grid), "grid free");
  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
