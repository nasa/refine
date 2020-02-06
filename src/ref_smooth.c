
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

#include "ref_smooth.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_geom.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"

static REF_STATUS ref_smooth_add_pliant_force(REF_NODE ref_node, REF_INT center,
                                              REF_INT neighbor,
                                              REF_DBL *total_force_vector) {
  REF_INT ixyz;
  REF_DBL norm[3], l4, force, ratio;
  for (ixyz = 0; ixyz < 3; ixyz++)
    norm[ixyz] = ref_node_xyz(ref_node, ixyz, center) -
                 ref_node_xyz(ref_node, ixyz, neighbor);
  RSS(ref_node_ratio(ref_node, center, neighbor, &ratio), "get r0");
  l4 = ratio * ratio * ratio * ratio;
  force = (1.0 - l4) * exp(-l4);
  if (ratio < 1.0) force += (1.0 - ratio) * (1.0 - ratio);
  if (ref_math_divisible(norm[0], ratio) &&
      ref_math_divisible(norm[1], ratio) &&
      ref_math_divisible(norm[2], ratio)) {
    for (ixyz = 0; ixyz < 3; ixyz++) norm[ixyz] /= ratio;
  } else {
    return REF_DIV_ZERO;
  }
  for (ixyz = 0; ixyz < 3; ixyz++)
    total_force_vector[ixyz] += force * norm[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_ratio_around(REF_GRID ref_grid, REF_INT node,
                                       REF_DBL *min_ratio, REF_DBL *max_ratio) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL ratio;

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (node != nodes[cell_node]) {
        RSS(ref_node_ratio(ref_node, node, nodes[cell_node], &ratio), "ratio");
        if (none_found) {
          none_found = REF_FALSE;
          *min_ratio = ratio;
          *max_ratio = ratio;
        } else {
          *min_ratio = MIN(*min_ratio, ratio);
          *max_ratio = MAX(*max_ratio, ratio);
        }
      }
    }
  }

  if (none_found) {
    *min_ratio = 2000.0;
    *max_ratio = -2.0;
    THROW("no triangle found, can not compute ratio");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_quality_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_quality) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL quality;

  *min_quality = 1.0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    none_found = REF_FALSE;
    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual");
    *min_quality = MIN(*min_quality, quality);
  }

  if (none_found) {
    *min_quality = -2.0;
    THROW("no triagle found, can not compute quality");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_normdev_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_normdev) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL normdev;

  *min_normdev = 2.0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    none_found = REF_FALSE;
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "qual");
    *min_normdev = MIN(*min_normdev, normdev);
  }

  if (none_found) {
    *min_normdev = -2.0;
    THROW("no triagle found, can not compute normdev");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_uv_area_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_uv_area) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT id, item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL sign_uv_area, uv_area;

  *min_uv_area = -2.0;

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    id = nodes[ref_cell_node_per(ref_cell)];
    RSS(ref_geom_uv_area_sign(ref_grid, id, &sign_uv_area), "sign");
    RSS(ref_geom_uv_area(ref_grid_geom(ref_grid), nodes, &uv_area), "uv area");
    uv_area *= sign_uv_area;
    if (none_found) {
      *min_uv_area = uv_area;
      none_found = REF_FALSE;
    } else {
      *min_uv_area = MIN(*min_uv_area, uv_area);
    }
  }

  if (none_found) THROW("no triagle found, can not compute min uv area");

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_valid_twod_tri(REF_GRID ref_grid, REF_INT node,
                                     REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL valid;

  if (!ref_grid_twod(ref_grid)) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    if (!valid) return REF_SUCCESS;
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_outward_norm(REF_GRID ref_grid, REF_INT node,
                                   REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL valid;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    if (!valid) return REF_SUCCESS;
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_ideal(REF_GRID ref_grid, REF_INT node, REF_INT tri,
                                REF_DBL *ideal_location) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT n0, n1;
  REF_INT ixyz, i;
  REF_DBL dn[3];
  REF_DBL dt[3];
  REF_DBL log_m0[6], log_m1[6], log_m2[6], log_m[6];
  REF_DBL m[6];
  REF_DBL tangent_length, projection, scale, length_in_metric;

  RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri, nodes), "get tri");
  n0 = REF_EMPTY;
  n1 = REF_EMPTY;
  if (node == nodes[0]) {
    n0 = nodes[1];
    n1 = nodes[2];
  }
  if (node == nodes[1]) {
    n0 = nodes[2];
    n1 = nodes[0];
  }
  if (node == nodes[2]) {
    n0 = nodes[0];
    n1 = nodes[1];
  }
  if (n0 == REF_EMPTY || n1 == REF_EMPTY) THROW("empty triangle side");

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal_location[ixyz] = 0.5 * (ref_node_xyz(ref_node, ixyz, n0) +
                                  ref_node_xyz(ref_node, ixyz, n1));
  for (ixyz = 0; ixyz < 3; ixyz++)
    dn[ixyz] = ref_node_xyz(ref_node, ixyz, node) - ideal_location[ixyz];
  for (ixyz = 0; ixyz < 3; ixyz++)
    dt[ixyz] =
        ref_node_xyz(ref_node, ixyz, n1) - ref_node_xyz(ref_node, ixyz, n0);

  tangent_length = ref_math_dot(dt, dt);
  projection = ref_math_dot(dn, dt);

  if (ref_math_divisible(projection, tangent_length)) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      dn[ixyz] -= (projection / tangent_length) * dt[ixyz];
  } else {
    printf("projection = %e tangent_length = %e\n", projection, tangent_length);
    return REF_DIV_ZERO;
  }

  RSS(ref_math_normalize(dn), "normalize direction");

  /* averaged metric */
  RSS(ref_node_metric_get_log(ref_node, n0, log_m0), "get n0 log m");
  RSS(ref_node_metric_get_log(ref_node, n1, log_m1), "get n1 log m");
  RSS(ref_node_metric_get_log(ref_node, node, log_m2), "get node log m");
  for (i = 0; i < 6; i++) log_m[i] = (log_m0[i] + log_m1[i] + log_m2[i]) / 3.0;
  RSS(ref_matrix_exp_m(log_m, m), "exp avg");

  length_in_metric = ref_matrix_sqrt_vt_m_v(m, dn);

  scale = 0.5 * sqrt(3.0); /* altitude of equilateral triangle */
  if (ref_math_divisible(scale, length_in_metric)) {
    scale = scale / length_in_metric;
  } else {
    printf(" length_in_metric = %e, not invertable\n", length_in_metric);
    return REF_DIV_ZERO;
  }

  for (ixyz = 0; ixyz < 3; ixyz++) ideal_location[ixyz] += scale * dn[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_quality(REF_GRID ref_grid, REF_INT node, REF_INT id,
                                  REF_INT *nodes, REF_DBL *uv, REF_DBL *dq_duv,
                                  REF_DBL step, REF_DBL *qnew) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DBL uvnew[2];
  uvnew[0] = uv[0] + step * dq_duv[0];
  uvnew[1] = uv[1] + step * dq_duv[1];

  RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, uvnew), "set uv");
  RSS(ref_geom_constrain(ref_grid, node), "constrain");
  RSS(ref_node_tri_quality(ref_node, nodes, qnew), "qual");

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_ideal_uv(REF_GRID ref_grid, REF_INT node, REF_INT tri,
                                   REF_DBL *ideal_uv) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT n0, n1;
  REF_INT id, geom;
  REF_DBL uv_orig[2];
  REF_DBL uv[2];
  REF_DBL q0, q;
  REF_DBL xyz[3], dxyz_duv[15], dq_dxyz[3], dq_duv[2], dq_duv0[2], dq_duv1[2];
  REF_DBL slope, beta, num, denom;
  REF_DBL step1, step2, step3, q1, q2, q3;
  REF_INT tries, search;
  REF_BOOL verbose = REF_FALSE;

  RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri, nodes), "get tri");
  n0 = REF_EMPTY;
  n1 = REF_EMPTY;
  if (node == nodes[0]) {
    n0 = nodes[1];
    n1 = nodes[2];
  }
  if (node == nodes[1]) {
    n0 = nodes[2];
    n1 = nodes[0];
  }
  if (node == nodes[2]) {
    n0 = nodes[0];
    n1 = nodes[1];
  }
  if (n0 == REF_EMPTY || n1 == REF_EMPTY) {
    THROW("empty triangle side");
  }
  nodes[0] = node;
  nodes[1] = n0;
  nodes[2] = n1;

  RSS(ref_geom_unique_id(ref_geom, node, REF_GEOM_FACE, &id), "id");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_orig), "uv");
  RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom), "geom");

  RSS(ref_node_tri_quality(ref_node, nodes, &q0), "qual");

  uv[0] = uv_orig[0];
  uv[1] = uv_orig[1];
  dq_duv0[0] = 0; /* uninit warning */
  dq_duv0[1] = 0;
  q = q0;
  for (tries = 0; tries < 30 && q < 0.99; tries++) {
    RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, uv), "set uv");
    RSS(ref_geom_constrain(ref_grid, node), "constrain");
    RSS(ref_node_tri_dquality_dnode0(ref_node, nodes, &q, dq_dxyz), "qual");
    RSS(ref_geom_eval(ref_geom, geom, xyz, dxyz_duv), "eval face");
    dq_duv1[0] = dq_dxyz[0] * dxyz_duv[0] + dq_dxyz[1] * dxyz_duv[1] +
                 dq_dxyz[2] * dxyz_duv[2];
    dq_duv1[1] = dq_dxyz[0] * dxyz_duv[3] + dq_dxyz[1] * dxyz_duv[4] +
                 dq_dxyz[2] * dxyz_duv[5];

    if (0 == tries) {
      beta = 0;
      dq_duv[0] = dq_duv1[0];
      dq_duv[1] = dq_duv1[1];
    } else {
      /* fletcher-reeves */
      num = dq_duv1[0] * dq_duv1[0] + dq_duv1[1] * dq_duv1[1];
      denom = dq_duv0[0] * dq_duv0[0] + dq_duv0[1] * dq_duv0[1];
      /* polak-ribiere */
      num = dq_duv1[0] * (dq_duv1[0] - dq_duv0[0]) +
            dq_duv1[1] * (dq_duv1[1] - dq_duv0[1]);
      denom = dq_duv0[0] * dq_duv0[0] + dq_duv0[1] * dq_duv0[1];
      beta = 0;
      if (ref_math_divisible(num, denom)) beta = num / denom;
      beta = MAX(0.0, beta);
      dq_duv[0] = dq_duv1[0] + beta * dq_duv[0];
      dq_duv[1] = dq_duv1[1] + beta * dq_duv[1];
    }
    dq_duv0[0] = dq_duv1[0];
    dq_duv0[1] = dq_duv1[1];

    slope = sqrt(dq_duv[0] * dq_duv[0] + dq_duv[1] * dq_duv[1]);
    step3 = (1.0 - q) / slope;
    step1 = 0;
    step2 = 0.5 * (step1 + step3);
    RSS(ref_smooth_tri_quality(ref_grid, node, id, nodes, uv, dq_duv, step1,
                               &q1),
        "set uv for q1");
    RSS(ref_smooth_tri_quality(ref_grid, node, id, nodes, uv, dq_duv, step2,
                               &q2),
        "set uv for q2");
    RSS(ref_smooth_tri_quality(ref_grid, node, id, nodes, uv, dq_duv, step3,
                               &q3),
        "set uv for q3");
    for (search = 0; search < 15; search++) {
      if (q1 > q3) {
        step3 = step2;
        q3 = q2;
      } else {
        step1 = step2;
        q1 = q2;
      }
      step2 = 0.5 * (step1 + step3);
      RSS(ref_smooth_tri_quality(ref_grid, node, id, nodes, uv, dq_duv, step2,
                                 &q2),
          "set uv for q2");
    }
    RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv), "uv");

    if (verbose && tries > 25) {
      printf(" slow conv %2d    q %f dq_duv1 %f %f\n", tries, q2, dq_duv1[0],
             dq_duv1[1]);
      printf("              step %f dq_duv  %f %f beta %f\n", step2, dq_duv[0],
             dq_duv[1], beta);
    }
  }

  if (verbose && q < 0.99) {
    printf(" bad ideal q %f dq_duv %f %f\n", q, dq_duv[0], dq_duv[1]);
  }

  RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, uv_orig), "set uv");
  RSS(ref_geom_constrain(ref_grid, node), "constrain");

  ideal_uv[0] = uv[0];
  ideal_uv[1] = uv[1];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_weighted_ideal(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *ideal_location) {
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ixyz;
  REF_DBL tri_ideal[3];
  REF_DBL quality, weight, normalization;

  normalization = 0.0;
  for (ixyz = 0; ixyz < 3; ixyz++) ideal_location[ixyz] = 0.0;

  each_ref_cell_having_node(ref_grid_tri(ref_grid), node, item, cell) {
    RSS(ref_smooth_tri_ideal(ref_grid, node, cell, tri_ideal), "tri ideal");
    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell, nodes), "nodes");
    RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
        "tri qual");
    quality = MAX(quality, ref_grid_adapt(ref_grid, smooth_min_quality));
    weight = 1.0 / quality;
    normalization += weight;
    for (ixyz = 0; ixyz < 3; ixyz++)
      ideal_location[ixyz] += weight * tri_ideal[ixyz];
  }

  if (ref_math_divisible(1.0, normalization)) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ideal_location[ixyz] = (1.0 / normalization) * ideal_location[ixyz];
  } else {
    printf("normalization = %e at %e %e %e\n", normalization,
           ref_node_xyz(ref_grid_node(ref_grid), 0, node),
           ref_node_xyz(ref_grid_node(ref_grid), 1, node),
           ref_node_xyz(ref_grid_node(ref_grid), 2, node));
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_pliant_uv(REF_GRID ref_grid, REF_INT node,
                                    REF_DBL *ideal_uv) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT ixyz, id;
  REF_INT max_node = 100, nnode;
  REF_INT node_list[100];
  REF_INT edge;
  REF_DBL total_force[3], dxyz[3], xyz_orig[3], dxyz_duv[15];
  REF_DBL dxyz_duvn[9], duvn_dxyz[9], r[3], s[3], n[3];
  REF_BOOL self_check = REF_FALSE;

  RSS(ref_geom_unique_id(ref_geom, node, REF_GEOM_FACE, &id), "get id");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, ideal_uv), "get uv_orig");

  RSS(ref_cell_node_list_around(ref_grid_tri(ref_grid), node, max_node, &nnode,
                                node_list),
      "node list for edges");

  for (ixyz = 0; ixyz < 3; ixyz++) total_force[ixyz] = 0.0;
  for (edge = 0; edge < nnode; edge++) {
    RSS(ref_smooth_add_pliant_force(ref_node, node, node_list[edge],
                                    total_force),
        "edge");
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    dxyz[ixyz] =
        ref_grid_adapt(ref_grid, smooth_pliant_alpha) * total_force[ixyz];

  RSS(ref_geom_eval_at(ref_geom, REF_GEOM_FACE, id, ideal_uv, xyz_orig,
                       dxyz_duv),
      "eval face derivatives");
  RSS(ref_geom_face_rsn(ref_geom, id, ideal_uv, r, s, n),
      "eval orthonormal face system");

  for (ixyz = 0; ixyz < 3; ixyz++) {
    dxyz_duvn[ixyz + 0] = dxyz_duv[ixyz + 0];
    dxyz_duvn[ixyz + 3] = dxyz_duv[ixyz + 3];
    dxyz_duvn[ixyz + 6] = n[ixyz];
  }

  RSS(ref_matrix_inv_gen(3, dxyz_duvn, duvn_dxyz), "inverse transformation");

  for (ixyz = 0; ixyz < 3; ixyz++) {
    ideal_uv[0] += duvn_dxyz[0 + 3 * ixyz] * dxyz[ixyz];
    ideal_uv[1] += duvn_dxyz[1 + 3 * ixyz] * dxyz[ixyz];
  }

  if (self_check) {
    REF_DBL du, dv, dn, check_dxyz[3];
    du = 0;
    dv = 0;
    dn = 0;
    for (ixyz = 0; ixyz < 3; ixyz++) {
      du += duvn_dxyz[0 + 3 * ixyz] * dxyz[ixyz];
      dv += duvn_dxyz[1 + 3 * ixyz] * dxyz[ixyz];
      dn += duvn_dxyz[2 + 3 * ixyz] * dxyz[ixyz];
    }
    check_dxyz[0] = du * dxyz_duv[0] + dv * dxyz_duv[3] + dn * n[0];
    check_dxyz[1] = du * dxyz_duv[1] + dv * dxyz_duv[4] + dn * n[1];
    check_dxyz[2] = du * dxyz_duv[2] + dv * dxyz_duv[5] + dn * n[2];
    if (ABS(check_dxyz[0] - dxyz[0]) > 1.0e-12 ||
        ABS(check_dxyz[1] - dxyz[1]) > 1.0e-12 ||
        ABS(check_dxyz[2] - dxyz[2]) > 1.0e-12) {
      printf("du %f dv %f dn %f dxyz %f %f %f\n", du, dv, dn, dxyz[0], dxyz[1],
             dxyz[2]);
      printf(" err %e %e %e\n", check_dxyz[0] - dxyz[0],
             check_dxyz[1] - dxyz[1], check_dxyz[2] - dxyz[2]);
      ref_node_location(ref_node, node);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_weighted_ideal_uv(REF_GRID ref_grid, REF_INT node,
                                            REF_DBL *ideal_uv) {
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT iuv;
  REF_DBL tri_uv[2];
  REF_DBL quality, weight, normalization;

  normalization = 0.0;
  for (iuv = 0; iuv < 2; iuv++) ideal_uv[iuv] = 0.0;

  each_ref_cell_having_node(ref_grid_tri(ref_grid), node, item, cell) {
    RSS(ref_smooth_tri_ideal_uv(ref_grid, node, cell, tri_uv), "tri ideal");
    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell, nodes), "nodes");
    RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
        "tri qual");
    quality = MAX(quality, ref_grid_adapt(ref_grid, smooth_min_quality));
    weight = 1.0 / quality;
    normalization += weight;
    for (iuv = 0; iuv < 2; iuv++) ideal_uv[iuv] += weight * tri_uv[iuv];
  }

  if (ref_math_divisible(1.0, normalization)) {
    for (iuv = 0; iuv < 2; iuv++)
      ideal_uv[iuv] = (1.0 / normalization) * ideal_uv[iuv];
  } else {
    printf("normalization = %e\n", normalization);
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_boundary_nodes(REF_GRID ref_grid, REF_INT node,
                                          REF_INT *node0, REF_INT *node1) {
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_INT item, cell, cell_edge, other;

  *node0 = REF_EMPTY;
  *node1 = REF_EMPTY;

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      if (node == ref_cell_e2n(ref_cell, 0, cell_edge, cell)) {
        other = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      } else if (node == ref_cell_e2n(ref_cell, 1, cell_edge, cell)) {
        other = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      } else {
        continue;
      }
      if (REF_EMPTY == *node0) {
        *node0 = other;
        continue;
      }
      if (REF_EMPTY == *node1) {
        *node1 = other;
        continue;
      }
      THROW("found more than two boundary edges");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_tri_improve(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality0, quality, min_ratio, max_ratio;
  REF_INT ixyz;
  REF_BOOL allowed;

  /* can't handle boundaries yet */
  if (!ref_cell_node_empty(ref_grid_qua(ref_grid), node)) return REF_SUCCESS;

  for (ixyz = 0; ixyz < 3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node, ixyz, node);

  RSS(ref_smooth_tri_weighted_ideal(ref_grid, node, ideal), "ideal");

  RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality0), "q");

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ref_node_xyz(ref_node, ixyz, node) =
          backoff * ideal[ixyz] + (1.0 - backoff) * original[ixyz];
    RSS(ref_smooth_outward_norm(ref_grid, node, &allowed), "normals");
    if (allowed) {
      RSS(ref_metric_interpolate_node(ref_grid, node), "interp node");
      RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality), "q");
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
      if ((quality > quality0) &&
          (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio)) &&
          (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio))) {
        return REF_SUCCESS;
      }
    }
    backoff *= 0.5;
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ref_node_xyz(ref_node, ixyz, node) = original[ixyz];
  RSS(ref_metric_interpolate_node(ref_grid, node), "interp");

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_node_same_tangent(REF_GRID ref_grid, REF_INT node,
                                               REF_INT node0, REF_INT node1,
                                               REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL tan0[3], tan1[3];
  REF_DBL dot;
  REF_INT i;

  *allowed = REF_TRUE;
  for (i = 0; i < 3; i++)
    tan0[i] =
        ref_node_xyz(ref_node, i, node) - ref_node_xyz(ref_node, i, node0);
  for (i = 0; i < 3; i++)
    tan1[i] =
        ref_node_xyz(ref_node, i, node1) - ref_node_xyz(ref_node, i, node);

  RSS(ref_math_normalize(tan0), "edge 0 zero length");
  RSB(ref_math_normalize(tan1), "edge 1 zero length",
      { printf("nodes %d %d %d\n", node0, node, node1); });

  dot = ref_math_dot(tan0, tan1);
  if (dot < ref_node_same_normal_tol(ref_node)) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_move_edge_to(REF_GRID ref_grid, REF_INT node1,
                                   REF_DBL *xyz1, REF_INT node2,
                                   REF_DBL *xyz2) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_BOOL geom_node, geom_edge;
  REF_INT id1, geom1, id2, geom2;
  REF_DBL t1, t_orig1, t_target1, t2, t_orig2, t_target2;
  REF_DBL q1, normdev1, normdev_orig1, min_uv_area1;
  REF_DBL q2, normdev2, normdev_orig2, min_uv_area2;
  REF_STATUS interp_status1, interp_status2;
  REF_INT interp_guess1, interp_guess2;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);
  REF_DBL backoff;
  REF_INT tries;

  RSS(ref_geom_is_a(ref_geom, node1, REF_GEOM_NODE, &geom_node), "node check");
  RSS(ref_geom_is_a(ref_geom, node1, REF_GEOM_EDGE, &geom_edge), "edge check");
  RAS(!geom_node, "geom node not allowed");
  RAS(geom_edge, "geom edge required");
  RSS(ref_geom_is_a(ref_geom, node2, REF_GEOM_NODE, &geom_node), "node check");
  RSS(ref_geom_is_a(ref_geom, node2, REF_GEOM_EDGE, &geom_edge), "edge check");
  RAS(!geom_node, "geom node not allowed");
  RAS(geom_edge, "geom edge required");

  RSS(ref_geom_unique_id(ref_geom, node1, REF_GEOM_EDGE, &id1), "get id");
  RSS(ref_geom_unique_id(ref_geom, node2, REF_GEOM_EDGE, &id2), "get id");
  RSS(ref_geom_find(ref_geom, node1, REF_GEOM_EDGE, id1, &geom1), "get geom");
  RSS(ref_geom_find(ref_geom, node2, REF_GEOM_EDGE, id2, &geom2), "get geom");

  RSS(ref_geom_tuv(ref_geom, node1, REF_GEOM_EDGE, id1, &t_orig1),
      "get t_orig");
  RSS(ref_geom_tuv(ref_geom, node2, REF_GEOM_EDGE, id2, &t_orig2),
      "get t_orig");
  RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id1, xyz1, &t_target1),
      "inv");
  RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id2, xyz2, &t_target2),
      "inv");

  RSS(ref_smooth_tri_normdev_around(ref_grid, node1, &normdev_orig1),
      "nd_orig");
  RSS(ref_smooth_tri_normdev_around(ref_grid, node2, &normdev_orig2),
      "nd_orig");

  interp_guess1 = REF_EMPTY;
  interp_guess2 = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess1 = ref_interp_cell(ref_interp, node1);
      interp_guess2 = ref_interp_cell(ref_interp, node2);
    }
  }

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    t1 = backoff * t_target1 + (1.0 - backoff) * t_orig1;
    t2 = backoff * t_target2 + (1.0 - backoff) * t_orig2;

    RSS(ref_geom_add(ref_geom, node1, REF_GEOM_EDGE, id1, &t1), "set t");
    RSS(ref_geom_constrain(ref_grid, node1), "constrain");
    interp_status1 = ref_metric_interpolate_node(ref_grid, node1);
    RXS(interp_status1, REF_NOT_FOUND, "ref_metric_interpolate_node failed");
    if (ref_grid_surf(ref_grid)) {
      q1 = 1.0;
    } else {
      RSS(ref_smooth_tet_quality_around(ref_grid, node1, &q1), "q");
    }
    RSS(ref_smooth_tri_normdev_around(ref_grid, node1, &normdev1), "nd");
    RSS(ref_smooth_tri_uv_area_around(ref_grid, node1, &min_uv_area1), "a");

    RSS(ref_geom_add(ref_geom, node2, REF_GEOM_EDGE, id2, &t2), "set t");
    RSS(ref_geom_constrain(ref_grid, node2), "constrain");
    interp_status2 = ref_metric_interpolate_node(ref_grid, node2);
    RXS(interp_status2, REF_NOT_FOUND, "ref_metric_interpolate_node failed");
    if (ref_grid_surf(ref_grid)) {
      q2 = 1.0;
    } else {
      RSS(ref_smooth_tet_quality_around(ref_grid, node2, &q2), "q");
    }
    RSS(ref_smooth_tri_normdev_around(ref_grid, node2, &normdev2), "nd");
    RSS(ref_smooth_tri_uv_area_around(ref_grid, node2, &min_uv_area2), "a");

    printf("boff %f nd %8.4f area %8.3e nd %8.4f area %8.3e\n", backoff,
           normdev1, min_uv_area1, normdev2, min_uv_area2);
    if ((q1 > 0.1 * ref_grid_adapt(ref_grid, smooth_min_quality)) &&
        (normdev1 > ref_grid_adapt(ref_grid, post_min_normdev) ||
         normdev1 > normdev_orig1) &&
        (min_uv_area1 > ref_node_min_uv_area(ref_node)) &&
        (q2 > 0.1 * ref_grid_adapt(ref_grid, smooth_min_quality)) &&
        (normdev2 > ref_grid_adapt(ref_grid, post_min_normdev) ||
         normdev2 > normdev_orig2) &&
        (min_uv_area2 > ref_node_min_uv_area(ref_node))) {
      return REF_SUCCESS;
    }
    backoff *= 0.5;
    if (REF_EMPTY != interp_guess1 && REF_SUCCESS != interp_status1)
      ref_interp_cell(ref_interp, node1) = interp_guess1;
    if (REF_EMPTY != interp_guess2 && REF_SUCCESS != interp_status2)
      ref_interp_cell(ref_interp, node2) = interp_guess2;
  }

  RSS(ref_geom_add(ref_geom, node1, REF_GEOM_EDGE, id1, &t_orig1), "set t");
  RSS(ref_geom_add(ref_geom, node2, REF_GEOM_EDGE, id2, &t_orig2), "set t");
  RSS(ref_geom_constrain(ref_grid, node1), "constrain");
  RSS(ref_geom_constrain(ref_grid, node2), "constrain");
  RXS(ref_metric_interpolate_node(ref_grid, node1), REF_NOT_FOUND, "interp");
  RXS(ref_metric_interpolate_node(ref_grid, node2), REF_NOT_FOUND, "interp");

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_bound_improve(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node0, node1;
  REF_INT tries;
  REF_DBL total_force[3];
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality, min_ratio, max_ratio;
  REF_INT ixyz;
  REF_BOOL allowed, geom_edge;

  /* boundaries only */
  if (ref_cell_node_empty(ref_grid_edg(ref_grid), node)) return REF_SUCCESS;
  /* protect mixed-element quads */
  if (!ref_cell_node_empty(ref_grid_qua(ref_grid), node)) return REF_SUCCESS;
  /* protect geometry */
  RSS(ref_geom_is_a(ref_grid_geom(ref_grid), node, REF_GEOM_EDGE, &geom_edge),
      "edge check");
  if (geom_edge) return REF_SUCCESS;

  RSS(ref_smooth_twod_boundary_nodes(ref_grid, node, &node0, &node1),
      "edge nodes");
  if (REF_EMPTY == node1) return REF_SUCCESS;
  RSS(ref_smooth_node_same_tangent(ref_grid, node, node0, node1, &allowed),
      "tan");
  if (!allowed) return REF_SUCCESS;

  for (ixyz = 0; ixyz < 3; ixyz++) total_force[ixyz] = 0.0;
  RSS(ref_smooth_add_pliant_force(ref_node, node, node0, total_force), "n0");
  RSS(ref_smooth_add_pliant_force(ref_node, node, node1, total_force), "n1");

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal[ixyz] =
        ref_node_xyz(ref_node, ixyz, node) +
        ref_grid_adapt(ref_grid, smooth_pliant_alpha) * total_force[ixyz];

  for (ixyz = 0; ixyz < 3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node, ixyz, node);

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ref_node_xyz(ref_node, ixyz, node) =
          backoff * ideal[ixyz] + (1.0 - backoff) * original[ixyz];
    RSS(ref_smooth_outward_norm(ref_grid, node, &allowed), "normals");
    if (allowed) {
      RSS(ref_metric_interpolate_node(ref_grid, node), "interp node");
      RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality), "q");
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
      if (quality > ref_grid_adapt(ref_grid, smooth_min_quality)) {
        return REF_SUCCESS;
      }
    }
    backoff *= 0.5;
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ref_node_xyz(ref_node, ixyz, node) = original[ixyz];
  RSS(ref_metric_interpolate_node(ref_grid, node), "interp");

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_tri_pliant(REF_GRID ref_grid, REF_INT node,
                                        REF_DBL *ideal_location) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT ixyz;
  REF_INT max_node = 100, nnode;
  REF_INT node_list[100];
  REF_INT edge;
  REF_DBL total_force[3];

  RSS(ref_cell_node_list_around(ref_grid_tri(ref_grid), node, max_node, &nnode,
                                node_list),
      "node list for edges");

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal_location[ixyz] = ref_node_xyz(ref_node, ixyz, node);

  for (ixyz = 0; ixyz < 3; ixyz++) total_force[ixyz] = 0.0;
  for (edge = 0; edge < nnode; edge++) {
    RSS(ref_smooth_add_pliant_force(ref_node, node, node_list[edge],
                                    total_force),
        "edge");
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal_location[ixyz] +=
        ref_grid_adapt(ref_grid, smooth_pliant_alpha) * total_force[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_tri_pliant(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality0, quality, min_ratio, max_ratio;
  REF_INT ixyz;
  REF_BOOL allowed;

  /* can't handle boundaries yet */
  if (!ref_cell_node_empty(ref_grid_edg(ref_grid), node)) return REF_SUCCESS;

  for (ixyz = 0; ixyz < 3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node, ixyz, node);

  RSS(ref_smooth_tri_pliant(ref_grid, node, ideal), "ideal");

  RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality0), "q");

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ref_node_xyz(ref_node, ixyz, node) =
          backoff * ideal[ixyz] + (1.0 - backoff) * original[ixyz];
    RSS(ref_smooth_outward_norm(ref_grid, node, &allowed), "normals");
    if (allowed) {
      RSS(ref_metric_interpolate_node(ref_grid, node), "interp node");
      RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality), "q");
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
      if ((quality > 0.9 * quality0 && quality > 0.4) &&
          (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio)) &&
          (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio))) {
        return REF_SUCCESS;
      }
    }
    backoff *= 0.5;
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ref_node_xyz(ref_node, ixyz, node) = original[ixyz];
  RSS(ref_metric_interpolate_node(ref_grid, node), "interp");

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_node_same_normal(REF_GRID ref_grid, REF_INT node,
                                              REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normal[3], first_normal[3];
  REF_DBL dot;
  REF_STATUS status;
  REF_BOOL first_tri;

  *allowed = REF_TRUE;

  first_tri = REF_TRUE;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_node_tri_normal(ref_node, nodes, normal), "orig normal");
    if (first_tri) {
      first_tri = REF_FALSE;
      first_normal[0] = normal[0];
      first_normal[1] = normal[1];
      first_normal[2] = normal[2];
      RSS(ref_math_normalize(first_normal), "original triangle has zero area");
    }
    status = ref_math_normalize(normal);
    if (REF_DIV_ZERO == status) { /* new triangle face has zero area */
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
    RSS(status, "new normal length");
    dot = ref_math_dot(first_normal, normal);
    if (dot < ref_node_same_normal_tol(ref_node)) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_no_geom_tri_improve(REF_GRID ref_grid,
                                                 REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, tri_quality0, tri_quality, tet_quality, min_ratio, max_ratio;
  REF_INT ixyz;
  REF_INT n_ids, ids[2];
  REF_BOOL allowed, geom_face;
  REF_STATUS interp_status;
  REF_INT interp_guess;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);

  /* can't handle mixed elements */
  if (!ref_cell_node_empty(ref_grid_qua(ref_grid), node)) return REF_SUCCESS;

  /* don't move edge nodes */
  RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node, 2, &n_ids, ids),
      REF_INCREASE_LIMIT, "count faceids");
  if (n_ids > 1) return REF_SUCCESS;

  /* not for nodes with geom support */
  RSS(ref_geom_is_a(ref_grid_geom(ref_grid), node, REF_GEOM_FACE, &geom_face),
      "face check");
  if (geom_face) return REF_SUCCESS;

  RSS(ref_smooth_node_same_normal(ref_grid, node, &allowed), "normal dev");
  if (!allowed) return REF_SUCCESS;

  for (ixyz = 0; ixyz < 3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node, ixyz, node);
  interp_guess = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess = ref_interp_cell(ref_interp, node);
    }
  }

  RSS(ref_smooth_tri_quality_around(ref_grid, node, &tri_quality0), "q");
  RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
      "ratio");

  if (tri_quality0 < 0.5 || min_ratio < 0.5 || max_ratio > 2.0) {
    RSS(ref_smooth_tri_weighted_ideal(ref_grid, node, ideal), "ideal");
  } else {
    RSS(ref_smooth_tri_pliant(ref_grid, node, ideal), "ideal");
  }

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ref_node_xyz(ref_node, ixyz, node) =
          backoff * ideal[ixyz] + (1.0 - backoff) * original[ixyz];
    interp_status = ref_metric_interpolate_node(ref_grid, node);
    RXS(interp_status, REF_NOT_FOUND, "ref_metric_interpolate_node failed");
    if (REF_SUCCESS == interp_status) {
      RSS(ref_smooth_valid_twod_tri(ref_grid, node, &allowed), "twod tri");
      RSS(ref_smooth_tri_quality_around(ref_grid, node, &tri_quality), "q");
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
      if (allowed && (tri_quality > tri_quality0) &&
          (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio)) &&
          (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio))) {
        if (ref_cell_node_empty(ref_grid_tet(ref_grid), node)) {
          return REF_SUCCESS;
        } else {
          RSS(ref_smooth_tet_quality_around(ref_grid, node, &tet_quality), "q");
          RSS(ref_smooth_tet_ratio_around(ref_grid, node, &min_ratio,
                                          &max_ratio),
              "ratio");
          if ((REF_SUCCESS == interp_status) &&
              (tet_quality > ref_grid_adapt(ref_grid, smooth_min_quality)) &&
              (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio)) &&
              (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio))) {
            return REF_SUCCESS;
          }
        }
      }
    }
    backoff *= 0.5;
    if (REF_EMPTY != interp_guess && REF_SUCCESS != interp_status)
      ref_interp_cell(ref_interp, node) = interp_guess;
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ref_node_xyz(ref_node, ixyz, node) = original[ixyz];
  RXS(ref_metric_interpolate_node(ref_grid, node), REF_NOT_FOUND, "interp");

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_local_cell_about(REF_CELL ref_cell,
                                              REF_NODE ref_node,
                                              REF_INT about_node,
                                              REF_BOOL *allowed) {
  REF_INT item, cell, node;

  *allowed = REF_FALSE;

  each_ref_cell_having_node(ref_cell, about_node, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }
  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_BOOL allowed;
  REF_DBL quality, min_ratio, max_ratio;

  /* boundary */
  each_ref_node_valid_node(ref_node, node) {
    /* boundaries only */
    allowed = ref_cell_node_empty(ref_grid_edg(ref_grid), node);
    if (allowed) continue;

    RSS(ref_smooth_local_cell_about(ref_grid_pri(ref_grid), ref_node, node,
                                    &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }

    ref_node_age(ref_node, node) = 0;
    RSS(ref_smooth_twod_bound_improve(ref_grid, node), "improve");
  }

  /* interior */
  each_ref_node_valid_node(ref_node, node) {
    /* already did boundaries */
    allowed = ref_cell_node_empty(ref_grid_edg(ref_grid), node);
    if (!allowed) continue;

    RSS(ref_smooth_local_cell_about(ref_grid_pri(ref_grid), ref_node, node,
                                    &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }

    ref_node_age(ref_node, node) = 0;
    RSS(ref_smooth_tri_quality_around(ref_grid, node, &quality), "q");
    RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
        "ratio");
    if (quality < 0.5 || min_ratio < 0.5 || max_ratio > 2.0) {
      RSS(ref_smooth_twod_tri_improve(ref_grid, node), "improve");
    } else {
      RSS(ref_smooth_twod_tri_pliant(ref_grid, node), "improve");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_quality_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_quality) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL quality;

  *min_quality = 1.0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    none_found = REF_FALSE;
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    *min_quality = MIN(*min_quality, quality);
  }

  if (none_found) {
    *min_quality = -2.0;
    return REF_NOT_FOUND;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_ratio_around(REF_GRID ref_grid, REF_INT node,
                                       REF_DBL *min_ratio, REF_DBL *max_ratio) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL ratio;

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (node != nodes[cell_node]) {
        RSS(ref_node_ratio(ref_node, node, nodes[cell_node], &ratio), "ratio");
        if (none_found) {
          none_found = REF_FALSE;
          *min_ratio = ratio;
          *max_ratio = ratio;
        } else {
          *min_ratio = MIN(*min_ratio, ratio);
          *max_ratio = MAX(*max_ratio, ratio);
        }
      }
    }
  }

  if (none_found) {
    *min_ratio = 2000.0;
    *max_ratio = -2.0;
    THROW("no tet found, can not compute ratio");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_ideal(REF_GRID ref_grid, REF_INT node, REF_INT tet,
                                REF_DBL *ideal_location) {
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_node, tri_node;
  REF_INT ixyz;
  REF_DBL dn[3];
  REF_DBL log_m0[6], log_m1[6], log_m2[6], log_m3[6], log_m[6];
  REF_DBL m[6];
  REF_DBL scale, length_in_metric;
  REF_INT i;

  RSS(ref_cell_nodes(ref_cell, tet, nodes), "get tri");
  tri_nodes[0] = REF_EMPTY;
  tri_nodes[1] = REF_EMPTY;
  tri_nodes[2] = REF_EMPTY;
  for (tet_node = 0; tet_node < 4; tet_node++)
    if (node == nodes[tet_node]) {
      for (tri_node = 0; tri_node < 3; tri_node++)
        tri_nodes[tri_node] =
            nodes[ref_cell_f2n_gen(ref_cell, tri_node, tet_node)];
    }
  if (tri_nodes[0] == REF_EMPTY || tri_nodes[1] == REF_EMPTY ||
      tri_nodes[2] == REF_EMPTY)
    THROW("empty tetrahedra face");

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal_location[ixyz] = (ref_node_xyz(ref_node, ixyz, tri_nodes[0]) +
                            ref_node_xyz(ref_node, ixyz, tri_nodes[1]) +
                            ref_node_xyz(ref_node, ixyz, tri_nodes[2])) /
                           3.0;

  RSS(ref_node_tri_normal(ref_node, tri_nodes, dn), "tri normal");

  RSS(ref_math_normalize(dn), "normalize direction");

  /* averaged metric */
  RSS(ref_node_metric_get_log(ref_node, nodes[0], log_m0), "get n0 log m");
  RSS(ref_node_metric_get_log(ref_node, nodes[1], log_m1), "get n1 log m");
  RSS(ref_node_metric_get_log(ref_node, nodes[2], log_m2), "get n2 log m");
  RSS(ref_node_metric_get_log(ref_node, nodes[3], log_m3), "get n3 log m");
  for (i = 0; i < 6; i++)
    log_m[i] = (log_m0[i] + log_m1[i] + log_m2[i] + log_m3[i]) * 0.25;
  RSS(ref_matrix_exp_m(log_m, m), "exp avg");

  length_in_metric = ref_matrix_sqrt_vt_m_v(m, dn);

  scale = sqrt(6.0) / 3.0; /* altitude of regular tetrahedra */
  if (ref_math_divisible(scale, length_in_metric)) {
    scale = scale / length_in_metric;
  } else {
    printf(" length_in_metric = %e, not invertable\n", length_in_metric);
    return REF_DIV_ZERO;
  }

  for (ixyz = 0; ixyz < 3; ixyz++) ideal_location[ixyz] += scale * dn[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_weighted_ideal(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *ideal_location) {
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ixyz;
  REF_DBL tet_ideal[3];
  REF_DBL quality, weight, normalization;

  normalization = 0.0;
  for (ixyz = 0; ixyz < 3; ixyz++) ideal_location[ixyz] = 0.0;

  each_ref_cell_having_node(ref_grid_tet(ref_grid), node, item, cell) {
    RSS(ref_smooth_tet_ideal(ref_grid, node, cell, tet_ideal), "tet ideal");
    RSS(ref_cell_nodes(ref_grid_tet(ref_grid), cell, nodes), "nodes");
    RSS(ref_node_tet_quality(ref_grid_node(ref_grid), nodes, &quality),
        "tet qual");
    quality = MAX(quality, ref_grid_adapt(ref_grid, smooth_min_quality));
    weight = 1.0 / quality;
    normalization += weight;
    for (ixyz = 0; ixyz < 3; ixyz++)
      ideal_location[ixyz] += weight * tet_ideal[ixyz];
  }

  if (ref_math_divisible(1.0, normalization)) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ideal_location[ixyz] = (1.0 / normalization) * ideal_location[ixyz];
  } else {
    printf("normalization = %e\n", normalization);
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_pliant(REF_GRID ref_grid, REF_INT node,
                                 REF_DBL *ideal_location) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT ixyz;
  REF_INT max_node = 100, nnode;
  REF_INT node_list[100];
  REF_INT edge;
  REF_DBL total_force[3];
  RSS(ref_cell_node_list_around(ref_grid_tet(ref_grid), node, max_node, &nnode,
                                node_list),
      "node list for pliant tet edges");

  for (ixyz = 0; ixyz < 3; ixyz++) total_force[ixyz] = 0.0;
  for (edge = 0; edge < nnode; edge++) {
    RSS(ref_smooth_add_pliant_force(ref_node, node, node_list[edge],
                                    total_force),
        "edge");
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ideal_location[ixyz] =
        ref_node_xyz(ref_node, ixyz, node) +
        ref_grid_adapt(ref_grid, smooth_pliant_alpha) * total_force[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_improve(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality0, quality, min_ratio, max_ratio;
  REF_INT ixyz;
  REF_STATUS interp_status;
  REF_INT interp_guess;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);
  REF_BOOL pliant_smoothing;
  REF_BOOL accept;

  /* can't handle boundaries yet */
  if (!ref_cell_node_empty(ref_grid_tri(ref_grid), node) ||
      !ref_cell_node_empty(ref_grid_qua(ref_grid), node))
    return REF_SUCCESS;

  for (ixyz = 0; ixyz < 3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node, ixyz, node);
  interp_guess = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess = ref_interp_cell(ref_interp, node);
    }
  }

  RSS(ref_smooth_tet_quality_around(ref_grid, node, &quality0), "q");
  RSS(ref_smooth_tet_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
      "ratio");
  pliant_smoothing = (quality0 > 0.5 && min_ratio > 0.5 && max_ratio < 2.0);

  if (pliant_smoothing) {
    RSS(ref_smooth_tet_pliant(ref_grid, node, ideal), "pliant");
  } else {
    RSS(ref_smooth_tet_weighted_ideal(ref_grid, node, ideal), "ideal");
  }
  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (ixyz = 0; ixyz < 3; ixyz++)
      ref_node_xyz(ref_node, ixyz, node) =
          backoff * ideal[ixyz] + (1.0 - backoff) * original[ixyz];
    interp_status = ref_metric_interpolate_node(ref_grid, node);
    RXS(interp_status, REF_NOT_FOUND, "ref_metric_interpolate_node failed");
    RSS(ref_smooth_tet_quality_around(ref_grid, node, &quality), "q");
    RSS(ref_smooth_tet_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
        "ratio");
    accept = (REF_SUCCESS == interp_status);
    accept = accept && (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio));
    accept = accept && (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio));
    if (pliant_smoothing) {
      accept = accept && (quality > 0.9 * quality0);
      accept = accept && (quality > 0.4);
    } else {
      accept = accept && (quality > quality0);
    }
    if (accept) {
      return REF_SUCCESS;
    }
    backoff *= 0.5;
    if (REF_EMPTY != interp_guess && REF_SUCCESS != interp_status)
      ref_interp_cell(ref_interp, node) = interp_guess;
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    ref_node_xyz(ref_node, ixyz, node) = original[ixyz];
  RXS(ref_metric_interpolate_node(ref_grid, node), REF_NOT_FOUND, "interp");

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_geom_edge(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL edg = ref_grid_edg(ref_grid);
  REF_BOOL geom_node, geom_node0, geom_node1, geom_edge;
  REF_INT id;
  REF_INT nodes[2], nnode;
  REF_DBL t_orig, t0, t1;
  REF_DBL r0, r1;
  REF_DBL q_orig;
  REF_DBL normdev_orig, normdev;
  REF_DBL min_uv_area;

  REF_INT ixyz;
  REF_DBL total_force[3];
  REF_DBL dxyz[3], dxyz_dt[6], xyz_orig[3], dt, dt_ds, tangent[3];

  REF_DBL t, q, backoff, t_target, min_ratio, max_ratio;
  REF_INT tries;
  REF_BOOL verbose = REF_FALSE;
  REF_INT edge_nodes[REF_CELL_MAX_SIZE_PER], sense;

  REF_STATUS interp_status;
  REF_INT interp_guess;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node check");
  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge check");
  RAS(!geom_node, "geom node not allowed");
  RAS(geom_edge, "geom edge required");

  RSS(ref_geom_unique_id(ref_geom, node, REF_GEOM_EDGE, &id), "get id");
  RSS(ref_cell_node_list_around(edg, node, 2, &nnode, nodes), "edge neighbors");
  REIS(2, nnode, "expected two nodes");

  /* if node0 or node1 constrainted by geoemetry, give precedence */
  RSS(ref_geom_is_a(ref_geom, nodes[0], REF_GEOM_NODE, &geom_node0), "node0");
  RSS(ref_geom_is_a(ref_geom, nodes[1], REF_GEOM_NODE, &geom_node1), "node1");

  for (ixyz = 0; ixyz < 3; ixyz++) total_force[ixyz] = 0.0;
  if (!geom_node1 || (geom_node0 && geom_node1))
    RSS(ref_smooth_add_pliant_force(ref_node, node, nodes[0], total_force),
        "n0");
  if (!geom_node0 || (geom_node0 && geom_node1))
    RSS(ref_smooth_add_pliant_force(ref_node, node, nodes[1], total_force),
        "n1");

  for (ixyz = 0; ixyz < 3; ixyz++)
    dxyz[ixyz] =
        ref_grid_adapt(ref_grid, smooth_pliant_alpha) * total_force[ixyz];

  edge_nodes[0] = nodes[0];
  edge_nodes[1] = node;
  edge_nodes[2] = id;
  RSB(ref_geom_cell_tuv(ref_geom, nodes[0], edge_nodes, REF_GEOM_EDGE, &t0,
                        &sense),
      "get t0", {
        ref_node_location(ref_node, nodes[0]);
        ref_geom_tattle(ref_geom, nodes[0]);
      });
  RSB(ref_geom_cell_tuv(ref_geom, node, edge_nodes, REF_GEOM_EDGE, &t_orig,
                        &sense),
      "get t_orig", {
        ref_node_location(ref_node, node);
        ref_geom_tattle(ref_geom, node);
      });
  edge_nodes[0] = nodes[1];
  edge_nodes[1] = node;
  edge_nodes[2] = id;
  RSB(ref_geom_cell_tuv(ref_geom, nodes[1], edge_nodes, REF_GEOM_EDGE, &t1,
                        &sense),
      "get t1", {
        ref_node_location(ref_node, nodes[1]);
        ref_geom_tattle(ref_geom, nodes[1]);
      });

  RSS(ref_geom_eval_at(ref_geom, REF_GEOM_EDGE, id, &t_orig, xyz_orig, dxyz_dt),
      "eval edge derivatives");
  for (ixyz = 0; ixyz < 3; ixyz++) tangent[ixyz] = dxyz_dt[ixyz];
  dt_ds = sqrt(ref_math_dot(dxyz_dt, dxyz_dt));
  RSS(ref_math_normalize(tangent), "form tangent");
  dt = ref_math_dot(dxyz, tangent);
  if (ref_math_divisible(dt, dt_ds)) {
    dt /= dt_ds;
  } else {
    RSS(REF_DIV_ZERO, "unable to invert dt/dxyz");
  }
  t_target = t_orig + dt;

  /* reject a smooth point canidate outside of trange */
  if (t_target < MIN(t0, t1) || MAX(t0, t1) < t_target) {
    return REF_SUCCESS;
  }

  if (verbose) {
    printf("dxyz %f %f %f\n", dxyz[0], dxyz[1], dxyz[2]);
    printf("dxyz_dt %f %f %f\n", dxyz_dt[0], dxyz_dt[1], dxyz_dt[2]);
  }

  if (ref_grid_surf(ref_grid)) {
    q_orig = 1.0;
  } else {
    RSS(ref_smooth_tet_quality_around(ref_grid, node, &q_orig), "q_orig");
  }
  RSS(ref_smooth_tri_normdev_around(ref_grid, node, &normdev_orig), "nd_orig");
  interp_guess = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess = ref_interp_cell(ref_interp, node);
    }
  }

  if (verbose) printf("edge %d t %f %f %f q %f\n", id, t0, t_orig, t1, q_orig);

  if (verbose)
    printf("t_target %f at %f %f %f\n", t_target,
           ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
           ref_node_xyz(ref_node, 2, node));

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    t = backoff * t_target + (1.0 - backoff) * t_orig;

    RSS(ref_geom_add(ref_geom, node, REF_GEOM_EDGE, id, &t), "set t");
    RSS(ref_geom_constrain(ref_grid, node), "constrain");
    interp_status = ref_metric_interpolate_node(ref_grid, node);
    RXS(interp_status, REF_NOT_FOUND, "ref_metric_interpolate_node failed");
    RSS(ref_node_ratio(ref_node, nodes[0], node, &r0), "get r0");
    RSS(ref_node_ratio(ref_node, nodes[1], node, &r1), "get r1");
    if (ref_grid_surf(ref_grid)) {
      q = 1.0;
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
    } else {
      RSS(ref_smooth_tet_quality_around(ref_grid, node, &q), "q");
      RSS(ref_smooth_tet_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
    }
    RSS(ref_smooth_tri_normdev_around(ref_grid, node, &normdev), "nd");
    RSS(ref_smooth_tri_uv_area_around(ref_grid, node, &min_uv_area), "a");

    if (verbose) printf("t %f r %f %f q %f \n", t, r0, r1, q);
    if ((q > ref_grid_adapt(ref_grid, smooth_min_quality)) &&
        (normdev > ref_grid_adapt(ref_grid, post_min_normdev) ||
         normdev > normdev_orig) &&
        (min_uv_area > ref_node_min_uv_area(ref_node))) {
      return REF_SUCCESS;
    }
    backoff *= 0.5;
    if (REF_EMPTY != interp_guess && REF_SUCCESS != interp_status)
      ref_interp_cell(ref_interp, node) = interp_guess;
  }

  RSS(ref_geom_add(ref_geom, node, REF_GEOM_EDGE, id, &t_orig), "set t");
  RSS(ref_geom_constrain(ref_grid, node), "constrain");
  RXS(ref_metric_interpolate_node(ref_grid, node), REF_NOT_FOUND, "interp");

  if (ref_grid_surf(ref_grid)) {
    q = 1.0;
  } else {
    RSS(ref_smooth_tet_quality_around(ref_grid, node, &q), "q");
  }
  if (verbose) printf("undo q %f\n", q);

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_geom_face(REF_GRID ref_grid, REF_INT node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_BOOL geom_node, geom_edge, geom_face, no_quads;
  REF_INT id;
  REF_DBL uv_orig[2], uv_ideal[2];
  REF_DBL qtet_orig, qtri_orig;
  REF_DBL qtri, qtet, min_uv_area, min_ratio, max_ratio;
  REF_DBL normdev_orig, normdev;
  REF_DBL backoff, uv[2];
  REF_INT tries, iuv;
  REF_DBL uv_min[2], uv_max[2];
  REF_STATUS interp_status;
  REF_INT interp_guess;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);
  REF_BOOL pliant_smoothing = REF_FALSE;
  REF_BOOL accept;

  REF_BOOL verbose = REF_FALSE;

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node check");
  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge check");
  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &geom_face), "face check");
  no_quads = ref_cell_node_empty(ref_grid_qua(ref_grid), node);
  RAS(!geom_node, "geom node not allowed");
  RAS(!geom_edge, "geom edge not allowed");
  RAS(geom_face, "geom face required");
  RAS(no_quads, "quads not allowed");

  RSB(ref_geom_unique_id(ref_geom, node, REF_GEOM_FACE, &id), "get-id",
      ref_clump_around(ref_grid, node, "get-id.tec"));
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_orig), "get uv_orig");
  if (ref_grid_surf(ref_grid)) {
    qtet_orig = 1.0;
  } else {
    RSS(ref_smooth_tet_quality_around(ref_grid, node, &qtet_orig), "q tet");
  }
  RSS(ref_smooth_tri_quality_around(ref_grid, node, &qtri_orig), "q tri");
  RSS(ref_smooth_tri_normdev_around(ref_grid, node, &normdev_orig), "nd_orig");
  RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
      "ratio");
  pliant_smoothing = (qtri_orig > 0.5 && min_ratio > 0.5 && max_ratio < 2.0);
  interp_guess = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess = ref_interp_cell(ref_interp, node);
    }
  }

  if (verbose)
    printf("uv %f %f tri %f tet %f\n", uv_orig[0], uv_orig[1], qtri_orig,
           qtet_orig);

  if (pliant_smoothing) {
    RSS(ref_smooth_tri_pliant_uv(ref_grid, node, uv_ideal), "ideal");
  } else {
    RSS(ref_smooth_tri_weighted_ideal_uv(ref_grid, node, uv_ideal), "ideal");
  }

  RSS(ref_geom_tri_uv_bounding_box(ref_grid, node, uv_min, uv_max), "bb");

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++) {
    for (iuv = 0; iuv < 2; iuv++)
      uv[iuv] = backoff * uv_ideal[iuv] + (1.0 - backoff) * uv_orig[iuv];

    RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, uv), "set uv");
    RSS(ref_geom_constrain(ref_grid, node), "constrain");
    interp_status = ref_metric_interpolate_node(ref_grid, node);
    RXS(interp_status, REF_NOT_FOUND, "ref_metric_interpolate_node failed");

    if (ref_grid_surf(ref_grid)) {
      qtet = 1.0;
      RSS(ref_smooth_tri_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
    } else {
      RSS(ref_smooth_tet_quality_around(ref_grid, node, &qtet), "q tet");
      RSS(ref_smooth_tet_ratio_around(ref_grid, node, &min_ratio, &max_ratio),
          "ratio");
    }
    RSS(ref_smooth_tri_quality_around(ref_grid, node, &qtri), "q tri");
    RSS(ref_smooth_tri_normdev_around(ref_grid, node, &normdev), "nd");
    RSS(ref_smooth_tri_uv_area_around(ref_grid, node, &min_uv_area), "a");

    accept = (REF_SUCCESS == interp_status);
    accept = accept && (normdev > ref_grid_adapt(ref_grid, post_min_normdev) ||
                        normdev > normdev_orig);
    accept = accept && (min_uv_area > ref_node_min_uv_area(ref_node));
    accept = accept && (uv_min[0] < uv[0]) && (uv[0] < uv_max[0]);
    accept = accept && (uv_min[1] < uv[1]) && (uv[1] < uv_max[1]);
    accept = accept && (qtet > ref_grid_adapt(ref_grid, smooth_min_quality));
    accept = accept && (min_ratio >= ref_grid_adapt(ref_grid, post_min_ratio));
    accept = accept && (max_ratio <= ref_grid_adapt(ref_grid, post_max_ratio));
    if (pliant_smoothing) {
      accept = accept && (qtri > 0.9 * qtri_orig);
      accept = accept && (qtri > 0.4);
    } else {
      accept = accept && (qtri > qtri_orig);
    }

    if (accept) {
      if (verbose) printf("better qtri %f qtet %f\n", qtri, qtet);
      return REF_SUCCESS;
    }
    backoff *= 0.5;
    if (REF_EMPTY != interp_guess && REF_SUCCESS != interp_status)
      ref_interp_cell(ref_interp, node) = interp_guess;
  }

  RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, uv_orig), "set t");
  RSS(ref_geom_constrain(ref_grid, node), "constrain");
  RXS(ref_metric_interpolate_node(ref_grid, node), REF_NOT_FOUND, "interp");

  if (verbose)
    printf("undo qtri %f qtet %f was %f %f\n", qtri, qtet, qtri_orig,
           qtet_orig);

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_threed_pass(REF_GRID ref_grid) {
  REF_CELL ref_cell;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom, node;
  REF_BOOL allowed, geom_node, geom_edge, geom_face, interior;

  if (ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  /* smooth edges first if we have geom */
  each_ref_geom_edge(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    /* don't move geom nodes */
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node check");
    if (geom_node) continue;
    /* next to ghost node, can't move */
    RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }
    RSS(ref_smooth_geom_edge(ref_grid, node), "ideal node for edge");
    ref_node_age(ref_node, node) = 0;
  }

  /* smooth edges first without geom, for 2D */
  each_ref_node_valid_node(ref_node, node) {
    /* boundaries only */
    allowed = ref_cell_node_empty(ref_grid_edg(ref_grid), node);
    if (allowed) continue;
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge check");
    if (geom_edge) continue;

    RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }

    ref_node_age(ref_node, node) = 0;
    RSS(ref_smooth_twod_bound_improve(ref_grid, node), "improve");
  }

  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "mov edge");

  /* smooth faces if we have geom, but skip edges */
  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    /* don't move geom nodes */
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge check");
    if (geom_edge) continue;
    /* next to ghost node, can't move */
    RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }
    RSS(ref_smooth_geom_face(ref_grid, node), "ideal node for face");
    ref_node_age(ref_node, node) = 0;
  }

  /* smooth faces without geom */
  each_ref_node_valid_node(ref_node, node) {
    if (ref_cell_node_empty(ref_grid_tri(ref_grid), node)) continue;
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &geom_face), "face check");
    if (geom_face) continue;

    RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }
    RSS(ref_smooth_no_geom_tri_improve(ref_grid, node), "no geom smooth");
  }

  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "mov face");

  /* smooth interior */
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
        "para");
    if (!allowed) {
      ref_node_age(ref_node, node)++;
      continue;
    }

    interior = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
               ref_cell_node_empty(ref_grid_qua(ref_grid), node);
    if (interior) {
      RSS(ref_smooth_tet_improve(ref_grid, node), "ideal tet node");
      ref_node_age(ref_node, node) = 0;
    }
  }

  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "mov int");

  /* smooth low quality tets */
  if (!ref_grid_surf(ref_grid)) {
    REF_DBL quality, min_quality = 0.10;
    REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
      if (quality < min_quality) {
        each_ref_cell_cell_node(ref_cell, cell_node) {
          node = nodes[cell_node];
          RSS(ref_smooth_local_cell_about(ref_cell, ref_node, node, &allowed),
              "para");
          if (!allowed) {
            ref_node_age(ref_node, node)++;
            continue;
          }

          interior = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
                     ref_cell_node_empty(ref_grid_qua(ref_grid), node);
          if (interior) {
            RSS(ref_smooth_tet_improve(ref_grid, node), "ideal");
            ref_node_age(ref_node, node) = 0;
          }
        }
      }
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_smooth_threed_post_edge_split(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_BOOL allowed, interior;

  RSS(ref_smooth_local_cell_about(ref_grid_tet(ref_grid), ref_node, node,
                                  &allowed),
      "para");
  if (!allowed) {
    ref_node_age(ref_node, node)++;
    return REF_SUCCESS;
  }

  interior = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
             ref_cell_node_empty(ref_grid_qua(ref_grid), node);
  if (interior) {
    RSS(ref_smooth_tet_improve(ref_grid, node), "ideal tet node");
    ref_node_age(ref_node, node) = 0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_threed_post_face_split(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_BOOL allowed, interior;

  RSS(ref_smooth_local_cell_about(ref_grid_pri(ref_grid), ref_node, node,
                                  &allowed),
      "para");
  if (!allowed) {
    ref_node_age(ref_node, node)++;
    return REF_SUCCESS;
  }

  interior = ref_cell_node_empty(ref_grid_qua(ref_grid), node);
  if (interior) {
    RSS(ref_smooth_twod_tri_improve(ref_grid, node), "ideal tri node");
    ref_node_age(ref_node, node) = 0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_report_quality_around(REF_GRID ref_grid,
                                                REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT item, cell, i;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality;

  i = 0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    printf(" %d:%5.3f", i, quality);
    i++;
  }
  printf("\n");

  return REF_SUCCESS;
}

/* does not have ratio limits */
REF_STATUS ref_smooth_nso_step(REF_GRID ref_grid, REF_INT node,
                               REF_BOOL *complete) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, d_quality[3];
  REF_INT *cells, *active;
  REF_DBL *quals, *grads;
  REF_INT worst, degree;
  REF_DBL min_qual;
  REF_INT i, j, ixyz;
  REF_DBL dir[3];
  REF_DBL m1, m0, alpha, min_alpha;
  REF_INT mate, reductions, max_reductions;
  REF_DBL requirement;
  REF_DBL xyz[3];
  REF_DBL active_tol = 1.0e-12;
  REF_INT nactive;
  REF_DBL last_alpha, last_qual;
  REF_BOOL verbose = REF_FALSE;
  REF_STATUS interp_status;
  REF_INT interp_guess;
  REF_INTERP ref_interp = ref_grid_interp(ref_grid);

  *complete = REF_FALSE;

  xyz[0] = ref_node_xyz(ref_node, 0, node);
  xyz[1] = ref_node_xyz(ref_node, 1, node);
  xyz[2] = ref_node_xyz(ref_node, 2, node);
  interp_guess = REF_EMPTY;
  if (NULL != ref_interp) {
    if (ref_interp_continuously(ref_interp)) {
      interp_guess = ref_interp_cell(ref_interp, node);
    }
  }
  interp_status = REF_SUCCESS;

  RSS(ref_adj_degree(ref_cell_adj(ref_cell), node, &degree), "deg");

  ref_malloc(cells, degree, REF_INT);
  ref_malloc(active, degree, REF_INT);
  ref_malloc(quals, degree, REF_DBL);
  ref_malloc(grads, 3 * degree, REF_DBL);

  min_qual = 1.0;
  worst = REF_EMPTY;
  degree = 0;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_cell_orient_node0(4, node, nodes), "orient");
    RSS(ref_node_tet_dquality_dnode0(ref_node, nodes, &quality, d_quality),
        "qual");
    if (quality < min_qual || REF_EMPTY == worst) {
      worst = degree;
      min_qual = quality;
    }
    cells[degree] = cell;
    quals[degree] = quality;
    grads[0 + 3 * degree] = d_quality[0];
    grads[1 + 3 * degree] = d_quality[1];
    grads[2 + 3 * degree] = d_quality[2];
    degree++;
  }
  active[0] = worst;
  nactive = 1;
  for (i = 0; i < degree; i++) {
    if (i == worst) continue;
    if ((quals[i] - quals[active[0]]) < active_tol) {
      active[nactive] = i;
      nactive++;
    }
  }
  if (verbose)
    for (i = 0; i < nactive; i++)
      printf("%2d: %10.5f %10.5f %10.5f\n", active[i], grads[0 + 3 * active[i]],
             grads[1 + 3 * active[i]], grads[2 + 3 * active[i]]);

  if (4 <= nactive) {
    *complete = REF_TRUE;
    goto success_clean_and_return;
  }

  if (1 == nactive) {
    dir[0] = grads[0 + 3 * worst];
    dir[1] = grads[1 + 3 * worst];
    dir[2] = grads[2 + 3 * worst];
  } else { /* Charalambous and Conn DOI:10.1137/0715011 equ (3.2)  */
    REF_INT k;
    REF_DBL N[16];
    REF_DBL NNt[16];
    REF_DBL invNNt[16];
    REF_DBL NtinvNNt[16];
    REF_DBL NtinvNNtN[16];
    REF_DBL P[16];

    for (i = 0; i < nactive; i++) {
      N[i + nactive * 0] = 1.0;
      for (ixyz = 0; ixyz < 3; ixyz++)
        N[i + nactive * (1 + ixyz)] = -grads[ixyz + 3 * active[i]];
    }
    /* P = I - Nt [N Nt]^-1 N */
    for (i = 0; i < nactive; i++)
      for (j = 0; j < nactive; j++) NNt[i + j * nactive] = 0.0;
    for (i = 0; i < nactive; i++)
      for (j = 0; j < nactive; j++)
        for (k = 0; k < 4; k++)
          NNt[i + j * nactive] += N[i + nactive * k] * N[j + nactive * k];
    RSS(ref_matrix_inv_gen(nactive, NNt, invNNt), "inv");

    for (i = 0; i < 4; i++)
      for (j = 0; j < nactive; j++) NtinvNNt[i + 4 * j] = 0.0;

    for (i = 0; i < 4; i++)
      for (j = 0; j < nactive; j++)
        for (k = 0; k < nactive; k++)
          NtinvNNt[i + 4 * j] += N[k + i * nactive] * invNNt[k + j * nactive];

    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) NtinvNNtN[i + 4 * j] = 0.0;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        for (k = 0; k < nactive; k++)
          NtinvNNtN[i + 4 * j] += NtinvNNt[i + k * 4] * N[k + j * nactive];

    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) P[i + j * 4] = -NtinvNNtN[i + j * 4];
    for (i = 0; i < 4; i++) P[i + i * 4] += 1.0;

    for (ixyz = 0; ixyz < 3; ixyz++) dir[ixyz] = P[1 + ixyz];
  }

  RSS(ref_math_normalize(dir), "norm");

  m0 = ref_math_dot(dir, &(grads[3 * worst]));
  if (verbose)
    for (i = 0; i < nactive; i++)
      printf("slope %e\n", ref_math_dot(dir, &(grads[3 * active[i]])));
  if (0.0 >= m0) {
    printf("%s: %d: %s: m0 not positive %e", __FILE__, __LINE__, __func__, m0);
    *complete = REF_TRUE;
    goto success_clean_and_return;
  }
  mate = REF_EMPTY;
  min_alpha = 1.0e10;
  for (i = 0; i < degree; i++) {
    REF_BOOL skip;
    skip = REF_FALSE;
    for (j = 0; j < nactive; j++)
      if (i == active[j]) skip = REF_TRUE;
    if (skip) continue;
    m1 = ref_math_dot(dir, &(grads[3 * i]));
    /*
      cost = quals[i]+alpha*m1;
      cost = quals[worst]+alpha*m0;
      quals[i]+alpha*m1 = quals[worst]+alpha*m0;
    */
    if (!ref_math_divisible((quals[worst] - quals[i]), (m1 - m0))) continue;
    alpha = (quals[worst] - quals[i]) / (m1 - m0);
    if ((alpha > 0.0 && alpha < min_alpha)) {
      min_alpha = alpha;
      mate = i;
    }
  }
  if (REF_EMPTY == mate) {
    if (verbose) {
      for (i = 0; i < nactive; i++)
        printf("active slope %f %d\n",
               ref_math_dot(dir, &(grads[3 * active[i]])), i);
      for (i = 0; i < degree; i++)
        printf("all %f slope %f %d %e\n", quals[i],
               ref_math_dot(dir, &(grads[3 * i])), i,
               (1.0 - quals[i]) / ref_math_dot(dir, &(grads[3 * i])));
    }
    min_alpha = (1.0 - quals[0]) / ref_math_dot(dir, &(grads[3 * 0]));
    for (i = 1; i < degree; i++)
      min_alpha =
          MIN(min_alpha, (1.0 - quals[i]) / ref_math_dot(dir, &(grads[3 * i])));
  }

  alpha = min_alpha;
  last_alpha = 0.0;
  last_qual = 0.0;
  max_reductions = 8;
  for (reductions = 0; reductions < max_reductions; reductions++) {
    ref_node_xyz(ref_node, 0, node) = xyz[0] + alpha * dir[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1] + alpha * dir[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2] + alpha * dir[2];
    interp_status = ref_metric_interpolate_node(ref_grid, node);
    RXS(interp_status, REF_NOT_FOUND, "ref_metric_interpolate_node failed");

    RSS(ref_smooth_tet_quality_around(ref_grid, node, &quality), "rep");
    requirement = 0.9 * alpha * m0 + quals[worst];
    if (verbose)
      printf(" %d alpha %e min %f required %f actual %f\n", nactive, alpha,
             min_qual, requirement, quality);
    if (REF_SUCCESS == interp_status && reductions > 0 && quality < last_qual &&
        quality > min_qual) {
      if (verbose)
        printf("use last alpha %e min %f last_qual %f actual %f\n", last_alpha,
               min_qual, last_qual, quality);
      alpha = last_alpha;
      quality = last_qual;
      ref_node_xyz(ref_node, 0, node) = xyz[0] + alpha * dir[0];
      ref_node_xyz(ref_node, 1, node) = xyz[1] + alpha * dir[1];
      ref_node_xyz(ref_node, 2, node) = xyz[2] + alpha * dir[2];
      RSS(ref_metric_interpolate_node(ref_grid, node), "interp");
      break;
    }
    if (quality > requirement || alpha < 1.0e-12) break;
    last_alpha = alpha;
    last_qual = quality;
    alpha *= 0.5;
  }

  if (max_reductions <= reductions) { /* used all the reductions, step is small,
                                         marginal gains remain */
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    if (REF_SUCCESS != interp_status && REF_EMPTY != interp_guess)
      ref_interp_cell(ref_interp, node) = interp_guess;
    RXS(ref_metric_interpolate_node(ref_grid, node), REF_NOT_FOUND, "interp");
    *complete = REF_TRUE;
  }

  if (3 == nactive &&
      (quality - min_qual) < 1.0e-5) { /* very small step toward forth active,
                                          marginal gains remain */
    *complete = REF_TRUE;
  }

success_clean_and_return:

  ref_free(grads);
  ref_free(quals);
  ref_free(active);
  ref_free(cells);
  return REF_SUCCESS;
}

REF_STATUS ref_smooth_nso(REF_GRID ref_grid, REF_INT node) {
  REF_BOOL allowed, interior;
  REF_BOOL complete = REF_FALSE;
  REF_INT step;

  RSS(ref_smooth_local_cell_about(ref_grid_tet(ref_grid),
                                  ref_grid_node(ref_grid), node, &allowed),
      "para");
  if (!allowed) return REF_SUCCESS;

  interior = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
             ref_cell_node_empty(ref_grid_qua(ref_grid), node);
  if (!interior) return REF_SUCCESS;

  for (step = 0; step < 100; step++) {
    RSS(ref_smooth_nso_step(ref_grid, node, &complete), "step");
    if (complete) break;
  }

  return REF_SUCCESS;
}
