
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

#include "ref_metric.h"

#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_grid.h"
#include "ref_node.h"

#include "ref_interp.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"

#include "ref_dict.h"

#define REF_METRIC_MAX_DEGREE (1000)

REF_STATUS ref_metric_show(REF_DBL *m) {
  printf(" %18.10e %18.10e %18.10e\n", m[0], m[1], m[2]);
  printf(" %18.10e %18.10e %18.10e\n", m[1], m[3], m[4]);
  printf(" %18.10e %18.10e %18.10e\n", m[2], m[4], m[5]);
  return REF_SUCCESS;
}

REF_STATUS ref_metric_inspect(REF_NODE ref_node) {
  REF_INT node;
  each_ref_node_valid_node(ref_node, node)
      RSS(ref_metric_show(ref_node_metric_ptr(ref_node, node)), "show it");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_from_node(REF_DBL *metric, REF_NODE ref_node) {
  REF_INT node, im;

  each_ref_node_valid_node(ref_node, node) for (im = 0; im < 6; im++)
      metric[im + 6 * node] = ref_node_metric(ref_node, im, node);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_to_node(REF_DBL *metric, REF_NODE ref_node) {
  REF_INT node, im;

  each_ref_node_valid_node(ref_node, node) for (im = 0; im < 6; im++)
      ref_node_metric(ref_node, im, node) = metric[im + 6 * node];

  return REF_SUCCESS;
}

REF_STATUS ref_metric_unit_node(REF_NODE ref_node) {
  REF_INT node;

  each_ref_node_valid_node(ref_node, node) {
    ref_node_metric(ref_node, 0, node) = 1.0;
    ref_node_metric(ref_node, 1, node) = 0.0;
    ref_node_metric(ref_node, 2, node) = 0.0;
    ref_node_metric(ref_node, 3, node) = 1.0;
    ref_node_metric(ref_node, 4, node) = 0.0;
    ref_node_metric(ref_node, 5, node) = 1.0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_olympic_node(REF_NODE ref_node, REF_DBL h) {
  REF_INT node;
  REF_DBL hh;

  each_ref_node_valid_node(ref_node, node) {
    ref_node_metric(ref_node, 0, node) = 1.0 / (0.1 * 0.1);
    ref_node_metric(ref_node, 1, node) = 0.0;
    ref_node_metric(ref_node, 2, node) = 0.0;
    ref_node_metric(ref_node, 3, node) = 1.0 / (0.1 * 0.1);
    ref_node_metric(ref_node, 4, node) = 0.0;
    hh = h + (0.1 - h) * ABS(ref_node_xyz(ref_node, 2, node) - 0.5) / 0.5;
    ref_node_metric(ref_node, 5, node) = 1.0 / (hh * hh);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_ring_node(REF_NODE ref_node) {
  REF_INT node;
  REF_DBL hh;
  REF_DBL h = 0.01;
  REF_DBL x;
  each_ref_node_valid_node(ref_node, node) {
    x = ref_node_xyz(ref_node, 0, node);
    hh = h + (0.1 - h) * MIN(2 * ABS(x - 1.0), 1);
    ref_node_metric(ref_node, 0, node) = 1.0 / (hh * hh);
    ref_node_metric(ref_node, 1, node) = 0.0;
    ref_node_metric(ref_node, 2, node) = 0.0;
    ref_node_metric(ref_node, 3, node) = 1.0 / (0.1 * 0.1);
    ref_node_metric(ref_node, 4, node) = 0.0;
    ref_node_metric(ref_node, 5, node) = 1.0 / (0.1 * 0.1);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_polar2d_node(REF_NODE ref_node) {
  REF_INT node;
  REF_DBL x, z, r, t;
  REF_DBL h_y, h_t, h_r, h0;
  REF_DBL d[12], m[6];

  each_ref_node_valid_node(ref_node, node) {
    x = ref_node_xyz(ref_node, 0, node);
    z = ref_node_xyz(ref_node, 2, node);
    r = sqrt(x * x + z * z);
    t = atan2(z, x);
    h_y = 1.0;
    h_t = 0.1;
    h0 = 0.001;
    h_r = h0 + 2 * (0.1 - h0) * ABS(r - 0.5);
    ref_matrix_eig(d, 0) = 1.0 / (h_r * h_r);
    ref_matrix_eig(d, 1) = 1.0 / (h_t * h_t);
    ref_matrix_eig(d, 2) = 1.0 / (h_y * h_y);
    ref_matrix_vec(d, 0, 0) = cos(t);
    ref_matrix_vec(d, 1, 0) = 0.0;
    ref_matrix_vec(d, 2, 0) = sin(t);
    ref_matrix_vec(d, 0, 1) = -sin(t);
    ref_matrix_vec(d, 1, 1) = 0.0;
    ref_matrix_vec(d, 2, 1) = cos(t);
    ref_matrix_vec(d, 0, 2) = 0.0;
    ref_matrix_vec(d, 1, 2) = 1.0;
    ref_matrix_vec(d, 2, 2) = 0.0;
    ref_matrix_form_m(d, m);
    ref_node_metric(ref_node, 0, node) = m[0];
    ref_node_metric(ref_node, 1, node) = m[1];
    ref_node_metric(ref_node, 2, node) = m[2];
    ref_node_metric(ref_node, 3, node) = m[3];
    ref_node_metric(ref_node, 4, node) = m[4];
    ref_node_metric(ref_node, 5, node) = m[5];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_ugawg_node(REF_NODE ref_node, REF_INT version) {
  REF_INT node;
  REF_DBL x, y, r, t;
  REF_DBL h_z, h_t, h_r, h0, d0;
  REF_DBL d[12], m[6];

  each_ref_node_valid_node(ref_node, node) {
    x = ref_node_xyz(ref_node, 0, node);
    y = ref_node_xyz(ref_node, 1, node);
    r = sqrt(x * x + y * y);
    t = atan2(y, x);
    h_z = 0.1;
    h_t = 0.1;
    h0 = 0.001;
    h_r = h0 + 2 * (0.1 - h0) * ABS(r - 0.5);
    if (2 == version) {
      d0 = MIN(10.0 * ABS(r - 0.5), 1.0);
      h_t = 0.1 * d0 + 0.025 * (1.0 - d0);
    }
    ref_matrix_eig(d, 0) = 1.0 / (h_r * h_r);
    ref_matrix_eig(d, 1) = 1.0 / (h_t * h_t);
    ref_matrix_eig(d, 2) = 1.0 / (h_z * h_z);
    ref_matrix_vec(d, 0, 0) = cos(t);
    ref_matrix_vec(d, 1, 0) = sin(t);
    ref_matrix_vec(d, 2, 0) = 0.0;
    ref_matrix_vec(d, 0, 1) = -sin(t);
    ref_matrix_vec(d, 1, 1) = cos(t);
    ref_matrix_vec(d, 2, 1) = 0.0;
    ref_matrix_vec(d, 0, 2) = 0.0;
    ref_matrix_vec(d, 1, 2) = 0.0;
    ref_matrix_vec(d, 2, 2) = 1.0;
    ref_matrix_form_m(d, m);
    ref_node_metric(ref_node, 0, node) = m[0];
    ref_node_metric(ref_node, 1, node) = m[1];
    ref_node_metric(ref_node, 2, node) = m[2];
    ref_node_metric(ref_node, 3, node) = m[3];
    ref_node_metric(ref_node, 4, node) = m[4];
    ref_node_metric(ref_node, 5, node) = m[5];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_masabl_node(REF_NODE ref_node) {
  REF_INT node;
  REF_DBL hx, hz, c, k1;

  each_ref_node_valid_node(ref_node, node) {
    hx =
        0.01 + 0.2 * cos(ref_math_pi * (ref_node_xyz(ref_node, 0, node) - 0.5));
    ref_node_metric(ref_node, 0, node) = 1.0 / (hx * hx);
    ref_node_metric(ref_node, 1, node) = 0.0;
    ref_node_metric(ref_node, 2, node) = 0.0;
    ref_node_metric(ref_node, 3, node) = 1.0 / (0.1 * 0.1);
    ref_node_metric(ref_node, 4, node) = 0.0;
    c = 0.001;
    k1 = 6.0;
    hz = c * exp(k1 * ref_node_xyz(ref_node, 2, node));
    ref_node_metric(ref_node, 5, node) = 1.0 / (hz * hz);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_twod_node(REF_NODE ref_node) {
  REF_INT node;

  each_ref_node_valid_node(ref_node, node) {
    ref_node_metric(ref_node, 1, node) = 0.0;
    ref_node_metric(ref_node, 3, node) = 1.0;
    ref_node_metric(ref_node, 4, node) = 0.0;
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_metric_interpolate_twod(REF_GRID ref_grid,
                                              REF_GRID parent_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_NODE parent_node = ref_grid_node(parent_grid);
  REF_INT node, tri, ixyz, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[3], interpolated_xyz[3], bary[3];
  REF_DBL tol = 1.0e-11;
  REF_DBL log_parent_m[3][6];
  REF_DBL log_interpolated_m[6];

  if (ref_mpi_para(ref_grid_mpi(ref_grid)))
    RSS(REF_IMPLEMENT, "twod ref_metric_interpolate not para");

  if (!ref_grid_twod(ref_grid))
    RSS(REF_IMPLEMENT, "ref_metric_interpolate only implemented for twod");

  each_ref_node_valid_node(ref_node, node) {
    /* skip mixed element nodes, they can't move */
    if (!ref_cell_node_empty(ref_grid_hex(ref_grid), node)) continue;
    for (ixyz = 0; ixyz < 3; ixyz++)
      xyz[ixyz] = ref_node_xyz(ref_node, ixyz, node);
    tri = REF_EMPTY;
    RSS(ref_grid_enclosing_tri(parent_grid, xyz, &tri, bary), "enclosing tri");
    RSS(ref_cell_nodes(ref_grid_tri(parent_grid), tri, nodes), "c2n");
    for (ixyz = 0; ixyz < 3; ixyz++) {
      interpolated_xyz[ixyz] = 0.0;
      for (ibary = 0; ibary < 3; ibary++)
        interpolated_xyz[ixyz] +=
            bary[ibary] * ref_node_real(parent_node, ixyz, nodes[ibary]);
    }
    /* override y for fake twod */
    interpolated_xyz[1] = ref_node_xyz(ref_node, 1, node);
    for (ixyz = 0; ixyz < 3; ixyz++)
      RWDS(xyz[ixyz], interpolated_xyz[ixyz], tol, "xyz check");
    for (ibary = 0; ibary < 3; ibary++)
      RSS(ref_matrix_log_m(ref_node_metric_ptr(parent_node, nodes[ibary]),
                           log_parent_m[ibary]),
          "log(parentM)");
    for (im = 0; im < 6; im++) {
      log_interpolated_m[im] = 0.0;
      for (ibary = 0; ibary < 3; ibary++)
        log_interpolated_m[im] += bary[ibary] * log_parent_m[ibary][im];
    }
    RSS(ref_matrix_exp_m(log_interpolated_m,
                         ref_node_metric_ptr(ref_node, node)),
        "exp(intrpM)");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_interpolate(REF_GRID to_grid, REF_GRID from_grid) {
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_MPI ref_mpi = ref_grid_mpi(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_INTERP ref_interp;
  REF_INT node, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL max_error, tol = 1.0e-8;
  REF_DBL log_parent_m[4][6];
  REF_DBL log_interpolated_m[6];
  REF_INT receptor, n_recept, donation, n_donor;
  REF_DBL *recept_m, *donor_m, *recept_bary, *donor_bary;
  REF_INT *donor_node, *donor_ret, *donor_cell;
  REF_INT *recept_proc, *recept_ret, *recept_node, *recept_cell;

  if (ref_grid_twod(to_grid)) {
    RSS(ref_metric_interpolate_twod(to_grid, from_grid), "2d version");
    return REF_SUCCESS;
  }

  RSS(ref_interp_create(&ref_interp, from_grid, to_grid), "make interp");
  RSS(ref_interp_locate(ref_interp), "map");

  RSS(ref_interp_max_error(ref_interp, &max_error), "err");
  if (max_error > tol && ref_mpi_once(ref_mpi)) {
    printf("warning: %e max_error greater than %e tol\n", max_error, tol);
  }

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      n_recept++;
    }
  }

  ref_malloc(recept_bary, 4 * n_recept, REF_DBL);
  ref_malloc(recept_cell, n_recept, REF_INT);
  ref_malloc(recept_node, n_recept, REF_INT);
  ref_malloc(recept_ret, n_recept, REF_INT);
  ref_malloc(recept_proc, n_recept, REF_INT);

  n_recept = 0;
  each_ref_node_valid_node(to_node, node) {
    if (ref_node_owned(to_node, node)) {
      RUS(REF_EMPTY, ref_interp->cell[node], "node needs to be localized");
      for (ibary = 0; ibary < 4; ibary++) {
        recept_bary[ibary + 4 * n_recept] = ref_interp->bary[ibary + 4 * node];
      }
      recept_proc[n_recept] = ref_interp->part[node];
      recept_cell[n_recept] = ref_interp->cell[node];
      recept_node[n_recept] = node;
      recept_ret[n_recept] = ref_mpi_rank(ref_mpi);
      n_recept++;
    }
  }

  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_cell, 1, n_recept,
                        (void **)(&donor_cell), &n_donor, REF_INT_TYPE),
      "blind send cell");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_ret, 1, n_recept,
                        (void **)(&donor_ret), &n_donor, REF_INT_TYPE),
      "blind send ret");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_node, 1, n_recept,
                        (void **)(&donor_node), &n_donor, REF_INT_TYPE),
      "blind send node");
  RSS(ref_mpi_blindsend(ref_mpi, recept_proc, (void *)recept_bary, 4, n_recept,
                        (void **)(&donor_bary), &n_donor, REF_DBL_TYPE),
      "blind send bary");

  ref_free(recept_proc);
  ref_free(recept_ret);
  ref_free(recept_node);
  ref_free(recept_cell);
  ref_free(recept_bary);

  ref_malloc(donor_m, 6 * n_donor, REF_DBL);

  for (donation = 0; donation < n_donor; donation++) {
    RSS(ref_cell_nodes(from_cell, donor_cell[donation], nodes),
        "node needs to be localized");
    for (ibary = 0; ibary < 4; ibary++)
      RSS(ref_matrix_log_m(ref_node_metric_ptr(from_node, nodes[ibary]),
                           log_parent_m[ibary]),
          "log(parentM)");
    for (im = 0; im < 6; im++) {
      log_interpolated_m[im] = 0.0;
      for (ibary = 0; ibary < 4; ibary++) {
        log_interpolated_m[im] +=
            donor_bary[ibary + 4 * donation] * log_parent_m[ibary][im];
      }
    }
    RSS(ref_matrix_exp_m(log_interpolated_m, &(donor_m[6 * donation])),
        "exp(intrpM)");
  }
  ref_free(donor_cell);
  ref_free(donor_bary);

  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_m, 6, n_donor,
                        (void **)(&recept_m), &n_recept, REF_DBL_TYPE),
      "blind send bary");
  RSS(ref_mpi_blindsend(ref_mpi, donor_ret, (void *)donor_node, 1, n_donor,
                        (void **)(&recept_node), &n_recept, REF_INT_TYPE),
      "blind send node");
  ref_free(donor_m);
  ref_free(donor_node);
  ref_free(donor_ret);

  for (receptor = 0; receptor < n_recept; receptor++) {
    node = recept_node[receptor];
    for (im = 0; im < 6; im++) {
      ref_node_metric(to_node, im, node) = recept_m[im + 6 * receptor];
    }
  }

  ref_free(recept_node);
  ref_free(recept_m);

  RSS(ref_interp_free(ref_interp), "interp free");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_gradation(REF_DBL *metric, REF_GRID ref_grid, REF_DBL r) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *metric_orig;
  REF_DBL ratio, enlarge, log_r;
  REF_DBL direction[3];
  REF_DBL limit_metric[6], limited[6];
  REF_INT node, i;
  REF_INT edge, node0, node1;

  log_r = log(r);

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc(metric_orig, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < 6; i++) metric_orig[i + 6 * node] = metric[i + 6 * node];
  }

  /* F. Alauzet doi:10.1016/j.finel.2009.06.028 equation (9) */

  each_ref_edge(ref_edge, edge) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    direction[0] =
        (ref_node_xyz(ref_node, 0, node1) - ref_node_xyz(ref_node, 0, node0));
    direction[1] =
        (ref_node_xyz(ref_node, 1, node1) - ref_node_xyz(ref_node, 1, node0));
    direction[2] =
        (ref_node_xyz(ref_node, 2, node1) - ref_node_xyz(ref_node, 2, node0));

    ratio = ref_matrix_sqrt_vt_m_v(&(metric_orig[6 * node1]), direction);
    enlarge = pow(1.0 + ratio * log_r, -2.0);
    for (i = 0; i < 6; i++)
      limit_metric[i] = metric_orig[i + 6 * node1] * enlarge;
    RSS(ref_matrix_intersect(&(metric_orig[6 * node0]), limit_metric, limited),
        "limit m0 with enlarged m1");
    RSS(ref_matrix_intersect(&(metric[6 * node0]), limited,
                             &(metric[6 * node0])),
        "update m0");

    ratio = ref_matrix_sqrt_vt_m_v(&(metric_orig[6 * node0]), direction);
    enlarge = pow(1.0 + ratio * log_r, -2.0);
    for (i = 0; i < 6; i++)
      limit_metric[i] = metric_orig[i + 6 * node0] * enlarge;
    RSS(ref_matrix_intersect(&(metric_orig[6 * node1]), limit_metric, limited),
        "limit m1 with enlarged m0");
    RSS(ref_matrix_intersect(&(metric[6 * node1]), limited,
                             &(metric[6 * node1])),
        "update m1");
  }

  ref_free(metric_orig);

  ref_edge_free(ref_edge);

  RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6), "update ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_surface_gradation(REF_DBL *metric, REF_GRID ref_grid,
                                        REF_DBL r) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *metric_orig;
  REF_DBL *metric_limit;
  REF_DBL ratio, lr;
  REF_DBL l0[6], l1[6];
  REF_DBL m0[6], m1[6];
  REF_INT node, i;
  REF_INT edge, node0, node1;
  REF_BOOL have_side;

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc(metric_orig, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(metric_limit, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < 6; i++) metric_orig[i + 6 * node] = metric[i + 6 * node];
  }
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < 6; i++)
      metric_limit[i + 6 * node] = metric[i + 6 * node] * (1.0 / r / r);
  }

  each_ref_edge(ref_edge, edge) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_cell_has_side(ref_cell, node0, node1, &have_side), "side");
    if (!have_side) continue;
    RSS(ref_node_ratio(ref_grid_node(ref_grid), node0, node1, &ratio), "ratio");
    lr = pow(r, ratio);
    for (i = 0; i < 6; i++)
      l0[i] = metric_limit[i + 6 * node0] * (1.0 / lr / lr);
    for (i = 0; i < 6; i++)
      l1[i] = metric_limit[i + 6 * node1] * (1.0 / lr / lr);
    RSS(ref_matrix_intersect(&(metric_orig[6 * node0]), l1, m0), "m0");
    RSS(ref_matrix_intersect(&(metric_orig[6 * node1]), l0, m1), "m1");
    RSS(ref_matrix_intersect(m0, &(metric[6 * node0]), &(metric[6 * node0])),
        "m0");
    RSS(ref_matrix_intersect(m1, &(metric[6 * node1]), &(metric[6 * node1])),
        "m0");
  }

  ref_free(metric_limit);
  ref_free(metric_orig);

  ref_edge_free(ref_edge);

  RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), metric, 6), "update ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_sanitize(REF_GRID ref_grid) {
  if (ref_grid_twod(ref_grid)) {
    RSS(ref_metric_sanitize_twod(ref_grid), "threed");
  } else {
    RSS(ref_metric_sanitize_threed(ref_grid), "threed");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_sanitize_threed(REF_GRID ref_grid) {
  REF_DBL *metric_orig;
  REF_DBL *metric_imply;
  REF_DBL *metric;

  ref_malloc(metric_orig, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(metric_imply, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  RSS(ref_metric_from_node(metric_orig, ref_grid_node(ref_grid)), "from");

  RSS(ref_metric_imply_from(metric_imply, ref_grid), "imply");

  RSS(ref_metric_smr(metric_imply, metric_orig, metric, ref_grid), "smr");

  RSS(ref_metric_imply_non_tet(metric, ref_grid), "imply non tet");

  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "to");

  ref_free(metric);
  ref_free(metric_imply);
  ref_free(metric_orig);

  return REF_SUCCESS;
}
REF_STATUS ref_metric_sanitize_twod(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *metric_orig;
  REF_DBL *metric_imply;
  REF_DBL *metric;
  REF_INT node;

  ref_malloc(metric_orig, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(metric_imply, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);
  ref_malloc(metric, 6 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL);

  RSS(ref_metric_from_node(metric_orig, ref_grid_node(ref_grid)), "from");
  for (node = 0; node < ref_node_max(ref_node); node++) {
    metric_orig[1 + 6 * node] = 0.0;
    metric_orig[3 + 6 * node] = 1.0;
    metric_orig[4 + 6 * node] = 0.0;
  }

  RSS(ref_metric_imply_from(metric_imply, ref_grid), "imply");
  for (node = 0; node < ref_node_max(ref_node); node++) {
    metric_imply[1 + 6 * node] = 0.0;
    metric_imply[3 + 6 * node] = 1.0;
    metric_imply[4 + 6 * node] = 0.0;
  }

  RSS(ref_metric_smr(metric_imply, metric_orig, metric, ref_grid), "smr");

  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "to");

  ref_free(metric);
  ref_free(metric_imply);
  ref_free(metric_orig);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_interpolated_curvature(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *metric;
  REF_INT gradation;

  ref_malloc(metric, 6 * ref_node_max(ref_node), REF_DBL);
  RSS(ref_metric_from_curvature(metric, ref_grid), "curve");
  for (gradation = 0; gradation < 10; gradation++) {
    RSS(ref_metric_gradation(metric, ref_grid, 1.2), "grad");
  }
  RSS(ref_metric_to_node(metric, ref_node), "to node");
  ref_free(metric);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_constrain_curvature(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *curvature_metric;
  REF_DBL m[6];
  REF_INT node, im, gradation;

  if (!ref_geom_model_loaded(ref_grid_geom(ref_grid))) {
    printf("No geometry model loaded, skipping curvature constraint.\n");
    return REF_SUCCESS;
  }

  ref_malloc(curvature_metric, 6 * ref_node_max(ref_node), REF_DBL);
  RSS(ref_metric_from_curvature(curvature_metric, ref_grid), "curve");
  for (gradation = 0; gradation < 5; gradation++) {
    RSS(ref_metric_gradation(curvature_metric, ref_grid, 1.5), "grad");
  }

  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_intersect(&(curvature_metric[6 * node]),
                             ref_node_metric_ptr(ref_node, node), m),
        "intersect");
    for (im = 0; im < 6; im++) ref_node_metric(ref_node, im, node) = m[im];
  }

  ref_free(curvature_metric);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_from_curvature(REF_DBL *metric, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom, node;
  REF_DBL kr, r[3], ks, s[3], n[3];
  REF_DBL diagonal_system[12];
  REF_INT i;
  REF_DBL drad;
  REF_DBL hmax;
  REF_DBL rlimit;
  REF_DBL h, hr, hs, hn;
  REF_DBL curvature_ratio, norm_ratio;

  if (!ref_geom_model_loaded(ref_geom)) {
    printf("\nNo geometry model, did you forget to load it?\n\n");
    RSS(REF_IMPLEMENT, "...or implement non-CAD curvature estimate")
  }

  /* 1/segments per radian */
  drad = 1.0 / ref_geom_segments_per_radian_of_curvature(ref_geom);
  RSS(ref_geom_egads_diagonal(ref_geom, &hmax), "bbox diag");
  hmax *= 0.1;          /* normal spacing and max tangential spacing */
  rlimit = hmax / drad; /* h = r*drad, r = h/drad */
  curvature_ratio = 1.0 / 20.0;
  norm_ratio = 5.0;

  each_ref_node_valid_node(ref_node, node) {
    h = hmax;
    metric[0 + 6 * node] = 1.0 / (h * h);
    metric[1 + 6 * node] = 0.0;
    metric[2 + 6 * node] = 0.0;
    metric[3 + 6 * node] = 1.0 / (h * h);
    metric[4 + 6 * node] = 0.0;
    metric[5 + 6 * node] = 1.0 / (h * h);
  }

  each_ref_geom_face(ref_geom, geom) {
    RSS(ref_geom_face_curvature(ref_geom, ref_geom_id(ref_geom, geom),
                                &(ref_geom_param(ref_geom, 0, geom)), &kr, r,
                                &ks, s),
        "curve");
    kr = ABS(kr);
    ks = ABS(ks);
    kr = MAX(kr, curvature_ratio * ks);
    ks = MAX(ks, curvature_ratio * kr);
    ref_math_cross_product(r, s, n);
    node = ref_geom_node(ref_geom, geom);
    for (i = 0; i < 3; i++) ref_matrix_vec(diagonal_system, i, 0) = r[i];
    hr = hmax;
    if (1.0 / rlimit < kr) hr = drad / kr;
    ref_matrix_eig(diagonal_system, 0) = 1.0 / hr / hr;
    for (i = 0; i < 3; i++) ref_matrix_vec(diagonal_system, i, 1) = s[i];
    hs = hmax;
    if (1.0 / rlimit < ks) hs = drad / ks;
    ref_matrix_eig(diagonal_system, 1) = 1.0 / hs / hs;
    for (i = 0; i < 3; i++) ref_matrix_vec(diagonal_system, i, 2) = n[i];
    hn = hmax;
    hn = MIN(hn, norm_ratio * hr);
    hn = MIN(hn, norm_ratio * hs);
    ref_matrix_eig(diagonal_system, 2) = 1.0 / hn / hn;
    RSS(ref_matrix_form_m(diagonal_system, &(metric[6 * node])), "reform m");
  }

  return REF_SUCCESS;
}

#define sub_tet_contribution(n0, n1, n2, n3)                            \
  {                                                                     \
    tet_nodes[0] = nodes[(n0)];                                         \
    tet_nodes[1] = nodes[(n1)];                                         \
    tet_nodes[2] = nodes[(n2)];                                         \
    tet_nodes[3] = nodes[(n3)];                                         \
    RSS(ref_matrix_imply_m(m, ref_node_xyz_ptr(ref_node, tet_nodes[0]), \
                           ref_node_xyz_ptr(ref_node, tet_nodes[1]),    \
                           ref_node_xyz_ptr(ref_node, tet_nodes[2]),    \
                           ref_node_xyz_ptr(ref_node, tet_nodes[3])),   \
        "impl");                                                        \
    RSS(ref_matrix_log_m(m, log_m), "log");                             \
    RSS(ref_node_tet_vol(ref_node, tet_nodes, &tet_volume), "vol");     \
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {        \
      total_node_volume[nodes[node]] += tet_volume;                     \
      for (im = 0; im < 6; im++)                                        \
        metric[im + 6 * nodes[node]] += tet_volume * log_m[im];         \
    }                                                                   \
  }

REF_STATUS ref_metric_imply_from(REF_DBL *metric, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *total_node_volume, tet_volume;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell;
  REF_CELL ref_cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init(total_node_volume, ref_node_max(ref_node), REF_DBL, 0.0);

  for (node = 0; node < ref_node_max(ref_node); node++)
    for (im = 0; im < 6; im++) metric[im + 6 * node] = 0.0;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 1, 2, 3);
  }

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 4, 5, 3);
    sub_tet_contribution(0, 1, 5, 4);
    sub_tet_contribution(0, 1, 2, 5);
  }

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 4, 1, 2);
    sub_tet_contribution(0, 3, 4, 2);
  }

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 5, 7, 4);
    sub_tet_contribution(0, 1, 7, 5);
    sub_tet_contribution(1, 6, 7, 5);

    sub_tet_contribution(0, 7, 2, 3);
    sub_tet_contribution(0, 7, 1, 2);
    sub_tet_contribution(1, 7, 6, 2);
  }

  each_ref_node_valid_node(ref_node, node) {
    RAS(0.0 < total_node_volume[node], "zero metric contributions");
    for (im = 0; im < 6; im++) {
      if (!ref_math_divisible(metric[im + 6 * node], total_node_volume[node]))
        RSS(REF_DIV_ZERO, "zero volume");
      log_m[im] = metric[im + 6 * node] / total_node_volume[node];
    }
    RSS(ref_matrix_exp_m(log_m, m), "exp");
    for (im = 0; im < 6; im++) metric[im + 6 * node] = m[im];
    total_node_volume[node] = 0.0;
  }

  ref_free(total_node_volume);

  RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");

  return REF_SUCCESS;
}
REF_STATUS ref_metric_imply_non_tet(REF_DBL *metric, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *total_node_volume, tet_volume;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell;
  REF_CELL ref_cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init(total_node_volume, ref_node_max(ref_node), REF_DBL, 0.0);

  each_ref_node_valid_node(ref_node, node) if (ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_pyr(ref_grid)),
                                                   node)) ||
                                               ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_pri(ref_grid)),
                                                   node)) ||
                                               ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_hex(ref_grid)),
                                                   node))) for (im = 0; im < 6;
                                                                im++)
      metric[im + 6 * node] = 0.0;

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 4, 5, 3);
    sub_tet_contribution(0, 1, 5, 4);
    sub_tet_contribution(0, 1, 2, 5);
  }

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 4, 1, 2);
    sub_tet_contribution(0, 3, 4, 2);
  }

  /*
VI1 VI6 VI8 VI5  VI1 VI2 VI8 VI6  VI2 VI7 VI8 VI6
VI1 VI8 VI3 VI4  VI1 VI8 VI2 VI3  VI2 VI8 VI7 VI3
  */

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    sub_tet_contribution(0, 5, 7, 4);
    sub_tet_contribution(0, 1, 7, 5);
    sub_tet_contribution(1, 6, 7, 5);

    sub_tet_contribution(0, 7, 2, 3);
    sub_tet_contribution(0, 7, 1, 2);
    sub_tet_contribution(1, 7, 6, 2);
  }

  each_ref_node_valid_node(ref_node, node) if (ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_pyr(ref_grid)),
                                                   node)) ||
                                               ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_pri(ref_grid)),
                                                   node)) ||
                                               ref_adj_valid(ref_adj_first(
                                                   ref_cell_adj(
                                                       ref_grid_hex(ref_grid)),
                                                   node))) {
    RAS(0.0 < total_node_volume[node], "zero metric contributions");
    for (im = 0; im < 6; im++) {
      if (!ref_math_divisible(metric[im + 6 * node], total_node_volume[node]))
        RSS(REF_DIV_ZERO, "zero volume");
      log_m[im] = metric[im + 6 * node] / total_node_volume[node];
    }
    RSS(ref_matrix_exp_m(log_m, m), "exp");
    for (im = 0; im < 6; im++) metric[im + 6 * node] = m[im];
    total_node_volume[node] = 0.0;
  }

  ref_free(total_node_volume);

  RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_smr(REF_DBL *metric0, REF_DBL *metric1, REF_DBL *metric,
                          REF_GRID ref_grid) {
  REF_INT node;
  REF_DBL metric_inv[6];
  REF_DBL inv_m1_m2[9];
  REF_DBL n_values[3], n_vectors[9], inv_n_vectors[9];
  REF_DBL diagonal_system[12];
  REF_DBL h0, h1, h, hmax, hmin, h2;
  REF_DBL eig;
  REF_INT i;

  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    RSS(ref_matrix_inv_m(&(metric0[6 * node]), metric_inv), "inv");
    RSS(ref_matrix_mult_m(metric_inv, &(metric1[6 * node]), inv_m1_m2), "mult");
    RSS(ref_matrix_diag_gen(3, inv_m1_m2, n_values, n_vectors), "gen eig");
    for (i = 0; i < 3; i++) {
      h0 = ref_matrix_sqrt_vt_m_v(&(metric0[6 * node]), &(n_vectors[i * 3]));
      if (!ref_math_divisible(1.0, h0)) RSS(REF_DIV_ZERO, "inf h0");
      h0 = 1.0 / h0;
      h1 = ref_matrix_sqrt_vt_m_v(&(metric1[6 * node]), &(n_vectors[i * 3]));
      if (!ref_math_divisible(1.0, h1)) RSS(REF_DIV_ZERO, "inf h1");
      h1 = 1.0 / h1;
      hmax = 4.00 * h0;
      hmin = 0.25 * h0;
      h = MIN(hmax, MAX(hmin, h1));
      h2 = h * h;
      if (!ref_math_divisible(1.0, h2)) RSS(REF_DIV_ZERO, "zero h^2");
      eig = 1.0 / h2;
      ref_matrix_eig(diagonal_system, i) = eig;
    }
    if (REF_SUCCESS != ref_matrix_inv_gen(3, n_vectors, inv_n_vectors)) {
      printf(" unable to invert eigenvectors:\n");
      printf(" %f %f %f\n", n_vectors[0], n_vectors[1], n_vectors[2]);
      printf(" %f %f %f\n", n_vectors[3], n_vectors[4], n_vectors[5]);
      printf(" %f %f %f\n", n_vectors[6], n_vectors[7], n_vectors[8]);
      RSS(REF_FAILURE, "gen eig");
    }
    RSS(ref_matrix_transpose_gen(3, inv_n_vectors, &(diagonal_system[3])),
        "gen eig");
    RSS(ref_matrix_form_m(diagonal_system, &(metric[6 * node])), "reform m");
  }

  return REF_SUCCESS;
}

/* Alauzet and A. Loseille doi:10.1016/j.jcp.2009.09.020
 * section 2.2.4.1. A double L2-projection */
REF_STATUS ref_metric_l2_projection_grad(REF_GRID ref_grid, REF_DBL *scalar,
                                         REF_DBL *grad) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT i, node, cell, group, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL div_by_zero;
  REF_DBL cell_vol, cell_grad[3];
  REF_DBL *vol;

  ref_malloc_init(vol, ref_node_max(ref_node), REF_DBL, 0.0);

  each_ref_node_valid_node(ref_node, node) for (i = 0; i < 3; i++)
      grad[i + 3 * node] = 0.0;

  each_ref_grid_ref_cell(ref_grid, group, ref_cell)
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    switch (ref_cell_node_per(ref_cell)) {
      case 4:
        RSS(ref_node_tet_vol(ref_node, nodes, &cell_vol), "vol");
        RSS(ref_node_tet_grad(ref_node, nodes, scalar, cell_grad), "grad");
        for (cell_node = 0; cell_node < 4; cell_node++)
          for (i = 0; i < 3; i++)
            grad[i + 3 * nodes[cell_node]] += cell_vol * cell_grad[i];
        for (cell_node = 0; cell_node < 4; cell_node++)
          vol[nodes[cell_node]] += cell_vol;
        break;
      default:
        RSS(REF_IMPLEMENT, "implement cell type");
        break;
    }
  }

  div_by_zero = REF_FALSE;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_math_divisible(grad[0 + 3 * node], vol[node]) &&
        ref_math_divisible(grad[1 + 3 * node], vol[node]) &&
        ref_math_divisible(grad[2 + 3 * node], vol[node])) {
      for (i = 0; i < 3; i++) grad[i + 3 * node] /= vol[node];
    } else {
      div_by_zero = REF_TRUE;
      for (i = 0; i < 3; i++) grad[i + 3 * node] = 0.0;
    }
  }
  RSS(ref_mpi_all_or(ref_grid_mpi(ref_grid), &div_by_zero), "mpi all or");
  RSS(ref_node_ghost_dbl(ref_node, grad, 3), "update ghosts");

  ref_free(vol);

  return (div_by_zero ? REF_DIV_ZERO : REF_SUCCESS);
}

REF_STATUS ref_metric_l2_projection_hessian(REF_GRID ref_grid, REF_DBL *scalar,
                                            REF_DBL *hessian) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i, node;
  REF_DBL *grad, *dsdx, *gradx, *grady, *gradz;
  REF_DBL diag_system[12];

  ref_malloc_init(grad, 3 * ref_node_max(ref_node), REF_DBL, 0.0);
  ref_malloc_init(dsdx, ref_node_max(ref_node), REF_DBL, 0.0);
  ref_malloc_init(gradx, 3 * ref_node_max(ref_node), REF_DBL, 0.0);
  ref_malloc_init(grady, 3 * ref_node_max(ref_node), REF_DBL, 0.0);
  ref_malloc_init(gradz, 3 * ref_node_max(ref_node), REF_DBL, 0.0);

  RSS(ref_metric_l2_projection_grad(ref_grid, scalar, grad), "l2 grad");

  i = 0;
  each_ref_node_valid_node(ref_node, node) dsdx[node] = grad[i + 3 * node];
  RSS(ref_metric_l2_projection_grad(ref_grid, dsdx, gradx), "gradx");

  i = 1;
  each_ref_node_valid_node(ref_node, node) dsdx[node] = grad[i + 3 * node];
  RSS(ref_metric_l2_projection_grad(ref_grid, dsdx, grady), "grady");

  i = 2;
  each_ref_node_valid_node(ref_node, node) dsdx[node] = grad[i + 3 * node];
  RSS(ref_metric_l2_projection_grad(ref_grid, dsdx, gradz), "gradz");

  /* average off-diagonals */
  each_ref_node_valid_node(ref_node, node) {
    hessian[0 + 6 * node] = gradx[0 + 3 * node];
    hessian[1 + 6 * node] = 0.5 * (gradx[1 + 3 * node] + grady[0 + 3 * node]);
    hessian[2 + 6 * node] = 0.5 * (gradx[2 + 3 * node] + gradz[0 + 3 * node]);
    hessian[3 + 6 * node] = grady[1 + 3 * node];
    hessian[4 + 6 * node] = 0.5 * (grady[2 + 3 * node] + gradz[1 + 3 * node]);
    hessian[5 + 6 * node] = gradz[2 + 3 * node];
  }

  /* positive eignevalues to make symmetric positive definite */
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_diag_m(&(hessian[6 * node]), diag_system), "eigen decomp");
    ref_matrix_eig(diag_system, 0) = ABS(ref_matrix_eig(diag_system, 0));
    ref_matrix_eig(diag_system, 1) = ABS(ref_matrix_eig(diag_system, 1));
    ref_matrix_eig(diag_system, 2) = ABS(ref_matrix_eig(diag_system, 2));
    RSS(ref_matrix_form_m(diag_system, &(hessian[6 * node])), "re-form hess");
  }

  ref_free(gradz);
  ref_free(grady);
  ref_free(gradx);
  ref_free(dsdx);
  ref_free(grad);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_kexact_hessian(REF_GRID ref_grid, REF_DBL *scalar,
                                     REF_DBL *hessian) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT node0, node1, node2, i1, i2, im, i, j;
  REF_INT nnode1, nnode2;
  REF_INT node_list1[REF_METRIC_MAX_DEGREE], node_list2[REF_METRIC_MAX_DEGREE],
      max_node = REF_METRIC_MAX_DEGREE;
  REF_DICT ref_dict;
  REF_DBL geom[9], ab[90];
  REF_DBL dx, dy, dz, dq;
  REF_DBL *a, *q, *r;
  REF_INT m, n;
  each_ref_node_valid_node(ref_node, node0) {
    /* use ref_dict to get a unique list of halo(2) nodes */
    RSS(ref_dict_create(&ref_dict), "create ref_dict");
    RXS(ref_cell_node_list_around(ref_cell, node0, max_node, &nnode1,
                                  node_list1),
        REF_INCREASE_LIMIT, "first halo of nodes");
    for (i1 = 0; i1 < nnode1; i1++) {
      node1 = node_list1[i1];
      RSS(ref_dict_store(ref_dict, node1, REF_EMPTY), "store node1");
      RXS(ref_cell_node_list_around(ref_cell, node1, max_node, &nnode2,
                                    node_list2),
          REF_INCREASE_LIMIT, "halo of halo of nodes");
      for (i2 = 0; i2 < nnode2; i2++) {
        node2 = node_list2[i2];
        if (node0 == node2) continue; /* skip self */
        RSS(ref_dict_store(ref_dict, node2, REF_EMPTY), "store node2");
      }
    }
    /* solve A with QR factorization size m x n */
    m = ref_dict_n(ref_dict);
    n = 9;
    ref_malloc(a, m * n, REF_DBL);
    ref_malloc(q, m * n, REF_DBL);
    ref_malloc(r, n * n, REF_DBL);
    i = 0;
    each_ref_dict_key(ref_dict, i2, node2) {
      dx = ref_node_xyz(ref_node, 0, node2) - ref_node_xyz(ref_node, 0, node0);
      dy = ref_node_xyz(ref_node, 1, node2) - ref_node_xyz(ref_node, 1, node0);
      dz = ref_node_xyz(ref_node, 2, node2) - ref_node_xyz(ref_node, 2, node0);
      geom[0] = 0.5 * dx * dx;
      geom[1] = dx * dy;
      geom[2] = dx * dz;
      geom[3] = 0.5 * dy * dy;
      geom[4] = dy * dz;
      geom[5] = 0.5 * dz * dz;
      geom[6] = dx;
      geom[7] = dy;
      geom[8] = dz;
      for (j = 0; j < n; j++) {
        a[i + m * j] = geom[j];
      }
      i++;
    }
    RSS(ref_matrix_qr(m, n, a, q, r), "kexact lsq hess qr");
    for (i = 0; i < 90; i++) ab[i] = 0.0;
    for (i = 0; i < 9; i++) {
      for (j = 0; j < 9; j++) {
        ab[i + 9 * j] += r[i + 9 * j];
      }
    }
    i = 0;
    each_ref_dict_key(ref_dict, i2, node2) {
      dq = scalar[node2] - scalar[node0];
      for (j = 0; j < 9; j++) {
        ab[j + 9 * 9] += q[i + m * j] * dq;
      }
      i++;
    }
    RSS(ref_matrix_solve_ab(9, 10, ab), "solve rx=qtb");
    j = 9;
    for (im = 0; im < 6; im++) {
      hessian[im + 6 * node0] = ab[im + 9 * j];
    }
    ref_free(r);
    ref_free(q);
    ref_free(a);
    RSS(ref_dict_free(ref_dict), "free ref_dict");
  }

  /* positive eignevalues to make symmetric positive definite */
  each_ref_node_valid_node(ref_node, node0) {
    REF_DBL diag_system[12];
    RSS(ref_matrix_diag_m(&(hessian[6 * node0]), diag_system), "eigen decomp");
    ref_matrix_eig(diag_system, 0) = ABS(ref_matrix_eig(diag_system, 0));
    ref_matrix_eig(diag_system, 1) = ABS(ref_matrix_eig(diag_system, 1));
    ref_matrix_eig(diag_system, 2) = ABS(ref_matrix_eig(diag_system, 2));
    RSS(ref_matrix_form_m(diag_system, &(hessian[6 * node0])), "re-form hess");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_extrapolate_boundary(REF_DBL *metric, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tris = ref_grid_tri(ref_grid);
  REF_CELL tets = ref_grid_tet(ref_grid);
  REF_INT node;
  REF_INT max_node = REF_METRIC_MAX_DEGREE, nnode;
  REF_INT node_list[REF_METRIC_MAX_DEGREE];
  REF_INT i, neighbor, nint;

  if (ref_grid_twod(ref_grid)) RSS(REF_IMPLEMENT, "2D not implmented");

  /* each boundary node */
  each_ref_node_valid_node(ref_node,
                           node) if (!ref_cell_node_empty(tris, node)) {
    RXS(ref_cell_node_list_around(tets, node, max_node, &nnode, node_list),
        REF_INCREASE_LIMIT, "unable to build neighbor list ");
    nint = 0;
    for (neighbor = 0; neighbor < nnode; neighbor++)
      if (ref_cell_node_empty(tris, node_list[neighbor])) nint++;
    if (0 < nint) {
      for (i = 0; i < 6; i++) metric[i + 6 * node] = 0.0;
      for (neighbor = 0; neighbor < nnode; neighbor++)
        if (ref_cell_node_empty(tris, node_list[neighbor])) {
          /* use Euclidean average, these are derivatives */
          for (i = 0; i < 6; i++)
            metric[i + 6 * node] += metric[i + 6 * node_list[neighbor]];
        }
      for (i = 0; i < 6; i++) metric[i + 6 * node] /= (REF_DBL)nint;
    }
  }

  RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_extrapolate_boundary_multipass(REF_DBL *metric,
                                                     REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tris = ref_grid_tri(ref_grid);
  REF_CELL tets = ref_grid_tet(ref_grid);
  REF_INT node;
  REF_INT max_node = REF_METRIC_MAX_DEGREE, nnode;
  REF_INT node_list[REF_METRIC_MAX_DEGREE];
  REF_INT i, neighbor, nint;
  REF_BOOL *needs_donor;
  REF_INT pass, remain;

  if (ref_grid_twod(ref_grid)) RSS(REF_IMPLEMENT, "2D not implmented");

  ref_malloc_init(needs_donor, ref_node_max(ref_node), REF_BOOL, REF_FALSE);

  /* each boundary node */
  each_ref_node_valid_node(ref_node,
                           node) if (!ref_cell_node_empty(tris, node)) {
    needs_donor[node] = REF_TRUE;
  }

  RSS(ref_node_ghost_int(ref_node, needs_donor), "update ghosts");

  for (pass = 0; pass < 10; pass++) {
    each_ref_node_valid_node(
        ref_node,
        node) if (ref_node_owned(ref_node, node) && needs_donor[node]) {
      RXS(ref_cell_node_list_around(tets, node, max_node, &nnode, node_list),
          REF_INCREASE_LIMIT, "unable to build neighbor list ");
      nint = 0;
      for (neighbor = 0; neighbor < nnode; neighbor++)
        if (!needs_donor[node_list[neighbor]]) nint++;
      if (0 < nint) {
        for (i = 0; i < 6; i++) metric[i + 6 * node] = 0.0;
        for (neighbor = 0; neighbor < nnode; neighbor++)
          if (!needs_donor[node_list[neighbor]]) {
            for (i = 0; i < 6; i++)
              metric[i + 6 * node] += metric[i + 6 * node_list[neighbor]];
          }
        for (i = 0; i < 6; i++) metric[i + 6 * node] /= (REF_DBL)nint;
        needs_donor[node] = REF_FALSE;
      }
    }

    RSS(ref_node_ghost_int(ref_node, needs_donor), "update ghosts");
    RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");

    remain = 0;
    each_ref_node_valid_node(ref_node, node) {
      if (ref_node_owned(ref_node, node) && needs_donor[node]) {
        remain++;
      }
    }
    RSS(ref_mpi_allsum(ref_grid_mpi(ref_grid), &remain, 1, REF_INT_TYPE),
        "sum updates");

    if (0 == remain) break;
  }

  ref_free(needs_donor);

  REIS(0, remain, "untouched boundary nodes remain");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_complexity(REF_DBL *metric, REF_GRID ref_grid,
                                 REF_DBL *complexity) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell_node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL volume, det;
  if (ref_grid_twod(ref_grid)) ref_cell = ref_grid_tri(ref_grid);
  *complexity = 0.0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_grid_twod(ref_grid)) {
      RSS(ref_node_tri_area(ref_node, nodes, &volume), "area");
    } else {
      RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");
    }
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (ref_node_owned(ref_node, nodes[cell_node])) {
        RSS(ref_matrix_det_m(&(metric[6 * nodes[cell_node]]), &det), "det");
        if (ABS(det) > 1.0e-15) {
          (*complexity) +=
              sqrt(det) * volume / ((REF_DBL)ref_cell_node_per(ref_cell));
        }
      }
    }
  }
  RSS(ref_mpi_allsum(ref_grid_mpi(ref_grid), complexity, 1, REF_DBL_TYPE),
      "dbl sum");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_lp(REF_DBL *metric, REF_GRID ref_grid, REF_DBL *scalar,
                         REF_METRIC_RECONSTRUCTION reconstruction,
                         REF_INT p_norm, REF_DBL gradation,
                         REF_DBL target_complexity) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i, node;
  REF_INT dimension = 3;
  REF_INT relaxations;
  REF_DBL det, exponent;
  REF_DBL current_complexity;
  if (ref_grid_twod(ref_grid)) RSS(REF_IMPLEMENT, "2D not implmented");
  switch (reconstruction) {
    case REF_METRIC_L2PROJECTION:
      RSS(ref_metric_l2_projection_hessian(ref_grid, scalar, metric), "l2");
      RSS(ref_metric_extrapolate_boundary_multipass(metric, ref_grid),
          "bound extrap");
      break;
    case REF_METRIC_KEXACT:
      RSS(ref_metric_kexact_hessian(ref_grid, scalar, metric), "k-exact");
      break;
    default:
      THROW("reconstruction not available");
  }

  RSS(ref_metric_roundoff_limit(metric, ref_grid),
      "floor metric eignvalues based on grid size and solution jitter");

  /* local scaling */
  exponent = -1.0 / ((REF_DBL)(2 * p_norm + dimension));
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_det_m(&(metric[6 * node]), &det), "det_m local hess scale");
    if (ABS(det) > 1.0e-15) {
      for (i = 0; i < 6; i++) metric[i + 6 * node] *= pow(det, exponent);
    }
  }
  /* global scaling and gradation limiting */
  for (relaxations = 0; relaxations < 10; relaxations++) {
    RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (!ref_math_divisible(target_complexity, current_complexity)) {
      return REF_DIV_ZERO;
    }
    each_ref_node_valid_node(ref_node, node) for (i = 0; i < 6; i++)
        metric[i + 6 * node] *=
        pow(target_complexity / current_complexity, 2.0 / 3.0);
    RSS(ref_metric_gradation(metric, ref_grid, gradation), "gradation");
  }
  RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
  if (!ref_math_divisible(target_complexity, current_complexity)) {
    return REF_DIV_ZERO;
  }
  each_ref_node_valid_node(ref_node, node) for (i = 0; i < 6; i++)
      metric[i + 6 * node] *=
      pow(target_complexity / current_complexity, 2.0 / 3.0);
  return REF_SUCCESS;
}

REF_STATUS ref_metric_set_zero_det(REF_DBL *metric_with_zeros,
                                   REF_DBL *metric_donor, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i, node;
  REF_DBL det;

  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_det_m(&(metric_with_zeros[6 * node]), &det),
        "det_m to check validity");
    if (ABS(det) < 1.0e-15) {
      for (i = 0; i < 6; i++) {
        metric_with_zeros[i + 6 * node] = metric_donor[i + 6 * node];
      }
    }
  }

  return REF_SUCCESS;
}
REF_STATUS ref_metric_roundoff_limit(REF_DBL *metric, REF_GRID ref_grid) {
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i, node;
  REF_DBL radius, dist;
  REF_DBL round_off_jitter = 1.0e-12;
  REF_DBL eig_floor;
  REF_INT nnode, node_list[REF_METRIC_MAX_DEGREE],
      max_node = REF_METRIC_MAX_DEGREE;
  REF_DBL diag_system[12];

  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_cell_node_list_around(ref_cell, node, max_node, &nnode, node_list),
        "first halo of nodes");
    radius = 0.0;
    for (i = 0; i < nnode; i++) {
      dist = sqrt(pow(ref_node_xyz(ref_node, 0, node_list[i]) -
                          ref_node_xyz(ref_node, 0, node),
                      2) +
                  pow(ref_node_xyz(ref_node, 1, node_list[i]) -
                          ref_node_xyz(ref_node, 1, node),
                      2) +
                  pow(ref_node_xyz(ref_node, 2, node_list[i]) -
                          ref_node_xyz(ref_node, 2, node),
                      2));
      if (i == 0) radius = dist;
      radius = MIN(radius, dist);
    }
    /* 2nd order central finite difference */
    eig_floor = 4 * round_off_jitter / radius / radius;

    RSS(ref_matrix_diag_m(&(metric[6 * node]), diag_system), "eigen decomp");
    ref_matrix_eig(diag_system, 0) =
        MAX(ref_matrix_eig(diag_system, 0), eig_floor);
    ref_matrix_eig(diag_system, 1) =
        MAX(ref_matrix_eig(diag_system, 1), eig_floor);
    ref_matrix_eig(diag_system, 2) =
        MAX(ref_matrix_eig(diag_system, 2), eig_floor);
    RSS(ref_matrix_form_m(diag_system, &(metric[6 * node])), "re-form hess");
  }

  RSS(ref_node_ghost_dbl(ref_node, metric, 6), "update ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_opt_goal(REF_DBL *metric, REF_GRID ref_grid,
                               REF_INT nequations, REF_DBL *solution,
                               REF_DBL target_complexity) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i, node;
  REF_INT dimension = 3;
  REF_INT p_norm = 1;
  REF_DBL det, exponent;
  REF_DBL current_complexity;
  REF_INT ldim;
  REF_INT var, dir;
  if (ref_grid_twod(ref_grid)) RSS(REF_IMPLEMENT, "2D not implmented");

  ldim = 4 * nequations;

  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < 6; i++) metric[i + 6 * node] = 0.0;
  }

  for (var = 0; var < nequations; var++) {
    REF_DBL *lam, *grad_lam, *flux, *hess_flux;
    ref_malloc_init(lam, ref_node_max(ref_node), REF_DBL, 0.0);
    ref_malloc_init(grad_lam, 3 * ref_node_max(ref_node), REF_DBL, 0.0);
    ref_malloc_init(flux, ref_node_max(ref_node), REF_DBL, 0.0);
    ref_malloc_init(hess_flux, 6 * ref_node_max(ref_node), REF_DBL, 0.0);
    each_ref_node_valid_node(ref_node, node) {
      lam[node] = solution[var + 5 + ldim * node];
    }
    RSS(ref_metric_l2_projection_grad(ref_grid, lam, grad_lam), "grad_lam");

    for (dir = 0; dir < 3; dir++) {
      each_ref_node_valid_node(ref_node, node) {
        flux[node] = solution[var + nequations * dir + ldim * node];
      }
      RSS(ref_metric_l2_projection_hessian(ref_grid, flux, hess_flux), "l2");
      each_ref_node_valid_node(ref_node, node) {
        for (i = 0; i < 6; i++)
          metric[i + 6 * node] +=
              ABS(grad_lam[dir + 3 * node]) * hess_flux[i + 6 * node];
      }
    }
    ref_free(hess_flux);
    ref_free(flux);
    ref_free(grad_lam);
    ref_free(lam);
  }

  /* local scaling */
  exponent = -1.0 / ((REF_DBL)(2 * p_norm + dimension));
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_det_m(&(metric[6 * node]), &det), "det_m local hess scale");
    if (ABS(det) > 1.0e-15) {
      for (i = 0; i < 6; i++) metric[i + 6 * node] *= pow(det, exponent);
    }
  }

  RSS(ref_metric_complexity(metric, ref_grid, &current_complexity), "cmp");
  if (!ref_math_divisible(target_complexity, current_complexity)) {
    return REF_DIV_ZERO;
  }
  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < 6; i++)
      metric[i + 6 * node] *=
          pow(target_complexity / current_complexity, 2.0 / 3.0);
  }

  return REF_SUCCESS;
}
