
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

#include "ref_adapt.h"
#include "ref_edge.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"

#include "ref_cavity.h"
#include "ref_collapse.h"
#include "ref_smooth.h"
#include "ref_split.h"
#include "ref_swap.h"

#include "ref_matrix.h"
#include "ref_node.h"

#include "ref_gather.h"
#include "ref_histogram.h"
#include "ref_metric.h"

REF_STATUS ref_adapt_create(REF_ADAPT *ref_adapt_ptr) {
  REF_ADAPT ref_adapt;

  ref_malloc(*ref_adapt_ptr, 1, REF_ADAPT_STRUCT);

  ref_adapt = *ref_adapt_ptr;

  ref_adapt->split_ratio_growth = REF_FALSE;
  ref_adapt->split_ratio = sqrt(2.0);
  ref_adapt->split_quality_absolute = 1.0e-3;
  ref_adapt->split_quality_relative = 0.1;

  ref_adapt->collapse_ratio = 1.0 / sqrt(2.0);
  ref_adapt->collapse_quality_absolute = 1.0e-3;

  ref_adapt->smooth_min_quality = 1.0e-3;

  ref_adapt->swap_max_degree = 10000;
  ref_adapt->swap_min_quality = 0.8;

  ref_adapt->post_min_normdev = 0.0;
  ref_adapt->post_min_ratio = 1.0e-3;
  ref_adapt->post_max_ratio = 3.0;

  ref_adapt->last_min_ratio = 0.5e-3;
  ref_adapt->last_max_ratio = 6.0;

  ref_adapt->instrument = REF_TRUE;
  ref_adapt->watch_param = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_adapt_deep_copy(REF_ADAPT *ref_adapt_ptr, REF_ADAPT original) {
  REF_ADAPT ref_adapt;

  ref_malloc(*ref_adapt_ptr, 1, REF_ADAPT_STRUCT);

  ref_adapt = *ref_adapt_ptr;

  ref_adapt->split_ratio_growth = original->split_ratio_growth;
  ref_adapt->split_ratio = original->split_ratio;
  ref_adapt->split_quality_absolute = original->split_quality_absolute;
  ref_adapt->split_quality_relative = original->split_quality_relative;

  ref_adapt->collapse_ratio = original->collapse_ratio;
  ref_adapt->collapse_quality_absolute = original->collapse_quality_absolute;

  ref_adapt->smooth_min_quality = original->smooth_min_quality;

  ref_adapt->swap_max_degree = original->swap_max_degree;
  ref_adapt->swap_min_quality = original->swap_min_quality;

  ref_adapt->post_min_normdev = original->post_min_normdev;
  ref_adapt->post_min_ratio = original->post_min_ratio;
  ref_adapt->post_max_ratio = original->post_max_ratio;

  ref_adapt->last_min_ratio = original->last_min_ratio;
  ref_adapt->last_max_ratio = original->last_max_ratio;

  ref_adapt->instrument = original->instrument;
  ref_adapt->watch_param = original->watch_param;

  return REF_SUCCESS;
}

REF_STATUS ref_adapt_free(REF_ADAPT ref_adapt) {
  ref_free(ref_adapt);

  return REF_SUCCESS;
}

static REF_STATUS ref_adapt_parameter(REF_GRID ref_grid, REF_BOOL *all_done) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_ADAPT ref_adapt = ref_grid->adapt;
  REF_CELL ref_cell;
  REF_INT cell;
  REF_LONG ncell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL det, max_det, complexity, min_metric_vol;
  REF_DBL quality, min_quality;
  REF_DBL normdev, min_normdev;
  REF_DBL volume, min_volume, max_volume;
  REF_BOOL active_twod;
  REF_DBL target_quality, target_normdev;
  REF_INT cell_node;
  REF_INT node, nnode;
  REF_DBL nodes_per_complexity;
  REF_INT degree, max_degree;
  REF_DBL ratio, min_ratio, max_ratio;
  REF_INT edge, part;
  REF_INT age, max_age;
  REF_BOOL active;
  REF_EDGE ref_edge;
  REF_DBL m[6];

  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  min_quality = 1.0;
  min_volume = REF_DBL_MAX;
  max_volume = REF_DBL_MIN;
  max_det = -1.0;
  complexity = 0.0;
  ncell = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_grid_twod(ref_grid)) {
      RSS(ref_node_node_twod(ref_grid_node(ref_grid), nodes[0], &active_twod),
          "active twod tri");
      if (!active_twod) continue;
      RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
      RSS(ref_node_tri_area(ref_grid_node(ref_grid), nodes, &volume), "vol");
    } else if (ref_grid_surf(ref_grid)) {
      REF_INT id;
      REF_DBL area_sign, uv_area;
      RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
      id = nodes[ref_cell_node_per(ref_cell)];
      RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign), "a sign");
      RSS(ref_geom_uv_area(ref_grid_geom(ref_grid), nodes, &uv_area),
          "uv area");
      volume = area_sign * uv_area;
    } else {
      RSS(ref_node_tet_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
      RSS(ref_node_tet_vol(ref_grid_node(ref_grid), nodes, &volume), "vol");
    }
    min_quality = MIN(min_quality, quality);
    min_volume = MIN(min_volume, volume);
    max_volume = MAX(max_volume, volume);

    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (ref_node_owned(ref_node, nodes[cell_node])) {
        RSS(ref_node_metric_get(ref_node, nodes[cell_node], m), "get");
        RSS(ref_matrix_det_m(m, &det), "det");
        max_det = MAX(max_det, det);
        if (ref_grid_surf(ref_grid)) {
          REF_DBL area, normal[3], normal_projection;
          RSS(ref_node_tri_area(ref_grid_node(ref_grid), nodes, &area), "vol");
          RSS(ref_node_tri_normal(ref_grid_node(ref_grid), nodes, normal),
              "norm");
          RSS(ref_math_normalize(normal), "normalize");
          normal_projection = ref_matrix_vt_m_v(m, normal);
          if (ref_math_divisible(det, normal_projection)) {
            if (det / normal_projection > 0.0) {
              complexity += sqrt(det / normal_projection) * area /
                            ((REF_DBL)ref_cell_node_per(ref_cell));
            }
          }
        } else {
          if (det > 0.0) {
            complexity +=
                sqrt(det) * volume / ((REF_DBL)ref_cell_node_per(ref_cell));
          }
        }
      }
    }
    RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "owner");
    if (part == ref_mpi_rank(ref_mpi)) ncell++;
  }
  quality = min_quality;
  RSS(ref_mpi_min(ref_mpi, &quality, &min_quality, REF_DBL_TYPE), "mpi min");
  RSS(ref_mpi_bcast(ref_mpi, &quality, 1, REF_DBL_TYPE), "mbast");
  volume = min_volume;
  RSS(ref_mpi_min(ref_mpi, &volume, &min_volume, REF_DBL_TYPE), "mpi min");
  RSS(ref_mpi_bcast(ref_mpi, &min_volume, 1, REF_DBL_TYPE), "bcast");
  volume = max_volume;
  RSS(ref_mpi_max(ref_mpi, &volume, &max_volume, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_volume, 1, REF_DBL_TYPE), "bcast");
  det = max_det;
  RSS(ref_mpi_max(ref_mpi, &det, &max_det, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_det, 1, REF_DBL_TYPE), "bcast");
  RAS(ref_math_divisible(1.0, sqrt(max_det)), "can not invert sqrt(max_det)");
  min_metric_vol = 1.0 / sqrt(max_det);

  RSS(ref_mpi_allsum(ref_mpi, &complexity, 1, REF_DBL_TYPE), "dbl sum");
  RSS(ref_mpi_allsum(ref_mpi, &ncell, 1, REF_LONG_TYPE), "cell int sum");

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      nnode++;
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, &nnode, 1, REF_INT_TYPE), "int sum");
  if (ref_grid_twod(ref_grid)) nnode = nnode / 2;

  nodes_per_complexity = (REF_DBL)nnode / complexity;

  max_degree = 0;
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_adj_degree(ref_cell_adj(ref_cell), node, &degree), "cell degree");
    max_degree = MAX(max_degree, degree);
  }
  degree = max_degree;
  RSS(ref_mpi_max(ref_mpi, &degree, &max_degree, REF_INT_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_degree, 1, REF_INT_TYPE), "min");

  max_age = 0;
  each_ref_node_valid_node(ref_node, node) {
    max_age = MAX(max_age, ref_node_age(ref_node, node));
  }
  age = max_age;
  RSS(ref_mpi_max(ref_mpi, &age, &max_age, REF_INT_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_age, 1, REF_INT_TYPE), "min");

  min_normdev = 2.0;
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid))) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "norm dev");
      min_normdev = MIN(min_normdev, normdev);
    }
  }
  normdev = min_normdev;
  RSS(ref_mpi_min(ref_mpi, &normdev, &min_normdev, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &min_normdev, 1, REF_DBL_TYPE), "min");

  min_ratio = REF_DBL_MAX;
  max_ratio = REF_DBL_MIN;
  RSS(ref_edge_create(&ref_edge, ref_grid), "make edges");
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    RSS(ref_node_edge_twod(ref_grid_node(ref_grid),
                           ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), &active),
        "twod edge");
    active = (active || !ref_grid_twod(ref_grid));
    if (part == ref_mpi_rank(ref_grid_mpi(ref_grid)) && active) {
      RSS(ref_node_ratio(ref_grid_node(ref_grid),
                         ref_edge_e2n(ref_edge, 0, edge),
                         ref_edge_e2n(ref_edge, 1, edge), &ratio),
          "rat");
      min_ratio = MIN(min_ratio, ratio);
      max_ratio = MAX(max_ratio, ratio);
    }
  }
  RSS(ref_edge_free(ref_edge), "free edge");
  ratio = min_ratio;
  RSS(ref_mpi_min(ref_mpi, &ratio, &min_ratio, REF_DBL_TYPE), "mpi min");
  RSS(ref_mpi_bcast(ref_mpi, &min_ratio, 1, REF_DBL_TYPE), "min");
  ratio = max_ratio;
  RSS(ref_mpi_max(ref_mpi, &ratio, &max_ratio, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_ratio, 1, REF_DBL_TYPE), "max");

  target_normdev = MAX(MIN(0.1, min_normdev), 1.0e-3);
  ref_adapt->post_min_normdev = target_normdev;

  target_quality = MAX(MIN(0.1, min_quality), 1.0e-3);
  ref_adapt->collapse_quality_absolute = target_quality;
  ref_adapt->smooth_min_quality = target_quality;

  ref_node->min_volume = MIN(1.0e-15, 0.01 * min_metric_vol);

  /* allow edge growth when interpolating metric continuously */
  ref_adapt->split_ratio_growth = REF_FALSE;
  if (NULL != ref_grid_interp(ref_grid)) {
    ref_adapt->split_ratio_growth =
        ref_interp_continuously(ref_grid_interp(ref_grid));
  }

  /* bound ratio to current range */
  ref_adapt->post_min_ratio = MIN(min_ratio, ref_adapt->collapse_ratio);
  ref_adapt->post_max_ratio = MAX(max_ratio, ref_adapt->split_ratio);

  if (ref_adapt->post_max_ratio > 4.0 && ref_adapt->post_min_ratio > 0.4) {
    ref_adapt->post_min_ratio =
        (4.0 / ref_adapt->post_max_ratio) * ref_adapt->post_min_ratio;
  }

  ref_adapt->split_ratio = sqrt(2.0);
  if (nodes_per_complexity > 3.0)
    ref_adapt->split_ratio = 0.5 * (sqrt(2.0) + max_ratio);

  if (ABS(ref_adapt->last_min_ratio - ref_adapt->post_min_ratio) <
          1e-2 * ref_adapt->post_min_ratio &&
      ABS(ref_adapt->last_max_ratio - ref_adapt->post_max_ratio) <
          1e-2 * ref_adapt->post_max_ratio &&
      (max_age < 50 ||
       (ref_adapt->post_min_ratio > 0.1 && ref_adapt->post_max_ratio < 3.0)) &&
      1.5 > ref_adapt->split_ratio) {
    *all_done = REF_TRUE;
    if (ref_grid_once(ref_grid)) {
      printf("termination recommended\n");
    }
  } else {
    *all_done = REF_FALSE;
  }
  RSS(ref_mpi_bcast(ref_mpi, all_done, 1, REF_INT_TYPE), "done");

  ref_adapt->last_min_ratio = ref_adapt->post_min_ratio;
  ref_adapt->last_max_ratio = ref_adapt->post_max_ratio;

  if (ref_grid_once(ref_grid)) {
    printf("limit quality %6.4f normdev %6.4f ratio %6.4f %6.2f split %6.2f\n",
           target_quality, ref_adapt->post_min_normdev,
           ref_adapt->post_min_ratio, ref_adapt->post_max_ratio,
           ref_adapt->split_ratio);
    printf("max degree %d max age %d normdev %7.4f\n", max_degree, max_age,
           min_normdev);
    printf("nnode %10d ncell %12ld complexity %12.1f ratio %5.2f\n", nnode,
           ncell, complexity, nodes_per_complexity);
    printf("volume range %e %e metric %e floor %e\n", max_volume, min_volume,
           min_metric_vol, ref_node->min_volume);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_adapt_tattle(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, min_quality;
  REF_DBL normdev, min_normdev;
  REF_BOOL active_twod;
  REF_INT node, nnode;
  REF_DBL ratio, min_ratio, max_ratio;
  REF_INT edge, part;
  REF_BOOL active;
  REF_EDGE ref_edge;
  char is_ok = ' ';
  char not_ok = '*';
  char quality_met, short_met, long_met, normdev_met;

  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  min_quality = 1.0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_grid_twod(ref_grid)) {
      RSS(ref_node_node_twod(ref_grid_node(ref_grid), nodes[0], &active_twod),
          "active twod tri");
      if (!active_twod) continue;
      RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
    } else if (ref_grid_surf(ref_grid)) {
      RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
    } else {
      RSS(ref_node_tet_quality(ref_grid_node(ref_grid), nodes, &quality),
          "qual");
    }
    min_quality = MIN(min_quality, quality);
  }
  quality = min_quality;
  RSS(ref_mpi_min(ref_mpi, &quality, &min_quality, REF_DBL_TYPE), "min");
  RSS(ref_mpi_bcast(ref_mpi, &quality, 1, REF_DBL_TYPE), "min");

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      nnode++;
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, &nnode, 1, REF_INT_TYPE), "int sum");
  if (ref_grid_twod(ref_grid)) nnode = nnode / 2;

  min_normdev = 2.0;
  if (ref_geom_model_loaded(ref_grid_geom(ref_grid))) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "norm dev");
      min_normdev = MIN(min_normdev, normdev);
    }
  }
  normdev = min_normdev;
  RSS(ref_mpi_min(ref_mpi, &normdev, &min_normdev, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &min_normdev, 1, REF_DBL_TYPE), "min");

  min_ratio = REF_DBL_MAX;
  max_ratio = REF_DBL_MIN;
  RSS(ref_edge_create(&ref_edge, ref_grid), "make edges");
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    RSS(ref_node_edge_twod(ref_grid_node(ref_grid),
                           ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), &active),
        "twod edge");
    active = (active || !ref_grid_twod(ref_grid));
    if (part == ref_mpi_rank(ref_grid_mpi(ref_grid)) && active) {
      RSS(ref_node_ratio(ref_grid_node(ref_grid),
                         ref_edge_e2n(ref_edge, 0, edge),
                         ref_edge_e2n(ref_edge, 1, edge), &ratio),
          "rat");
      min_ratio = MIN(min_ratio, ratio);
      max_ratio = MAX(max_ratio, ratio);
    }
  }
  RSS(ref_edge_free(ref_edge), "free edge");
  ratio = min_ratio;
  RSS(ref_mpi_min(ref_mpi, &ratio, &min_ratio, REF_DBL_TYPE), "mpi min");
  RSS(ref_mpi_bcast(ref_mpi, &min_ratio, 1, REF_DBL_TYPE), "min");
  ratio = max_ratio;
  RSS(ref_mpi_max(ref_mpi, &ratio, &max_ratio, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_ratio, 1, REF_DBL_TYPE), "max");

  if (ref_grid_once(ref_grid)) {
    quality_met = is_ok;
    short_met = is_ok;
    long_met = is_ok;
    normdev_met = is_ok;
    if (min_quality < ref_grid_adapt(ref_grid, smooth_min_quality))
      quality_met = not_ok;
    if (min_ratio < ref_grid_adapt(ref_grid, post_min_ratio))
      short_met = not_ok;
    if (max_ratio > ref_grid_adapt(ref_grid, post_max_ratio)) long_met = not_ok;
    if (min_normdev < ref_grid_adapt(ref_grid, post_min_normdev))
      normdev_met = not_ok;

    printf(
        "quality %c %6.4f ratio %c %6.4f %6.2f %c normdev %6.4f %c nnode %d\n",
        quality_met, min_quality, short_met, min_ratio, max_ratio, long_met,
        min_normdev, normdev_met, nnode);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_adapt_threed_swap(REF_GRID ref_grid) {
  REF_INT pass;
  RSS(ref_cavity_pass(ref_grid), "cavity pass");
  if (ref_grid_surf(ref_grid)) {
    for (pass = 0; pass < 3; pass++) {
      RSS(ref_swap_surf_pass(ref_grid), "swap pass");
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_adapt_threed_pass(REF_GRID ref_grid, REF_BOOL *all_done) {
  REF_INT ngeom;
  REF_BOOL all_done0, all_done1;

  RSS(ref_adapt_parameter(ref_grid, &all_done0), "param");

  RSS(ref_gather_ngeom(ref_grid_node(ref_grid), ref_grid_geom(ref_grid),
                       REF_GEOM_FACE, &ngeom),
      "count ngeom");
  if (ngeom > 0) RSS(ref_geom_verify_topo(ref_grid), "adapt preflight check");

  ref_gather_blocking_frame(ref_grid, "threed pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt start");

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  RSS(ref_collapse_pass(ref_grid), "col pass");
  ref_gather_blocking_frame(ref_grid, "collapse");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt col");

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  ref_grid_adapt(ref_grid, post_max_ratio) = sqrt(2.0);

  RSS(ref_collapse_pass(ref_grid), "col pass");
  ref_gather_blocking_frame(ref_grid, "collapse");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt col");

  ref_grid_adapt(ref_grid, post_max_ratio) =
      ref_grid_adapt(ref_grid, last_max_ratio);

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  RSS(ref_smooth_threed_pass(ref_grid), "smooth pass");
  ref_gather_blocking_frame(ref_grid, "smooth");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt move");

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  RSS(ref_adapt_parameter(ref_grid, &all_done1), "param");

  if (ref_grid_surf(ref_grid)) {
    RSS(ref_split_surf_pass(ref_grid), "split surfpass");
  } else {
    RSS(ref_split_pass(ref_grid), "split pass");
  }
  ref_gather_blocking_frame(ref_grid, "split");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt spl");

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  RSS(ref_smooth_threed_pass(ref_grid), "smooth pass");
  ref_gather_blocking_frame(ref_grid, "smooth");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt move");

  RSS(ref_adapt_threed_swap(ref_grid), "swap pass");
  ref_gather_blocking_frame(ref_grid, "swap");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  if (ref_grid_adapt(ref_grid, instrument))
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt swap");

  if (ngeom > 0)
    RSS(ref_geom_verify_topo(ref_grid), "geom topo postflight check");

  *all_done = (all_done0 && all_done1);

  return REF_SUCCESS;
}

static REF_STATUS ref_adapt_twod_pass(REF_GRID ref_grid, REF_BOOL *all_done) {
  RSS(ref_adapt_parameter(ref_grid, all_done), "param");

  ref_gather_blocking_frame(ref_grid, "twod pass");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  RSS(ref_smooth_twod_pass(ref_grid), "smooth pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "smooth");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  RSS(ref_collapse_twod_pass(ref_grid), "collapse pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "collapse");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  RSS(ref_smooth_twod_pass(ref_grid), "smooth pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "smooth");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  RSS(ref_adapt_parameter(ref_grid, all_done), "param");

  RSS(ref_split_twod_pass(ref_grid), "split pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "split");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  RSS(ref_smooth_twod_pass(ref_grid), "smooth pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "smooth");

  RSS(ref_swap_twod_pass(ref_grid), "swap pass");
  if (ref_grid_adapt(ref_grid, watch_param))
    RSS(ref_adapt_tattle(ref_grid), "tattle");
  ref_gather_blocking_frame(ref_grid, "swap");

  return REF_SUCCESS;
}

REF_STATUS ref_adapt_pass(REF_GRID ref_grid, REF_BOOL *all_done) {
  if (ref_grid_twod(ref_grid)) {
    RSS(ref_adapt_twod_pass(ref_grid, all_done), "2D pass");
  } else {
    RSS(ref_adapt_threed_pass(ref_grid, all_done), "3D pass");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_adapt_surf_to_geom(REF_GRID ref_grid) {
  REF_BOOL all_done = REF_FALSE;
  int passes = 15, pass;

  if (ref_mpi_para(ref_grid_mpi(ref_grid))) RSS(REF_IMPLEMENT, "seq only");

  RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "curvature");

  for (pass = 0; !all_done && pass < passes; pass++) {
    if (ref_grid_once(ref_grid))
      printf("\n pass %d of %d with %d ranks\n", pass + 1, passes,
             ref_mpi_n(ref_grid_mpi(ref_grid)));
    RSS(ref_adapt_pass(ref_grid, &all_done), "pass");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pass");
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "curvature");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "pack");
  }

  return REF_SUCCESS;
}
