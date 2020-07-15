
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

#include "ref_histogram.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"

REF_STATUS ref_histogram_create(REF_HISTOGRAM *ref_histogram_ptr) {
  REF_HISTOGRAM ref_histogram;

  ref_malloc(*ref_histogram_ptr, 1, REF_HISTOGRAM_STRUCT);
  ref_histogram = (*ref_histogram_ptr);

  ref_histogram_nbin(ref_histogram) = 18;
  ref_histogram_exp(ref_histogram) = 6.0;

  ref_malloc_init(ref_histogram->bins, ref_histogram_nbin(ref_histogram),
                  REF_INT, 0);

  ref_histogram_max(ref_histogram) = -1.0e20;
  ref_histogram_min(ref_histogram) = 1.0e20;
  ref_histogram_log_total(ref_histogram) = 0.0;
  ref_histogram_log_mean(ref_histogram) = 0.0;

  ref_histogram_nstat(ref_histogram) = 5;

  ref_malloc_init(ref_histogram->stats, ref_histogram_nstat(ref_histogram),
                  REF_DBL, 0.0);

  ref_histogram->df = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_resolution(REF_HISTOGRAM ref_histogram, REF_INT nbin,
                                    REF_DBL exp) {
  ref_histogram_nbin(ref_histogram) = nbin;
  ref_histogram_exp(ref_histogram) = exp;

  ref_free(ref_histogram->bins);
  ref_malloc_init(ref_histogram->bins, ref_histogram_nbin(ref_histogram),
                  REF_INT, 0);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_free(REF_HISTOGRAM ref_histogram) {
  if (NULL == (void *)ref_histogram) return REF_NULL;
  if (NULL != ref_histogram->df) fclose(ref_histogram->df);
  ref_free(ref_histogram->stats);
  ref_free(ref_histogram->bins);
  ref_free(ref_histogram);
  return REF_SUCCESS;
}

REF_INT ref_histogram_to_bin(REF_HISTOGRAM ref_histogram, REF_DBL observation) {
  REF_DBL dbin;
  REF_INT ibin;
  dbin = ref_histogram_exp(ref_histogram) * log2((observation));
  dbin += ref_histogram_nbin(ref_histogram) / 2;
  ibin = (REF_INT)dbin;
  ibin = MAX(0, ibin);
  ibin = MIN(ref_histogram_nbin(ref_histogram) - 1, ibin);
  return ibin;
}

REF_STATUS ref_histogram_add(REF_HISTOGRAM ref_histogram, REF_DBL observation) {
  REF_INT i;

  if (observation <= 0.0) return REF_INVALID;

  ref_histogram_max(ref_histogram) =
      MAX(ref_histogram_max(ref_histogram), observation);
  ref_histogram_min(ref_histogram) =
      MIN(ref_histogram_min(ref_histogram), observation);
  ref_histogram_log_total(ref_histogram) += log2(observation);

  i = ref_histogram_to_bin(ref_histogram, observation);
  i = MIN(i, ref_histogram_nbin(ref_histogram) - 1);
  i = MAX(i, 0);

  ref_histogram_bin(ref_histogram, i)++;

  if (NULL != ref_histogram->df)
    fprintf(ref_histogram->df, "%f\n", observation);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_gather(REF_HISTOGRAM ref_histogram, REF_MPI ref_mpi) {
  REF_INT *bins;
  REF_DBL min, max, log_total;
  REF_INT i, observations;

  ref_malloc(bins, ref_histogram_nbin(ref_histogram), REF_INT);

  RSS(ref_mpi_sum(ref_mpi, ref_histogram->bins, bins,
                  ref_histogram_nbin(ref_histogram), REF_INT_TYPE),
      "sum");
  RSS(ref_mpi_max(ref_mpi, &ref_histogram_max(ref_histogram), &max,
                  REF_DBL_TYPE),
      "max");
  RSS(ref_mpi_min(ref_mpi, &ref_histogram_min(ref_histogram), &min,
                  REF_DBL_TYPE),
      "min");
  RSS(ref_mpi_sum(ref_mpi, &ref_histogram_log_total(ref_histogram), &log_total,
                  1, REF_DBL_TYPE),
      "log_total");

  if (ref_mpi_once(ref_mpi)) {
    ref_histogram_max(ref_histogram) = max;
    ref_histogram_min(ref_histogram) = min;
    ref_histogram_log_total(ref_histogram) = log_total;
    for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
      ref_histogram->bins[i] = bins[i];
    observations = 0;
    for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
      observations += ref_histogram->bins[i];
    ref_histogram_log_mean(ref_histogram) = 0.0;
    if (observations > 0)
      ref_histogram_log_mean(ref_histogram) = log_total / (REF_DBL)observations;
  } else {
    ref_histogram_max(ref_histogram) = -1.0e20;
    ref_histogram_min(ref_histogram) = 1.0e20;
    ref_histogram_log_total(ref_histogram) = 0.0;
    for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
      ref_histogram->bins[i] = 0;
    ref_histogram_log_mean(ref_histogram) = 0.0;
  }

  RSS(ref_mpi_bcast(ref_mpi, &ref_histogram_log_mean(ref_histogram), 1,
                    REF_DBL_TYPE),
      "log_mean");

  ref_free(bins);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_print(REF_HISTOGRAM ref_histogram, REF_GRID ref_grid,
                               const char *description) {
  REF_INT i, sum;
  REF_DBL log_mean;

  sum = 0;
  for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
    sum += ref_histogram_bin(ref_histogram, i);

  printf("%7.3f %10.3e min %s\n", ref_histogram_min(ref_histogram),
         ref_histogram_min(ref_histogram), description);

  for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
    if (ref_histogram_to_obs(i + 1) > ref_histogram_min(ref_histogram) &&
        ref_histogram_to_obs(i - 1) < ref_histogram_max(ref_histogram)) {
      if ((ref_histogram_to_obs(i) > ref_grid_adapt(ref_grid, split_ratio) ||
           ref_histogram_to_obs(i - 1) <
               ref_grid_adapt(ref_grid, collapse_ratio)) &&
          ref_histogram_bin(ref_histogram, i) > 0) {
        printf("%7.3f:%10d *\n", ref_histogram_to_obs(i),
               ref_histogram_bin(ref_histogram, i));
      } else {
        printf("%7.3f:%10d\n", ref_histogram_to_obs(i),
               ref_histogram_bin(ref_histogram, i));
      }
    }

  printf("%7.3f %10.3e:%10d max %s\n", ref_histogram_max(ref_histogram),
         ref_histogram_max(ref_histogram), sum, description);
  log_mean = ref_histogram_log_total(ref_histogram) / (REF_DBL)sum;
  printf("%18.10f mean %s\n", pow(2.0, log_mean), description);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_gnuplot(REF_HISTOGRAM ref_histogram,
                                 const char *description) {
  REF_INT i;
  FILE *f;
  char filename[1024];
  REF_INT sum;
  REF_DBL norm, portion;

  sum = 0;
  for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
    sum += ref_histogram_bin(ref_histogram, i);
  norm = 0.0;
  if (0 < sum) norm = 1.0 / (REF_DBL)sum;

  sprintf(filename, "ref_histogram_%s.gnuplot", description);
  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "set terminal postscript eps enhanced color\n");
  sprintf(filename, "ref_histogram_%s.eps", description);
  fprintf(f, "set output '%s'\n", filename);

  fprintf(f, "set style data histogram\n");
  fprintf(f, "set style histogram cluster gap 1\n");
  fprintf(f, "set style fill solid border -1\n");
  fprintf(f, "set boxwidth 0.9\n");
  fprintf(f, "set xtic rotate by -45 scale 0\n");
  fprintf(f, "set yrange [0:0.4]\n");
  fprintf(f, "set format y \"%%.2f\"\n");
  fprintf(
      f,
      "plot '-' using 2:xticlabels(1) title '%s' linecolor rgb \"#FF0000\"\n",
      description);

  for (i = 0; i < ref_histogram_nbin(ref_histogram) - 1; i++) {
    portion = (REF_DBL)ref_histogram_bin(ref_histogram, i) * norm;
    fprintf(f, "%.2f-%.2f %.5f\n", ref_histogram_to_obs(i),
            ref_histogram_to_obs(i + 1), portion);
  }

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_tec(REF_HISTOGRAM ref_histogram,
                             const char *description) {
  REF_INT i;
  FILE *f;
  char filename[1024];
  REF_INT sum;
  REF_DBL norm, area, portion;

  sum = 0;
  for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
    sum += ref_histogram_bin(ref_histogram, i);
  norm = 0.0;
  if (0 < sum) norm = 1.0 / (REF_DBL)sum;

  sprintf(filename, "ref_histogram_%s.tec", description);
  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  for (i = 0; i < ref_histogram_nbin(ref_histogram) - 2; i++) {
    area = ref_histogram_to_obs(i + 1) - ref_histogram_to_obs(i);
    portion = (REF_DBL)ref_histogram_bin(ref_histogram, i) / area * norm;
    portion = MAX(portion, 1.0e-20);
    fprintf(f, "%.8e %.8e\n%.8e %.8e\n", ref_histogram_to_obs(i), portion,
            ref_histogram_to_obs(i + 1), portion);
  }

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_zone(REF_HISTOGRAM ref_histogram, FILE *file,
                              const char *zone_title, REF_DBL time) {
  REF_INT i;
  REF_INT sum;
  REF_DBL norm, area, portion, center;

  sum = 0;
  for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
    sum += ref_histogram_bin(ref_histogram, i);
  norm = 0.0;
  if (0 < sum) norm = 1.0 / (REF_DBL)sum;

  if (NULL == zone_title) {
    fprintf(file, "ZONE I=%d, DATAPACKING=POINT, SOLUTIONTIME=%f\n",
            ref_histogram_nbin(ref_histogram) - 2, time);
  } else {
    fprintf(file, "ZONE T=\"%s\", I=%d, DATAPACKING=POINT, SOLUTIONTIME=%f\n",
            zone_title, ref_histogram_nbin(ref_histogram) - 2, time);
  }
  for (i = 0; i < ref_histogram_nbin(ref_histogram) - 2; i++) {
    area = ref_histogram_to_obs(i + 1) - ref_histogram_to_obs(i);
    portion = (REF_DBL)ref_histogram_bin(ref_histogram, i) / area * norm;
    portion = MAX(portion, 1.0e-20);
    center = 0.5 * (ref_histogram_to_obs(i) + ref_histogram_to_obs(i + 1));
    fprintf(file, "%.8e %.8e\n", center, portion);
  }
  REIS(0, fflush(file), "tec zone fflush");

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_add_stat(REF_HISTOGRAM ref_histogram,
                                  REF_DBL observation) {
  REF_INT i;
  REF_DBL log_obs;

  if (observation <= 0.0) return REF_INVALID;

  log_obs = log2(observation);

  for (i = 0; i < ref_histogram_nstat(ref_histogram); i++)
    ref_histogram_stat(ref_histogram, i) +=
        pow(log_obs - ref_histogram_log_mean(ref_histogram), i);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_gather_stat(REF_HISTOGRAM ref_histogram,
                                     REF_MPI ref_mpi) {
  REF_DBL *stats;
  REF_INT i, observations;

  ref_malloc(stats, ref_histogram_nstat(ref_histogram), REF_DBL);

  RSS(ref_mpi_sum(ref_mpi, ref_histogram->stats, stats,
                  ref_histogram_nstat(ref_histogram), REF_DBL_TYPE),
      "sum");

  if (ref_mpi_once(ref_mpi)) {
    observations = 0;
    for (i = 0; i < ref_histogram_nbin(ref_histogram); i++)
      observations += ref_histogram->bins[i];
    if (observations > 0)
      for (i = 0; i < ref_histogram_nstat(ref_histogram); i++)
        ref_histogram_stat(ref_histogram, i) = stats[i] / (REF_DBL)observations;
  } else {
    for (i = 0; i < ref_histogram_nstat(ref_histogram); i++)
      ref_histogram_stat(ref_histogram, i) = 0.0;
  }

  ref_free(stats);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_print_stat(REF_HISTOGRAM ref_histogram) {
  REF_INT i;
  REF_DBL skewness, kurtosis, n, d;

  printf("mom ");
  for (i = 0; i < ref_histogram_nstat(ref_histogram); i++)
    printf("%12.8f(%d)", ref_histogram_stat(ref_histogram, i), i);
  printf("\n");

  n = ref_histogram_stat(ref_histogram, 3);
  d = pow(ref_histogram_stat(ref_histogram, 2), 1.5);
  skewness = 0.0;
  if (ref_math_divisible(n, d)) skewness = n / d;

  n = ref_histogram_stat(ref_histogram, 4);
  d = pow(ref_histogram_stat(ref_histogram, 2), 2);
  kurtosis = 0.0;
  if (ref_math_divisible(n, d)) kurtosis = n / d - 3.0;

  printf("skewness %12.8f kurtosis %12.8f\n", skewness, kurtosis);

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_add_ratio(REF_HISTOGRAM ref_histogram,
                                   REF_GRID ref_grid) {
  REF_EDGE ref_edge;
  REF_INT edge, part;
  REF_DBL ratio;

  RSS(ref_edge_create(&ref_edge, ref_grid), "make edges");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part == ref_mpi_rank(ref_grid_mpi(ref_grid))) {
      RSS(ref_node_ratio(ref_grid_node(ref_grid),
                         ref_edge_e2n(ref_edge, 0, edge),
                         ref_edge_e2n(ref_edge, 1, edge), &ratio),
          "rat");
      RSB(ref_histogram_add(ref_histogram, ratio), "add", {
        printf("ratio %e at %f %f %f\n", ratio,
               ref_node_xyz(ref_grid_node(ref_grid), 0,
                            ref_edge_e2n(ref_edge, 0, edge)),
               ref_node_xyz(ref_grid_node(ref_grid), 1,
                            ref_edge_e2n(ref_edge, 0, edge)),
               ref_node_xyz(ref_grid_node(ref_grid), 2,
                            ref_edge_e2n(ref_edge, 0, edge)));
      });
    }
  }

  RSS(ref_histogram_gather(ref_histogram, ref_grid_mpi(ref_grid)), "gather");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part == ref_mpi_rank(ref_grid_mpi(ref_grid))) {
      RSS(ref_node_ratio(ref_grid_node(ref_grid),
                         ref_edge_e2n(ref_edge, 0, edge),
                         ref_edge_e2n(ref_edge, 1, edge), &ratio),
          "rat");
      RSS(ref_histogram_add_stat(ref_histogram, ratio), "add");
    }
  }
  RSS(ref_histogram_gather_stat(ref_histogram, ref_grid_mpi(ref_grid)),
      "gather");

  RSS(ref_edge_free(ref_edge), "free edge");

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_add_quality(REF_HISTOGRAM ref_histogram,
                                     REF_GRID ref_grid) {
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality;

  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_node_part(ref_grid_node(ref_grid), nodes[0]) ==
        ref_mpi_rank(ref_grid_mpi(ref_grid))) {
      if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid)) {
        RSS(ref_node_tri_quality(ref_grid_node(ref_grid), nodes, &quality),
            "qual");
      } else {
        RSS(ref_node_tet_quality(ref_grid_node(ref_grid), nodes, &quality),
            "qual");
      }
      if (quality > 0.0) RSS(ref_histogram_add(ref_histogram, quality), "add");
    }
  }

  RSS(ref_histogram_gather(ref_histogram, ref_grid_mpi(ref_grid)), "gather");

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_ratio(REF_GRID ref_grid) {
  REF_HISTOGRAM ref_histogram;

  RSS(ref_histogram_create(&ref_histogram), "create");

  if (REF_FALSE)
    RSS(ref_histogram_debug(ref_histogram, "ref_histogram.len"), "dbug");

  RSS(ref_histogram_add_ratio(ref_histogram, ref_grid), "add ratio");

  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_print(ref_histogram, ref_grid, "edge ratio"), "print");
  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_print_stat(ref_histogram), "pr stat");

  RSS(ref_histogram_free(ref_histogram), "free gram");
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_quality(REF_GRID ref_grid) {
  REF_HISTOGRAM ref_histogram;

  RSS(ref_histogram_create(&ref_histogram), "create");

  RSS(ref_histogram_add_quality(ref_histogram, ref_grid), "add quality");

  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_print(ref_histogram, ref_grid, "quality"), "print");

  RSS(ref_histogram_free(ref_histogram), "free gram");
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_ratio_tec(REF_GRID ref_grid) {
  REF_HISTOGRAM ref_histogram;

  RSS(ref_histogram_create(&ref_histogram), "create");
  RSS(ref_histogram_resolution(ref_histogram, 288, 12.0), "res");

  RSS(ref_histogram_add_ratio(ref_histogram, ref_grid), "add ratio");

  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_tec(ref_histogram, "ratio"), "tec");

  RSS(ref_histogram_free(ref_histogram), "free gram");
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_quality_tec(REF_GRID ref_grid) {
  REF_HISTOGRAM ref_histogram;

  RSS(ref_histogram_create(&ref_histogram), "create");
  RSS(ref_histogram_resolution(ref_histogram, 288, 12.0), "res");

  RSS(ref_histogram_add_quality(ref_histogram, ref_grid), "add ratio");

  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_tec(ref_histogram, "quality"), "tec");

  RSS(ref_histogram_free(ref_histogram), "free gram");
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_node_tec(REF_GRID ref_grid, REF_DBL *observations) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_HISTOGRAM ref_histogram;
  REF_INT node;

  RSS(ref_histogram_create(&ref_histogram), "create");
  RSS(ref_histogram_resolution(ref_histogram, 288, 12.0), "res");

  each_ref_node_valid_node(ref_node, node) {
    if (observations[node] > 0.0)
      RSS(ref_histogram_add(ref_histogram, observations[node]), "add");
  }

  RSS(ref_histogram_gather(ref_histogram, ref_grid_mpi(ref_grid)), "gather");

  if (ref_grid_once(ref_grid))
    RSS(ref_histogram_tec(ref_histogram, "node"), "tec");

  RSS(ref_histogram_free(ref_histogram), "free gram");
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_debug(REF_HISTOGRAM ref_histogram,
                               const char *filename) {
  ref_histogram->df = fopen(filename, "w");
  if (NULL == (void *)ref_histogram->df)
    printf("unable to open %s\n", filename);
  RNS(ref_histogram->df, "unable to open file");

  return REF_SUCCESS;
}
