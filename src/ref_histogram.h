
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

#ifndef REF_HISTOGRAM_H
#define REF_HISTOGRAM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_HISTOGRAM_STRUCT REF_HISTOGRAM_STRUCT;
typedef REF_HISTOGRAM_STRUCT *REF_HISTOGRAM;
END_C_DECLORATION

#include <stdio.h>

#include "ref_grid.h"
#include "ref_mpi.h"

BEGIN_C_DECLORATION
struct REF_HISTOGRAM_STRUCT {
  REF_DBL max, min;
  REF_DBL log_total;
  REF_DBL log_mean;
  REF_DBL exp;
  REF_INT nbin;
  REF_INT *bins;
  REF_INT nstat;
  REF_DBL *stats;
  FILE *df;
};

REF_STATUS ref_histogram_create(REF_HISTOGRAM *ref_histogram);
REF_STATUS ref_histogram_free(REF_HISTOGRAM ref_histogram);

#define ref_histogram_max(ref_histogram) ((ref_histogram)->max)
#define ref_histogram_min(ref_histogram) ((ref_histogram)->min)
#define ref_histogram_log_total(ref_histogram) ((ref_histogram)->log_total)
#define ref_histogram_log_mean(ref_histogram) ((ref_histogram)->log_mean)
#define ref_histogram_exp(ref_histogram) ((ref_histogram)->exp)
#define ref_histogram_nbin(ref_histogram) ((ref_histogram)->nbin)
#define ref_histogram_bin(ref_histogram, i) ((ref_histogram)->bins[(i)])
#define ref_histogram_nstat(ref_histogram) ((ref_histogram)->nstat)
#define ref_histogram_stat(ref_histogram, i) ((ref_histogram)->stats[(i)])

REF_INT ref_histogram_to_bin(REF_HISTOGRAM ref_histogram, REF_DBL observation);

#define ref_histogram_to_obs(i)                        \
  (pow(2.0, (1.0 / ref_histogram_exp(ref_histogram)) * \
                ((REF_DBL)((i)-ref_histogram_nbin(ref_histogram) / 2))))

/* defaults are 18, 6.0 resolved 72, 12.0 */
REF_STATUS ref_histogram_resolution(REF_HISTOGRAM ref_histogram, REF_INT nbin,
                                    REF_DBL exp);
REF_STATUS ref_histogram_add(REF_HISTOGRAM ref_histogram, REF_DBL observation);
REF_STATUS ref_histogram_gather(REF_HISTOGRAM ref_histogram, REF_MPI ref_mpi);

REF_STATUS ref_histogram_print(REF_HISTOGRAM ref_histogram, REF_GRID ref_grid,
                               const char *description);
REF_STATUS ref_histogram_gnuplot(REF_HISTOGRAM ref_histogram,
                                 const char *description);
REF_STATUS ref_histogram_tec(REF_HISTOGRAM ref_histogram,
                             const char *description);
REF_STATUS ref_histogram_zone(REF_HISTOGRAM ref_histogram, FILE *file,
                              const char *zone_title, REF_DBL time);

REF_STATUS ref_histogram_add_stat(REF_HISTOGRAM ref_histogram,
                                  REF_DBL observation);
REF_STATUS ref_histogram_gather_stat(REF_HISTOGRAM ref_histogram,
                                     REF_MPI ref_mpi);
REF_STATUS ref_histogram_print_stat(REF_HISTOGRAM ref_histogram);

REF_STATUS ref_histogram_add_ratio(REF_HISTOGRAM ref_histogram,
                                   REF_GRID ref_grid);
REF_STATUS ref_histogram_add_quality(REF_HISTOGRAM ref_histogram,
                                     REF_GRID ref_grid);
REF_STATUS ref_histogram_ratio(REF_GRID ref_grid);
REF_STATUS ref_histogram_quality(REF_GRID ref_grid);

REF_STATUS ref_histogram_ratio_tec(REF_GRID ref_grid);
REF_STATUS ref_histogram_quality_tec(REF_GRID ref_grid);
REF_STATUS ref_histogram_node_tec(REF_GRID ref_grid, REF_DBL *observations);

REF_STATUS ref_histogram_debug(REF_HISTOGRAM ref_histogram,
                               const char *filename);

END_C_DECLORATION

#endif /* REF_HISTOGRAM_H */
