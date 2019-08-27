
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

#ifndef REF_GATHER_H
#define REF_GATHER_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_GATHER_STRUCT REF_GATHER_STRUCT;
typedef REF_GATHER_STRUCT *REF_GATHER;
END_C_DECLORATION

#include "ref_cell.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_node.h"

#include "ref_mpi.h"

BEGIN_C_DECLORATION
struct REF_GATHER_STRUCT {
  REF_BOOL recording;
  FILE *grid_file;
  FILE *hist_file;
  REF_DBL time;
  REF_BOOL low_quality_zone;
  REF_DBL min_quality;
};

#define ref_gather_low_quality_zone(ref_gather) ((ref_gather)->low_quality_zone)
#define ref_gather_min_quality(ref_gather) ((ref_gather)->min_quality)

REF_STATUS ref_gather_create(REF_GATHER *ref_gather);
REF_STATUS ref_gather_free(REF_GATHER ref_gather);

#define ref_gather_blocking_frame(ref_grid, zone_title) \
  RSS(ref_gather_tec_movie_frame(ref_grid, zone_title), "movie frame")

REF_STATUS ref_gather_tec_movie_record_button(REF_GATHER ref_gather,
                                              REF_BOOL on_or_off);
REF_STATUS ref_gather_tec_movie_frame(REF_GRID ref_grid,
                                      const char *zone_title);

REF_STATUS ref_gather_tec_part(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_gather_by_extension(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_gather_metric(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_gather_scalar(REF_GRID ref_grid, REF_INT ldim, REF_DBL *scalar,
                             const char *filename);

REF_STATUS ref_gather_ncell(REF_NODE ref_node, REF_CELL ref_cell,
                            REF_LONG *ncell);
REF_STATUS ref_gather_ngeom(REF_NODE ref_node, REF_GEOM ref_geom, REF_INT type,
                            REF_INT *ngeom);

REF_STATUS ref_gather_scalar_tec(REF_GRID ref_grid, REF_INT ldim,
                                 REF_DBL *scalar, const char **scalar_names,
                                 const char *filename);
REF_STATUS ref_gather_scalar_surf_tec(REF_GRID ref_grid, REF_INT ldim,
                                      REF_DBL *scalar,
                                      const char **scalar_names,
                                      const char *filename);

REF_STATUS ref_gather_scalar_by_extension(REF_GRID ref_grid, REF_INT ldim,
                                          REF_DBL *scalar,
                                          const char **scalar_names,
                                          const char *filename);
REF_STATUS ref_gather_surf_status_tec(REF_GRID ref_grid, const char *filename);

END_C_DECLORATION

#endif /* REF_GATHER_H */
