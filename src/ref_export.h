
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

#ifndef REF_EXPORT_H
#define REF_EXPORT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
/* these values requested to emulate Feflo.a behavior, should be all one? */
#define REF_EXPORT_MESHB_VERTEX_ID (1)
#define REF_EXPORT_MESHB_2D_ID (1)
#define REF_EXPORT_MESHB_3D_ID (0)
#define REF_EXPORT_MESHB_VERTEX_3 (10000000)
#define REF_EXPORT_MESHB_VERTEX_4 (200000000)
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_export_tec_surf(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_tec_metric_axis(REF_GRID ref_grid,
                                      const char *root_filename);
REF_STATUS ref_export_tec_metric_ellipse(REF_GRID ref_grid,
                                         const char *root_filename);
REF_STATUS ref_export_tec_ratio(REF_GRID ref_grid, const char *root_filename);

REF_STATUS ref_export_meshb_next_position(FILE *file, REF_INT version,
                                          REF_FILEPOS next_position);

REF_STATUS ref_export_order_segments(REF_INT n, REF_INT *c2n, REF_INT *order);

REF_STATUS ref_export_by_extension(REF_GRID ref_grid, const char *filename);

END_C_DECLORATION

#endif /* REF_EXPORT_H */
