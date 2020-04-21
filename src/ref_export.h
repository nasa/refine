
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
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_export_by_extension(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_vtk(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_tec(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_tec_surf(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_tec_cubic_edge_zone(REF_GRID ref_grid, FILE *file);
REF_STATUS ref_export_tec_edge_zone(REF_GRID ref_grid, FILE *file);
REF_STATUS ref_export_tec_surf_zone(REF_GRID ref_grid, FILE *file);
REF_STATUS ref_export_tec_vol_zone(REF_GRID ref_grid, FILE *file);

REF_STATUS ref_export_tec_int(REF_GRID ref_grid, REF_INT *scalar,
                              const char *filename);
REF_STATUS ref_export_tec_dbl(REF_GRID ref_grid, REF_INT ldim, REF_DBL *scalar,
                              const char *filename);

REF_STATUS ref_export_tec_part(REF_GRID ref_grid, const char *root_filename);
REF_STATUS ref_export_metric_xyzdirlen(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_tec_metric_axis(REF_GRID ref_grid,
                                      const char *root_filename);
REF_STATUS ref_export_tec_metric_ellipse(REF_GRID ref_grid,
                                         const char *root_filename);
REF_STATUS ref_export_tec_ratio(REF_GRID ref_grid, const char *root_filename);

REF_STATUS ref_export_poly(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_smesh(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_fgrid(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_su2(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_cogsg(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_c(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_eps(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_pdf(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_html(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_meshb_next_position(FILE *file, REF_INT version,
                                          REF_FILEPOS next_position);
REF_STATUS ref_export_meshb(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_twod_msh(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_twod_sol(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_export_plt(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_export_plt_tet_zone(REF_GRID ref_grid, FILE *file);
REF_STATUS ref_export_plt_surf_zone(REF_GRID ref_grid, FILE *file);

REF_STATUS ref_export_order_segments(REF_INT n, REF_INT *c2n, REF_INT *order);

END_C_DECLORATION

#endif /* REF_EXPORT_H */
