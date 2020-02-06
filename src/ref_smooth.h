
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

#ifndef REF_SMOOTH_H
#define REF_SMOOTH_H

#include "ref_defs.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_smooth_tri_ratio_around(REF_GRID ref_grid, REF_INT node,
                                       REF_DBL *min_ratio, REF_DBL *max_ratio);
REF_STATUS ref_smooth_tri_quality_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_quality);
REF_STATUS ref_smooth_tri_normdev_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_normdev);

REF_STATUS ref_smooth_outward_norm(REF_GRID ref_grid, REF_INT node,
                                   REF_BOOL *allowed);

REF_STATUS ref_smooth_tri_ideal(REF_GRID ref_grid, REF_INT node, REF_INT tri,
                                REF_DBL *ideal_location);

REF_STATUS ref_smooth_tri_ideal_uv(REF_GRID ref_grid, REF_INT node, REF_INT tri,
                                   REF_DBL *ideal_uv);

REF_STATUS ref_smooth_tri_weighted_ideal(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *ideal_location);

REF_STATUS ref_smooth_tri_weighted_ideal_uv(REF_GRID ref_grid, REF_INT node,
                                            REF_DBL *ideal_uv);

REF_STATUS ref_smooth_twod_boundary_nodes(REF_GRID ref_grid, REF_INT node,
                                          REF_INT *node0, REF_INT *node1);

REF_STATUS ref_smooth_twod_tri_improve(REF_GRID ref_grid, REF_INT node);
REF_STATUS ref_smooth_twod_tri_pliant(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_smooth_twod_pass(REF_GRID ref_grid);

REF_STATUS ref_smooth_tet_quality_around(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *min_quality);
REF_STATUS ref_smooth_tet_ratio_around(REF_GRID ref_grid, REF_INT node,
                                       REF_DBL *min_ratio, REF_DBL *max_ratio);

REF_STATUS ref_smooth_tet_ideal(REF_GRID ref_grid, REF_INT node, REF_INT tet,
                                REF_DBL *ideal_location);

REF_STATUS ref_smooth_tet_weighted_ideal(REF_GRID ref_grid, REF_INT node,
                                         REF_DBL *ideal_location);

REF_STATUS ref_smooth_tet_improve(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_smooth_geom_edge(REF_GRID ref_grid, REF_INT node);
REF_STATUS ref_smooth_geom_face(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_smooth_pass(REF_GRID ref_grid);

REF_STATUS ref_smooth_threed_post_edge_split(REF_GRID ref_grid, REF_INT node);
REF_STATUS ref_smooth_threed_post_face_split(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_smooth_nso(REF_GRID ref_grid, REF_INT node);

END_C_DECLORATION

#endif /* REF_SMOOTH_H */
