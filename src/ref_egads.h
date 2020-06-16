
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

#ifndef REF_EGADS_H
#define REF_EGADS_H

#include "ref_defs.h"
BEGIN_C_DECLORATION

#define REF_EGADS_MISSING_TPARAM (1)
#define REF_EGADS_CHORD_TPARAM (2)
#define REF_EGADS_WIDTH_TPARAM (4)
#define REF_EGADS_ALL_TPARAM (4 + 2 + 1)

END_C_DECLORATION

#include "ref_geom.h"

BEGIN_C_DECLORATION

REF_STATUS ref_egads_open(REF_GEOM ref_geom);
REF_STATUS ref_egads_close(REF_GEOM ref_geom);
REF_STATUS ref_egads_out_level(REF_GEOM ref_geom, REF_INT out_level);
REF_STATUS ref_egads_load(REF_GEOM ref_geom, const char *filename);

REF_STATUS ref_egads_tess(REF_GRID ref_grid, REF_INT auto_tparams,
                          REF_DBL *global_params);

REF_STATUS ref_egads_mark_jump_degen(REF_GRID ref_grid);

REF_STATUS ref_egads_recon(REF_GRID ref_grid);

REF_STATUS ref_egads_diagonal(REF_GEOM ref_geom, REF_INT geom, REF_DBL *diag);
REF_STATUS ref_egads_tolerance(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                               REF_DBL *tolerance);

REF_STATUS ref_egads_edge_curvature(REF_GEOM ref_geom, REF_INT geom, REF_DBL *k,
                                    REF_DBL *normal);
REF_STATUS ref_egads_face_curvature(REF_GEOM ref_geom, REF_INT geom,
                                    REF_DBL *kr, REF_DBL *r, REF_DBL *ks,
                                    REF_DBL *s);
REF_STATUS ref_egads_edge_trange(REF_GEOM ref_geom, REF_INT id,
                                 REF_DBL *trange);
REF_STATUS ref_egads_edge_face_uv(REF_GEOM ref_geom, REF_INT edgeid,
                                  REF_INT faceid, REF_INT sense, REF_DBL t,
                                  REF_DBL *uv);

REF_STATUS ref_egads_eval(REF_GEOM ref_geom, REF_INT geom, REF_DBL *xyz,
                          REF_DBL *dxyz_dtuv);
REF_STATUS ref_egads_eval_at(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                             REF_DBL *params, REF_DBL *xyz, REF_DBL *dxyz_dtuv);
REF_STATUS ref_egads_inverse_eval(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                                  REF_DBL *xyz, REF_DBL *param);

REF_STATUS ref_egads_gap(REF_GEOM ref_geom, REF_INT node, REF_DBL *gap);

REF_STATUS ref_egads_feature_size(REF_GRID ref_grid, REF_INT node, REF_DBL *h0,
                                 REF_DBL *dir0, REF_DBL *h1, REF_DBL *dir1,
                                 REF_DBL *h2, REF_DBL *dir2);

END_C_DECLORATION

#endif /* REF_EGADS_H */
