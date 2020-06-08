
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

#ifndef REF_MESHLINK_H
#define REF_MESHLINK_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_meshlink_open(REF_GRID ref_grid, const char *xml_filename);
REF_STATUS ref_meshlink_parse(REF_GRID ref_grid, const char *geom_filename);
REF_STATUS ref_meshlink_link(REF_GRID ref_grid, const char *base_name);
REF_STATUS ref_meshlink_constrain(REF_GRID ref_grid, REF_INT node);
REF_STATUS ref_meshlink_gap(REF_GRID ref_grid, REF_INT node, REF_DBL *gap);
REF_STATUS ref_meshlink_tri_norm_deviation(REF_GRID ref_grid, REF_INT *nodes,
                                           REF_DBL *dot_product);
REF_STATUS ref_meshlink_edge_curvature(REF_GRID ref_grid, REF_INT geom,
                                       REF_DBL *k, REF_DBL *normal);
REF_STATUS ref_meshlink_face_curvature(REF_GRID ref_grid, REF_INT geom,
                                       REF_DBL *kr, REF_DBL *r, REF_DBL *ks,
                                       REF_DBL *s);
REF_STATUS ref_meshlink_close(REF_GRID ref_grid);

REF_STATUS ref_meshlink_infer_orientation(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_MESHLINK_H */
