
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

#ifndef REF_INFLATE_H
#define REF_INFLATE_H

#include "ref_defs.h"
#include "ref_dict.h"
#include "ref_grid.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

REF_STATUS ref_inflate_pri_min_dot(REF_NODE ref_node, REF_INT *nodes,
                                   REF_DBL *min_dot);

REF_STATUS ref_inflate_face(REF_GRID ref_grid, REF_DICT faceids,
                            REF_DBL *origin, REF_DBL thickness, REF_DBL xshift);

REF_STATUS ref_inflate_radially(REF_GRID ref_grid, REF_DICT faceids,
                                REF_DBL *origin, REF_DBL thickness,
                                REF_DBL mach_angle_rad, REF_DBL alpha_rad);

REF_STATUS ref_inflate_rate(REF_INT nlayers, REF_DBL first_thickness,
                            REF_DBL total_thickness, REF_DBL *rate);
REF_STATUS ref_inflate_total_thickness(REF_INT nlayers, REF_DBL first_thickness,
                                       REF_DBL rate, REF_DBL *total_thickness);
REF_STATUS ref_inflate_dthickness(REF_INT nlayers, REF_DBL first_thickness,
                                  REF_DBL rate, REF_DBL *dHdr);

REF_STATUS ref_inflate_origin(REF_GRID ref_grid, REF_DICT faceids,
                              REF_DBL *origin);
REF_STATUS ref_inflate_read_usm3d_mapbc(REF_DICT faceids,
                                        const char *mapbc_file_name,
                                        const char *family_name,
                                        REF_INT boundary_condition);
END_C_DECLORATION

#endif /* REF_INFLATE_H */
