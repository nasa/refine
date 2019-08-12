
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

#ifndef REF_COLLAPSE_H
#define REF_COLLAPSE_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

/* node1 is removed */

REF_STATUS ref_collapse_pass(REF_GRID ref_grid);

REF_STATUS ref_collapse_to_remove_node1(REF_GRID ref_grid, REF_INT *node0,
                                        REF_INT node1);

REF_STATUS ref_collapse_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1);
/*                               keep node0,  remove node1 */

REF_STATUS ref_collapse_edge_geometry(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_manifold(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_BOOL *allowed);
REF_STATUS ref_collapse_edge_chord_height(REF_GRID ref_grid, REF_INT node0,
                                          REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_same_normal(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_mixed(REF_GRID ref_grid, REF_INT node0,
                                   REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_local_cell(REF_GRID ref_grid, REF_INT node0,
                                        REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_cad_constrained(REF_GRID ref_grid, REF_INT node0,
                                             REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_quality(REF_GRID ref_grid, REF_INT node0,
                                     REF_INT node1, REF_BOOL *allowed);
REF_STATUS ref_collapse_edge_ratio(REF_GRID ref_grid, REF_INT node0,
                                   REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_edge_normdev(REF_GRID ref_grid, REF_INT node0,
                                     REF_INT node1, REF_BOOL *allowed);

REF_STATUS ref_collapse_face(REF_GRID ref_grid, REF_INT keep0, REF_INT remove0,
                             REF_INT keep1, REF_INT remove1);

REF_STATUS ref_collapse_face_local_pris(REF_GRID ref_grid, REF_INT keep,
                                        REF_INT remove, REF_BOOL *allowed);

REF_STATUS ref_collapse_face_quality(REF_GRID ref_grid, REF_INT keep,
                                     REF_INT remove, REF_BOOL *allowed);
REF_STATUS ref_collapse_face_ratio(REF_GRID ref_grid, REF_INT keep,
                                   REF_INT remove, REF_BOOL *allowed);

REF_STATUS ref_collapse_face_outward_norm(REF_GRID ref_grid, REF_INT keep,
                                          REF_INT remove, REF_BOOL *allowed);

REF_STATUS ref_collapse_face_geometry(REF_GRID ref_grid, REF_INT keep,
                                      REF_INT remove, REF_BOOL *allowed);

REF_STATUS ref_collapse_face_same_tangent(REF_GRID ref_grid, REF_INT keep,
                                          REF_INT remove, REF_BOOL *allowed);

REF_STATUS ref_collapse_twod_pass(REF_GRID ref_grid);

REF_STATUS ref_collapse_face_remove_node1(REF_GRID ref_grid,
                                          REF_INT *actual_node0, REF_INT node1);

END_C_DECLORATION

#endif /* REF_COLLAPSE_H */
