
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

#ifndef REF_SPLIT_H
#define REF_SPLIT_H

#include "ref_defs.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_split_surf_pass(REF_GRID ref_grid);
REF_STATUS ref_split_pass(REF_GRID ref_grid);

REF_STATUS ref_split_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_INT new_node);
REF_STATUS ref_split_face(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_INT node2, REF_INT new_node);

REF_STATUS ref_split_edge_mixed(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_BOOL *allowed);

REF_STATUS ref_split_edge_tet_quality(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_INT new_node,
                                      REF_BOOL *allowed);
REF_STATUS ref_split_edge_tet_ratio(REF_GRID ref_grid, REF_INT node0,
                                    REF_INT node1, REF_INT new_node,
                                    REF_BOOL *allowed);
REF_STATUS ref_split_edge_tri_conformity(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_INT new_node,
                                         REF_BOOL *allowed);

REF_STATUS ref_split_twod_pass(REF_GRID ref_grid);

REF_STATUS ref_split_twod_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                               REF_INT new_node);

REF_STATUS ref_split_prism_tri_quality(REF_GRID ref_grid, REF_INT node0,
                                       REF_INT node1, REF_INT new_node,
                                       REF_BOOL *allowed);
REF_STATUS ref_split_prism_tri_ratio(REF_GRID ref_grid, REF_INT node0,
                                     REF_INT node1, REF_INT new_node,
                                     REF_BOOL *allowed);

REF_STATUS ref_split_edge_pattern(REF_GRID ref_grid, REF_INT first,
                                  REF_INT skip);

REF_STATUS ref_split_edge_geometry(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_SPLIT_H */
