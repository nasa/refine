
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

#ifndef REF_SWAP_H
#define REF_SWAP_H

#include "ref_defs.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_swap_remove_two_face_cell(REF_GRID ref_grid, REF_INT cell);
REF_STATUS ref_swap_remove_three_face_cell(REF_GRID ref_grid, REF_INT cell);
REF_STATUS ref_swap_pass(REF_GRID ref_grid);

/*   orig    swap
 *   0 - 2   0 - 2
 *   | \ |   | / |
 *   3 - 1   3 - 1
 */
REF_STATUS ref_swap_node23(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                           REF_INT *node2, REF_INT *node3);
REF_STATUS ref_swap_same_faceid(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_BOOL *allowed);
REF_STATUS ref_swap_manifold(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                             REF_BOOL *allowed);
REF_STATUS ref_swap_geom_topo(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                              REF_BOOL *allowed);
REF_STATUS ref_swap_local_cell(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                               REF_BOOL *allowed);
REF_STATUS ref_swap_conforming(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                               REF_BOOL *allowed);
REF_STATUS ref_swap_ratio(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_BOOL *allowed);
REF_STATUS ref_swap_quality(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                            REF_BOOL *allowed);

REF_STATUS ref_swap_surf_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1);

REF_STATUS ref_swap_surf_pass(REF_GRID ref_grid);
REF_STATUS ref_swap_twod_pass(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_SWAP_H */
