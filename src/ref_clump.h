
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

#ifndef REF_CLUMP_H
#define REF_CLUMP_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_clump_around(REF_GRID ref_grid, REF_INT node,
                            const char *filename);
REF_STATUS ref_clump_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                             const char *filename);
REF_STATUS ref_clump_between_export_to(REF_GRID ref_grid, REF_INT node0,
                                       REF_INT node1, const char *filename);
REF_STATUS ref_clump_tri_around(REF_GRID ref_grid, REF_INT node,
                                const char *filename);
REF_STATUS ref_clump_short_edges(REF_GRID ref_grid, REF_DBL ratio_tol);
REF_STATUS ref_clump_short_edges_twod(REF_GRID ref_grid);

REF_STATUS ref_clump_long_edges(REF_GRID ref_grid, REF_DBL ratio_tol);

REF_STATUS ref_clump_tet_quality(REF_GRID ref_grid, REF_DBL min_quality,
                                 const char *filename);

END_C_DECLORATION

#endif /* REF_CLUMP_H */
