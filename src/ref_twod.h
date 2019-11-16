
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

#ifndef REF_TWOD_H
#define REF_TWOD_H

#include "ref_cell.h"
#include "ref_defs.h"

BEGIN_C_DECLORATION

REF_STATUS ref_twod_opposite_node(REF_CELL pri, REF_INT node,
                                  REF_INT *opposite);

REF_STATUS ref_twod_opposite_edge(REF_CELL pri, REF_INT node0, REF_INT node1,
                                  REF_INT *node2, REF_INT *node3);

REF_STATUS ref_twod_tri_pri_tri(REF_CELL tri, REF_CELL pri, REF_INT cell,
                                REF_INT *pri_cell, REF_INT *tri_cell);

END_C_DECLORATION

#endif /* REF_TWOD_H */
