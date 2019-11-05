
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

#include <stdio.h>
#include <stdlib.h>

#include "ref_dist.h"

#include "ref_grid.h"
#include "ref_cell.h"

#include "ref_search.h"

REF_STATUS ref_dist_collisions(REF_GRID ref_grid) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_SEARCH ref_search;
  RSS(ref_search_create(&ref_search, ref_cell_n(ref_cell)), "create search");

  RSS(ref_search_free(ref_search), "free");    
  return REF_SUCCESS;
}
