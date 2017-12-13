

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

#ifndef REF_INTERP_H
#define REF_INTERP_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_INTERP_STRUCT REF_INTERP_STRUCT;
typedef REF_INTERP_STRUCT * REF_INTERP;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_list.h"

BEGIN_C_DECLORATION
struct REF_INTERP_STRUCT {
  REF_GRID from_grid;
  REF_GRID to_grid;
  REF_BOOL instrument;
  REF_INT n_walk;
  REF_INT walk_steps;
  REF_INT n_geom;
  REF_INT n_geom_fail;
  REF_INT n_tree;
  REF_INT tree_cells;
  REF_INT *guess;
  REF_INT *cell;
  REF_INT *part;
  REF_DBL *bary;
  REF_DBL inside;
  REF_DBL bound;
  REF_LIST ref_list;
  REF_LIST visualize;
};

#define ref_interp_from_grid(ref_interp) ( (ref_interp)->from_grid )
#define ref_interp_to_grid(ref_interp) ( (ref_interp)->to_grid )

REF_STATUS ref_interp_create( REF_INTERP *ref_interp, 
			      REF_GRID from_grid, REF_GRID to_grid );

REF_STATUS ref_interp_free( REF_INTERP ref_interp );

REF_STATUS ref_interp_locate( REF_INTERP ref_interp);

REF_STATUS ref_interp_min_bary( REF_INTERP ref_interp, REF_DBL *min_bary );
REF_STATUS ref_interp_max_error( REF_INTERP ref_interp, REF_DBL *max_error );
REF_STATUS ref_interp_stats( REF_INTERP ref_interp );

REF_STATUS ref_interp_tec( REF_INTERP ref_interp, const char *filename );

#endif /* REF_INTERP_H */
