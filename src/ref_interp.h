

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

BEGIN_C_DECLORATION
struct REF_INTERP_STRUCT {
  REF_INT steps;
};

#define ref_interp_steps(ref_interp) ( (ref_interp)->steps )

REF_STATUS ref_interp_create( REF_INTERP *ref_interp );

REF_STATUS ref_interp_free( REF_INTERP ref_interp );

REF_STATUS ref_interp_map( REF_INTERP ref_interp, 
			   REF_GRID from, REF_GRID to );

#endif /* REF_INTERP_H */
