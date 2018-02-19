
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

#ifndef REF_ELAST_H
#define REF_ELAST_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ELAST_STRUCT REF_ELAST_STRUCT;
typedef REF_ELAST_STRUCT * REF_ELAST;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_comprow.h"

BEGIN_C_DECLORATION

struct REF_ELAST_STRUCT {
  REF_GRID ref_grid;
  REF_COMPROW ref_comprow;
  REF_DBL *a;
  REF_DBL *displacement;
  REF_INT *bc;
};

#define ref_elast_grid(ref_elast) ((ref_elast)->ref_grid)
#define ref_elast_comprow(ref_elast) ((ref_elast)->ref_comprow)

REF_STATUS ref_elast_create( REF_ELAST *ref_elast, REF_GRID ref_grid );
REF_STATUS ref_elast_free( REF_ELAST ref_elast );

REF_STATUS ref_elast_inspect( REF_ELAST ref_elast );

REF_STATUS ref_elast_displace( REF_ELAST ref_elast,
                               REF_INT node, REF_DBL *dxyz );
REF_STATUS ref_elast_assemble( REF_ELAST ref_elast );

REF_STATUS ref_elast_relax( REF_ELAST ref_elast, REF_DBL *l2norm );

END_C_DECLORATION

#endif /* REF_ELAST_H */
