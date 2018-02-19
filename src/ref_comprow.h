
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

#ifndef REF_COMPROW_H
#define REF_COMPROW_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_COMPROW_STRUCT REF_COMPROW_STRUCT;
typedef REF_COMPROW_STRUCT * REF_COMPROW;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_COMPROW_STRUCT {
  REF_INT max, nnz;
  REF_INT *first;
  REF_INT *col;
};

#define ref_comprow_max(ref_comprow) ((ref_comprow)->max)
#define ref_comprow_nnz(ref_comprow) ((ref_comprow)->nnz)

#define each_ref_comprow_row_entry( ref_comprow, row, entry )           \
  for ( (entry) = ref_comprow->first[(row)];                            \
        (entry) < ref_comprow->first[(row)+1];                          \
        (entry)++ )                                                     \


REF_STATUS ref_comprow_create( REF_COMPROW *ref_comprow, REF_GRID ref_grid );
REF_STATUS ref_comprow_free( REF_COMPROW ref_comprow );

REF_STATUS ref_comprow_inspect( REF_COMPROW ref_comprow );

REF_STATUS ref_comprow_entry( REF_COMPROW ref_comprow,
                              REF_INT row, REF_INT col, REF_INT *entry );

END_C_DECLORATION

#endif /* REF_COMPROW_H */
