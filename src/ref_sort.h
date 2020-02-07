
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

#ifndef REF_SORT_H
#define REF_SORT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

REF_STATUS ref_sort_insertion_int(REF_INT n, REF_INT *original,
                                  REF_INT *sorted);

REF_STATUS ref_sort_heap_int(REF_INT n, REF_INT *original,
                             REF_INT *sorted_index);
REF_STATUS ref_sort_heap_glob(REF_INT n, REF_GLOB *original,
                              REF_INT *sorted_index);
REF_STATUS ref_sort_heap_dbl(REF_INT n, REF_DBL *original,
                             REF_INT *sorted_index);

REF_STATUS ref_sort_unique_int(REF_INT n, REF_INT *original, REF_INT *nunique,
                               REF_INT *unique);

REF_STATUS ref_sort_search_int(REF_INT n, REF_INT *ascending_list,
                               REF_INT target, REF_INT *position);
REF_STATUS ref_sort_search_glob(REF_INT n, REF_GLOB *ascending_list,
                                REF_GLOB target, REF_INT *position);

END_C_DECLORATION

#endif /* REF_SORT_H */
