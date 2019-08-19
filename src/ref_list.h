
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

#ifndef REF_LIST_H
#define REF_LIST_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_LIST_STRUCT REF_LIST_STRUCT;
typedef REF_LIST_STRUCT *REF_LIST;
END_C_DECLORATION

BEGIN_C_DECLORATION
struct REF_LIST_STRUCT {
  REF_INT n, max;
  REF_INT *value;
};

REF_STATUS ref_list_create(REF_LIST *ref_list);
REF_STATUS ref_list_free(REF_LIST ref_list);

REF_STATUS ref_list_deep_copy(REF_LIST *ref_list, REF_LIST original);

REF_STATUS ref_list_inspect(REF_LIST ref_list);

#define ref_list_n(ref_list) ((ref_list)->n)
#define ref_list_max(ref_list) ((ref_list)->max)
#define ref_list_value(ref_list, i) ((ref_list)->value[(i)])

#define each_ref_list_item(ref_list, item) \
  for ((item) = 0; (item) < ref_list_n(ref_list); (item)++)

REF_STATUS ref_list_push(REF_LIST ref_list, REF_INT last);
REF_STATUS ref_list_pop(REF_LIST ref_list, REF_INT *last);
REF_STATUS ref_list_shift(REF_LIST ref_list, REF_INT *first);
REF_STATUS ref_list_delete(REF_LIST ref_list, REF_INT value);

REF_STATUS ref_list_erase(REF_LIST ref_list);

REF_STATUS ref_list_contains(REF_LIST ref_list, REF_INT item,
                             REF_BOOL *contains);

END_C_DECLORATION

#endif /* REF_LIST_H */
