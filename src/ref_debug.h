
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

#ifndef REF_DEBUG_H
#define REF_DEBUG_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_debug_print_int(i) printf(" " #i " = %d\n", i);
#define ref_debug_print_dbl(d) printf(" " #d " = %e\n", i);
#define ref_debug_print_str(s) printf(" " #s "\n");

#define ref_debug_print_int_array(a, n)                                      \
  {                                                                          \
    REF_INT ref_debug_print_int_array_i;                                     \
    for (ref_debug_print_int_array_i = 0; ref_debug_print_int_array_i < (n); \
         ref_debug_print_int_array_i++)                                      \
      printf(" " #a "[%d] = %d\n", ref_debug_print_int_array_i,              \
             (a)[ref_debug_print_int_array_i]);                              \
  }

#define ref_debug_print_dbl_array(a, n)                                      \
  {                                                                          \
    REF_INT ref_debug_print_dbl_array_i;                                     \
    for (ref_debug_print_dbl_array_i = 0; ref_debug_print_dbl_array_i < (n); \
         ref_debug_print_dbl_array_i++)                                      \
      printf(" " #a "[%d] = %e\n", ref_debug_print_dbl_array_i,              \
             (a)[ref_debug_print_dbl_array_i]);                              \
  }

#define ref_debug_print_dbl2_array(a, l, n)                                    \
  {                                                                            \
    REF_INT ref_debug_print_int_array_i;                                       \
    REF_INT ref_debug_print_int_array_j;                                       \
    for (ref_debug_print_int_array_j = 0; ref_debug_print_int_array_j < (n);   \
         ref_debug_print_int_array_j++) {                                      \
      printf(" " #a "[%d] =", ref_debug_print_int_array_j);                    \
      for (ref_debug_print_int_array_i = 0; ref_debug_print_int_array_i < (l); \
           ref_debug_print_int_array_i++)                                      \
        printf(" %f", (a)[ref_debug_print_int_array_i +                        \
                          (l)*ref_debug_print_int_array_j]);                   \
      printf("\n");                                                            \
    }                                                                          \
  }

END_C_DECLORATION

#endif /* REF_DEBUG_H */
