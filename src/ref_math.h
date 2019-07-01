
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

#ifndef REF_MATH_H
#define REF_MATH_H

#include <math.h>

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_math_dot(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])
#define ref_math_pi (3.14159265358979)
#define ref_math_in_degrees(radians) ((radians)*180.0 / 3.14159265358979)
#define ref_math_in_radians(degrees) ((degrees) / 180.0 * 3.14159265358979)

/* tester : printf(" %e / %e = %e \n",n,d,n/d) */
#define ref_math_divisible(n, d) (ABS(1.0e20 * d) > ABS(n))

#define ref_math_int_addable(a, b)                   \
  ((long)(REF_INT_MAX) >= ((long)(a) + (long)(b)) && \
   ((long)(a) + (long)(b)) >= (long)(REF_INT_MIN))

#define ref_math_int_multipliable(a, b)              \
  ((long)(REF_INT_MAX) >= ((long)(a) * (long)(b)) && \
   ((long)(a) * (long)(b)) >= (long)(REF_INT_MIN))

#define ref_math_cross_product(v0, v1, product)         \
  (product)[0] = (v0)[1] * (v1)[2] - (v0)[2] * (v1)[1]; \
  (product)[1] = (v0)[2] * (v1)[0] - (v0)[0] * (v1)[2]; \
  (product)[2] = (v0)[0] * (v1)[1] - (v0)[1] * (v1)[0];

REF_STATUS ref_math_normalize(REF_DBL *normal);

REF_STATUS ref_math_orthonormal_system(REF_DBL *orth0, REF_DBL *orth1,
                                       REF_DBL *orth2);

END_C_DECLORATION

#endif /* REF_MATH_H */
