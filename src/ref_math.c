
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

#include "ref_math.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

REF_STATUS ref_math_normalize(REF_DBL *normal) {
  REF_DBL length;

  length = sqrt(ref_math_dot(normal, normal));

  if (!ref_math_divisible(normal[0], length) ||
      !ref_math_divisible(normal[1], length) ||
      !ref_math_divisible(normal[2], length))
    return REF_DIV_ZERO;

  normal[0] /= length;
  normal[1] /= length;
  normal[2] /= length;

  length = ref_math_dot(normal, normal);
  RAS((ABS(length - 1.0) < 1.0e-13), "vector length not unity");

  return REF_SUCCESS;
}

REF_STATUS ref_math_orthonormal_system(REF_DBL *orth0, REF_DBL *orth1,
                                       REF_DBL *orth2) {
  REF_DBL dot;
  RSS(ref_math_normalize(orth0), "orthonormalize input orth0");
  if (ABS(orth0[0]) >= ABS(orth0[1]) && ABS(orth0[0]) >= ABS(orth0[2])) {
    orth1[0] = 0.0;
    orth1[1] = 1.0;
    orth1[2] = 0.0;
  }
  if (ABS(orth0[1]) >= ABS(orth0[0]) && ABS(orth0[1]) >= ABS(orth0[2])) {
    orth1[0] = 0.0;
    orth1[1] = 0.0;
    orth1[2] = 1.0;
  }
  if (ABS(orth0[2]) >= ABS(orth0[0]) && ABS(orth0[2]) >= ABS(orth0[1])) {
    orth1[0] = 1.0;
    orth1[1] = 0.0;
    orth1[2] = 0.0;
  }
  dot = ref_math_dot(orth0, orth1);
  orth1[0] -= dot * orth0[0];
  orth1[1] -= dot * orth0[1];
  orth1[2] -= dot * orth0[2];
  RSS(ref_math_normalize(orth1), "orthonormalize orth1");
  ref_math_cross_product(orth0, orth1, orth2);
  return REF_SUCCESS;
}
