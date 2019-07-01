
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_math.h"

int main(void) {
  { /* same 1 */
    REF_DBL vect1[3] = {1.0, 0.0, 0.0};
    REF_DBL vect2[3] = {1.0, 0.0, 0.0};
    RWDS(1.0, ref_math_dot(vect1, vect2), -1.0, "dot");
  }

  { /* orth 1 */
    REF_DBL vect1[3] = {1.0, 0.0, 0.0};
    REF_DBL vect2[3] = {0.0, 1.0, 0.0};
    RWDS(0.0, ref_math_dot(vect1, vect2), -1.0, "dot");
  }

  { /* two */
    REF_DBL vect1[3] = {1.0, 0.0, 0.0};
    REF_DBL vect2[3] = {2.0, 1.0, 0.0};
    RWDS(2.0, ref_math_dot(vect1, vect2), -1.0, "dot");
  }

  { /* normal zero */
    REF_DBL vect[3] = {0.0, 0.0, 0.0};
    REIS(REF_DIV_ZERO, ref_math_normalize(vect), "expect fail");
    RWDS(0.0, vect[0], -1.0, "same");
    RWDS(0.0, vect[1], -1.0, "same");
    RWDS(0.0, vect[2], -1.0, "same");
  }

  { /* normal one */
    REF_DBL vect[3] = {1.0, 0.0, 0.0};
    RSS(ref_math_normalize(vect), "expect success");
    RWDS(1.0, vect[0], -1.0, "same");
    RWDS(0.0, vect[1], -1.0, "same");
    RWDS(0.0, vect[2], -1.0, "same");
  }

  { /* normal two */
    REF_DBL vect[3] = {2.0, 0.0, 0.0};
    RSS(ref_math_normalize(vect), "expect success");
    RWDS(1.0, vect[0], -1.0, "same");
    RWDS(0.0, vect[1], -1.0, "same");
    RWDS(0.0, vect[2], -1.0, "same");
  }

  { /* orthonormal system x */
    REF_DBL orth0[3] = {0.1, 0.0, 0.0};
    REF_DBL orth1[3], orth2[3];
    RSS(ref_math_orthonormal_system(orth0, orth1, orth2), "ortho");
    RWDS(1.0, orth0[0], -1.0, "same");
    RWDS(0.0, orth0[1], -1.0, "same");
    RWDS(0.0, orth0[2], -1.0, "same");
    RWDS(0.0, orth1[0], -1.0, "o1");
    RWDS(1.0, orth1[1], -1.0, "o1");
    RWDS(0.0, orth1[2], -1.0, "o1");
    RWDS(0.0, orth2[0], -1.0, "o2");
    RWDS(0.0, orth2[1], -1.0, "o2");
    RWDS(1.0, orth2[2], -1.0, "o2");
  }

  { /* orthonormal system y */
    REF_DBL orth0[3] = {0.0, 1.0, 0.0};
    REF_DBL orth1[3], orth2[3];
    RSS(ref_math_orthonormal_system(orth0, orth1, orth2), "ortho");
    RWDS(0.0, orth0[0], -1.0, "same");
    RWDS(1.0, orth0[1], -1.0, "same");
    RWDS(0.0, orth0[2], -1.0, "same");
    RWDS(0.0, orth1[0], -1.0, "o1");
    RWDS(0.0, orth1[1], -1.0, "o1");
    RWDS(1.0, orth1[2], -1.0, "o1");
    RWDS(1.0, orth2[0], -1.0, "o2");
    RWDS(0.0, orth2[1], -1.0, "o2");
    RWDS(0.0, orth2[2], -1.0, "o2");
  }

  { /* orthonormal system z */
    REF_DBL orth0[3] = {0.0, 0.0, 10.0};
    REF_DBL orth1[3], orth2[3];
    RSS(ref_math_orthonormal_system(orth0, orth1, orth2), "ortho");
    RWDS(0.0, orth0[0], -1.0, "same");
    RWDS(0.0, orth0[1], -1.0, "same");
    RWDS(1.0, orth0[2], -1.0, "same");
    RWDS(1.0, orth1[0], -1.0, "o1");
    RWDS(0.0, orth1[1], -1.0, "o1");
    RWDS(0.0, orth1[2], -1.0, "o1");
    RWDS(0.0, orth2[0], -1.0, "o2");
    RWDS(1.0, orth2[1], -1.0, "o2");
    RWDS(0.0, orth2[2], -1.0, "o2");
  }

  { /* acos */
    RWDS(0.0, acos(1.0), -1.0, "acos");
    RWDS(1.04719755119660, acos(0.5), -1.0, "acos");
    RWDS(1.57079632679490, acos(0.0), -1.0, "acos");
    RWDS(2.09439510239320, acos(-0.5), -1.0, "acos");
    RWDS(3.14159265358979, acos(-1.0), -1.0, "acos");
  }

  { /* divisible */
    REIS(REF_FALSE, ref_math_divisible(1.0, 0.0), "1/0");
    REIS(REF_TRUE, ref_math_divisible(1.0, 1.0), "1/1");
    REIS(REF_TRUE, ref_math_divisible(1.0, 1.0e-15), "1/1e-15");
    REIS(REF_FALSE, ref_math_divisible(1.0e40, 1.0e-40), "1e20/1e-20");
  }

  { /* addable overflow protection */
    REIS(REF_TRUE, ref_math_int_addable(0, 0), "0 + 0");
    REIS(REF_TRUE, ref_math_int_addable(2, 5), "2 + 5");
    REIS(REF_TRUE, ref_math_int_addable(-10, -20), "-10 + -20");
    REIS(REF_TRUE, ref_math_int_addable(10, -20), "10 + -20");
    REIS(REF_FALSE, ref_math_int_addable(2147483640, 2147483641),
         "add two INT_MAX - eps");
    REIS(REF_FALSE, ref_math_int_addable(-2147483640, -2147483641),
         "add two -(INT_MAX - eps)");
  }

  { /* multipliable overflow */
    REIS(REF_TRUE, ref_math_int_multipliable(0, 0), "0 * 0");
    REIS(REF_TRUE, ref_math_int_multipliable(3, 7), "3 * 7");
    REIS(REF_TRUE, ref_math_int_multipliable(5, -10), "5 * -10");
    REIS(REF_TRUE, ref_math_int_multipliable(-5, -10), "-5 * -10");
    REIS(REF_FALSE, ref_math_int_multipliable(47000, 47000),
         "multiply two ++ eps + sqrt(INT_MAX)");
    REIS(REF_FALSE, ref_math_int_multipliable(47000, -47000),
         "multiply two +- eps + sqrt(INT_MAX)");
    REIS(REF_FALSE, ref_math_int_multipliable(-47000, -47000),
         "multiply two -- eps + sqrt(INT_MAX)");
  }

  return 0;
}
