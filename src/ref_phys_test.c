
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

#include <stdio.h>
#include <stdlib.h>

#include "ref_phys.h"

int main(void) {
  { /* x-Euler flux */
    REF_DBL state[5], direction[3];
    REF_DBL flux[5];
    state[0] = 1.0;
    state[1] = 0.2;
    state[2] = 0.0;
    state[3] = 0.0;
    state[4] = 1.0 / 1.4;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    RSS(ref_phys_euler(state, direction, flux), "euler");
    RWDS(0.2, flux[0], -1, "mass flux");
    RWDS(0.04 + 1.0 / 1.4, flux[1], -1, "x mo flux");
    RWDS(0.0, flux[2], -1, "y mo flux");
    RWDS(0.0, flux[3], -1, "z mo flux");
    RWDS(0.504, flux[4], -1, "energy flux");
  }

  { /* Couette laminar flux */
    REF_DBL state[5], gradient[15], direction[3];
    REF_DBL flux[5];
    REF_DBL mach = 0.1, re = 10.0, temp = 273.0;
    REF_DBL dudy = 1.0, mu = 1.0;
    REF_INT i;
    for (i = 0; i < 15; i++) gradient[i] = 0.0;
    gradient[1 + 3 * 1] = dudy;
    state[0] = 1.0;
    state[1] = 0.1;
    state[2] = 0.0;
    state[3] = 0.0;
    state[4] = 1.0 / 1.4;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    RSS(ref_phys_laminar(state, gradient, mach, re, temp, direction, flux),
        "euler");
    RWDS(0.0, flux[0], -1, "mass flux");
    RWDS(0.0, flux[1], -1, "x mo flux");
    RWDS(mach/re*mu*dudy*direction[0], flux[2], -1, "y mo flux");
    RWDS(0.0, flux[3], -1, "z mo flux");
  }

  return 0;
}
