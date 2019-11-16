
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

#include "ref_html.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  if (2 == argc) {
    REF_HTML ref_html;
    REF_DBL origin[3];
    REF_DBL sys[12];
    RSS(ref_html_create(&ref_html, argv[1]), "open html");

    origin[0] = 0.0;
    origin[1] = 0.0;
    origin[2] = 0.0;
    sys[0] = 0.25;
    sys[1] = 0.5;
    sys[2] = 1.0;
    sys[3] = 1.0;
    sys[4] = 0.0;
    sys[5] = 0.0;
    sys[6] = 0.0;
    sys[7] = 1.0;
    sys[8] = 0.0;
    sys[9] = 0.0;
    sys[10] = 0.0;
    sys[11] = 1.0;
    RSS(ref_html_diagonal_system(ref_html, origin, sys), "sys html");

    RSS(ref_html_free(ref_html), "close html");
    return 0;
  }

  { /* html .meshb */
    REF_HTML ref_html;
    char file[] = "ref_html_test.html";
    RSS(ref_html_create(&ref_html, file), "open html");
    RSS(ref_html_free(ref_html), "close html");
    REIS(0, remove(file), "test clean up");
  }

  return 0;
}
