
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

#include "ref_args.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  if (argc > 1) {
    RSS(ref_args_inspect(argc, argv), "echo");
  }

  {
    REF_INT n = 3;
    char *a0 = "program";
    char *a1 = "-1";
    char *a2 = "-2";
    char *as[3];
    REF_INT pos;
    as[0] = a0;
    as[1] = a1;
    as[2] = a2;

    RSS(ref_args_find(n, as, "-1", &pos), "echo");
    REIS(1, pos, "location");

    REIS(REF_NOT_FOUND, ref_args_find(n, as, "-h", &pos), "not found");
    REIS(REF_EMPTY, pos, "location");

    REIS(REF_NOT_FOUND, ref_args_find(n, as, "--long", &pos), "not found");
    REIS(REF_EMPTY, pos, "location");
  }

  return 0;
}
