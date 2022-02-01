
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

REF_STATUS ref_args_inspect(REF_INT n, char **args) {
  REF_INT i;

  for (i = 0; i < n; i++) printf("%d : '%s'\n", i, args[i]);

  return REF_SUCCESS;
}

REF_STATUS ref_args_find(REF_INT n, char **args, const char *target,
                         REF_INT *pos) {
  REF_INT i;

  *pos = REF_EMPTY;

  for (i = 0; i < n; i++) {
    if (0 == strcmp(target, args[i])) {
      *pos = i;
      return REF_SUCCESS;
    }
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_args_char(REF_INT n, char **args, const char *long_target,
                         const char *short_target, char **value) {
  REF_INT pos;
  *value = NULL;
  if (REF_SUCCESS == ref_args_find(n, args, long_target, &pos)) {
    RAB(pos < n - 1, "missing value",
        { printf("for option %s", long_target); });
    *value = args[pos + 1];
    return REF_SUCCESS;
  }
  if (REF_SUCCESS == ref_args_find(n, args, short_target, &pos)) {
    RAB(pos < n - 1, "missing value",
        { printf("for option %s", short_target); });
    *value = args[pos + 1];
    return REF_SUCCESS;
  }
  return REF_NOT_FOUND;
}
