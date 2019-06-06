
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

#ifndef REF_ENDIAN_H
#define REF_ENDIAN_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define SWAP_INT(x)          \
  {                          \
    int y;                   \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp + 3) = *(xp + 0);   \
    *(yp + 2) = *(xp + 1);   \
    *(yp + 1) = *(xp + 2);   \
    *(yp + 0) = *(xp + 3);   \
    (x) = y;                 \
  }

#define SWAP_LONG(x)         \
  {                          \
    long y;                  \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp + 7) = *(xp + 0);   \
    *(yp + 6) = *(xp + 1);   \
    *(yp + 5) = *(xp + 2);   \
    *(yp + 4) = *(xp + 3);   \
    *(yp + 3) = *(xp + 4);   \
    *(yp + 2) = *(xp + 5);   \
    *(yp + 1) = *(xp + 6);   \
    *(yp + 0) = *(xp + 7);   \
    (x) = y;                 \
  }

#define SWAP_DBL(x)          \
  {                          \
    double y;                \
    char *xp = (char *)&(x); \
    char *yp = (char *)&(y); \
    *(yp + 7) = *(xp + 0);   \
    *(yp + 6) = *(xp + 1);   \
    *(yp + 5) = *(xp + 2);   \
    *(yp + 4) = *(xp + 3);   \
    *(yp + 3) = *(xp + 4);   \
    *(yp + 2) = *(xp + 5);   \
    *(yp + 1) = *(xp + 6);   \
    *(yp + 0) = *(xp + 7);   \
    (x) = y;                 \
  }

END_C_DECLORATION

#endif /* REF_ENDIAN_H */
