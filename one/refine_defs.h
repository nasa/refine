
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

#ifndef MASTER_HEADER_H
#define MASTER_HEADER_H

#ifdef __cplusplus
#  define BEGIN_C_DECLORATION extern "C" {
#  define END_C_DECLORATION }
#else
#  define BEGIN_C_DECLORATION
#  define END_C_DECLORATION
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

BEGIN_C_DECLORATION

#define EMPTY (-1)

#if !defined(ABS)
#define ABS(a)   ((a)>0?(a):-(a))
#endif

#if !defined(SIGN)
#define SIGN(a)   ( (a)==0 ? 0 : ((a)>0?1:-1) )
#endif

typedef int REF_STATUS;

#define REF_SUCCESS       (0)
#define REF_FAILURE       (1)

#define TSS(fcn,msg)							\
  {									\
    REF_STATUS code;							\
    code = (fcn);							\
    if (REF_SUCCESS != code){						\
      printf("%s: %d: %s: %d %s\n",__FILE__,__LINE__,__func__,code,(msg)); \
      return code;							\
    }									\
  }

#define SUPRESS_UNUSED_COMPILER_WARNING(ptr)                    \
  if (NULL == (&(ptr)+1)) printf("unused macro failed\n");

END_C_DECLORATION

#ifdef HAVE_SDK
#include <MeatLib/Common.h>
typedef GeoBool GridBool;
#else
/* lifted defs from the SDK/MeatLib/Common.h */

BEGIN_C_DECLORATION

typedef short   GridBool;

#undef TRUE
#undef FALSE
#define TRUE    ((GridBool)1)
#define FALSE   ((GridBool)0)

#if !defined(MIN)
#define MIN(a,b) ((a)<(b)?(a):(b)) 
#endif
#if !defined(MAX)
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

END_C_DECLORATION

#endif /* HAVE_SDK */
#endif /* MASTER_HEADER_H */
