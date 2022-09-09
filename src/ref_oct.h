
/* Copyright 2006, 2014, 2021 United States Government as represented
 * by the Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine version 3 unstructured grid adaptation platform is
 * licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#ifndef REF_OCT_H
#define REF_OCT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_OCT_STRUCT REF_OCT_STRUCT;
typedef REF_OCT_STRUCT *REF_OCT;
END_C_DECLORATION

BEGIN_C_DECLORATION
struct REF_OCT_STRUCT {
  REF_DBL bbox[6];
  REF_INT n, max;
  REF_INT *children;
};

REF_FCN REF_STATUS ref_oct_create(REF_OCT *ref_oct);

REF_FCN REF_STATUS ref_oct_free(REF_OCT ref_oct);

REF_FCN REF_STATUS ref_oct_split(REF_OCT ref_oct, REF_INT node);

REF_FCN REF_STATUS ref_oct_contains(REF_OCT ref_oct, REF_DBL *xyz,
                                    REF_INT *node);

REF_FCN REF_STATUS ref_oct_tec(REF_OCT ref_oct, const char *filename);

END_C_DECLORATION

#endif /* REF_OCT_H */
