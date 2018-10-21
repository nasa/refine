
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

#ifndef UNRAVEL_FORTRAN_H
#define UNRAVEL_FORTRAN_H

#include "refine_defs.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

BEGIN_C_DECLORATION
void REF_FORT_(unravel_start,UNRAVEL_START)( int *unravel_api_version, int *status );
void REF_FORT_(unravel_tet,UNRAVEL_TET)( int *c2n, double *x, double *y, double *z, int *status );
void REF_FORT_(unravel_thaw,UNRAVEL_THAW)( int *nodeid, int *status );
void REF_FORT_(unravel_it,UNRAVEL_IT)( int *status );
void REF_FORT_(unravel_xyz,UNRAVEL_XYZ)( int *nodeid, double *x, double *y, double *z, int *status );
void REF_FORT_(unravel_cleanup,UNRAVEL_CLEANUP)( int *status );
END_C_DECLORATION

#endif /* UNRAVEL_FORTRAN_H */

