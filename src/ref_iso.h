
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

#ifndef REF_ISO_H
#define REF_ISO_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_iso_insert(REF_GRID *iso_grid, REF_GRID ref_grid, REF_DBL *field,
                          REF_INT ldim, REF_DBL *in, REF_DBL **out);

REF_STATUS ref_iso_signed_distance(REF_GRID ref_grid, REF_DBL *field,
                                   REF_DBL *distance);

REF_STATUS ref_iso_triangle_segment(REF_DBL *triangle0, REF_DBL *triangle1,
                                    REF_DBL *triangle2, REF_DBL *segment0,
                                    REF_DBL *segment1, REF_DBL *tuvw);
REF_STATUS ref_iso_segment_segment(REF_DBL *candidate0, REF_DBL *candidate1,
                                   REF_DBL *segment0, REF_DBL *segment1,
                                   REF_DBL *tt);

REF_STATUS ref_iso_cast(REF_GRID *iso_grid, REF_DBL **iso_field,
                        REF_GRID ref_grid, REF_DBL *field, REF_INT ldim,
                        REF_DBL *segment0, REF_DBL *segment1);

REF_STATUS ref_iso_segment(REF_GRID ref_grid, REF_DBL *center, REF_DBL aoa,
                           REF_DBL phi, REF_DBL h, REF_DBL *segment0,
                           REF_DBL *segment1);

REF_STATUS ref_iso_boom_header(FILE **file, REF_INT ldim,
                               const char **scalar_names, const char *filename);
REF_STATUS ref_iso_boom_zone(FILE *file, REF_GRID ref_grid, REF_DBL *field,
                             REF_INT ldim, REF_DBL *center, REF_DBL aoa,
                             REF_DBL phi, REF_DBL h);

REF_STATUS ref_iso_slice(REF_GRID *iso_grid, REF_GRID ref_grid, REF_DBL *normal,
                         REF_DBL offset, REF_INT ldim, REF_DBL *in,
                         REF_DBL **out);

END_C_DECLORATION

#endif /* REF_ISO_H */
