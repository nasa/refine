
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

#ifndef REF_FORTRAN_H
#define REF_FORTRAN_H

#include "ref_defs.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

BEGIN_C_DECLORATION

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran */
#define REF_FORT(name, NAME) name##_
/* but for C identifiers containing underscores. */
#define REF_FORT_(name, NAME) name##_

extern REF_BOOL ref_fortran_allow_screen_output;

REF_STATUS REF_FORT_(ref_fortran_init,
                     REF_FORTRAN_INIT)(REF_INT *nnodes, REF_GLOB *nnodesg,
                                       REF_GLOB *l2g, REF_INT *part,
                                       REF_INT *partition, REF_DBL *x,
                                       REF_DBL *y, REF_DBL *z);

REF_STATUS REF_FORT_(ref_fortran_import_cell,
                     REF_FORTRAN_IMPORT_CELL)(REF_INT *node_per_cell,
                                              REF_INT *ncell, REF_INT *c2n);

REF_STATUS REF_FORT_(ref_fortran_import_face,
                     REF_FORTRAN_IMPORT_FACE)(REF_INT *face_index,
                                              REF_INT *node_per_face,
                                              REF_INT *nface, REF_INT *f2n);

REF_STATUS REF_FORT_(ref_fortran_import_metric,
                     REF_FORTRAN_IMPORT_METRIC)(REF_INT *nnodes,
                                                REF_DBL *metric);

REF_STATUS REF_FORT_(ref_fortran_import_ratio,
                     REF_FORTRAN_IMPORT_RATIO)(REF_INT *nnodes, REF_DBL *ratio);

REF_STATUS REF_FORT_(ref_fortran_viz, REF_FORTRAN_VIZ)(void);

REF_STATUS REF_FORT_(ref_fortran_adapt, REF_FORTRAN_ADAPT)(void);

REF_STATUS REF_FORT_(ref_fortran_size_node,
                     REF_FORTRAN_SIZE_NODE)(REF_INT *nnodes0, REF_INT *nnodes,
                                            REF_GLOB *nnodesg);

REF_STATUS REF_FORT_(ref_fortran_node,
                     REF_FORTRAN_NODE)(REF_INT *nnodes, REF_GLOB *l2g,
                                       REF_DBL *x, REF_DBL *y, REF_DBL *z);

REF_STATUS REF_FORT_(ref_fortran_size_cell,
                     REF_FORTRAN_SIZE_CELL)(REF_INT *node_per_cell,
                                            REF_INT *ncell);

REF_STATUS REF_FORT_(ref_fortran_cell, REF_FORTRAN_CELL)(REF_INT *node_per_cell,
                                                         REF_INT *ncell,
                                                         REF_INT *c2n);

REF_STATUS REF_FORT_(ref_fortran_size_face,
                     REF_FORTRAN_SIZE_FACE)(REF_INT *ibound,
                                            REF_INT *node_per_face,
                                            REF_INT *nface);

REF_STATUS REF_FORT_(ref_fortran_face,
                     REF_FORTRAN_FACE)(REF_INT *ibound, REF_INT *node_per_face,
                                       REF_INT *nface, REF_INT *f2n);

REF_STATUS REF_FORT_(ref_fortran_naux, REF_FORTRAN_NAUX)(REF_INT *naux);
REF_STATUS REF_FORT_(ref_fortran_import_aux,
                     REF_FORTRAN_IMPORT_AUX)(REF_INT *ldim, REF_INT *nnodes,
                                             REF_INT *offset, REF_DBL *aux);
REF_STATUS REF_FORT_(ref_fortran_aux,
                     REF_FORTRAN_AUX)(REF_INT *ldim, REF_INT *nnodes,
                                      REF_INT *offset, REF_DBL *aux);

REF_STATUS REF_FORT_(ref_fortran_free, REF_FORTRAN_FREE)(void);

END_C_DECLORATION

#endif /* REF_FORTRAN_H */
