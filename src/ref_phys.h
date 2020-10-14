
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

#ifndef REF_PHYS_H
#define REF_PHYS_H

#include <stdio.h>

#include "ref_defs.h"
#include "ref_dict.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_phys_make_primitive(REF_DBL *conserved, REF_DBL *primitive);
REF_STATUS ref_phys_make_conserved(REF_DBL *primitive, REF_DBL *conserved);
REF_STATUS ref_phys_entropy_adjoint(REF_DBL *primitive, REF_DBL *dual);

REF_STATUS ref_phys_euler(REF_DBL *state, REF_DBL *direction, REF_DBL *flux);
REF_STATUS ref_phys_euler_jac(REF_DBL *state, REF_DBL *direction,
                              REF_DBL *dflux_dcons);
REF_STATUS ref_phys_viscous(REF_DBL *state, REF_DBL *grad, REF_DBL turb,
                            REF_DBL mach, REF_DBL re, REF_DBL reference_temp,
                            REF_DBL *dir, REF_DBL *flux);
REF_STATUS ref_phys_mut_sa(REF_DBL turb, REF_DBL rho, REF_DBL nu,
                           REF_DBL *mut_sa);
REF_STATUS ref_phys_convdiff(REF_DBL *state, REF_DBL *grad, REF_DBL diffusivity,
                             REF_DBL *dir, REF_DBL *flux);

REF_STATUS ref_phys_euler_dual_flux(REF_GRID ref_grid, REF_INT ldim,
                                    REF_DBL *primitive_dual,
                                    REF_DBL *dual_flux);
REF_STATUS ref_phys_mask_strong_bcs(REF_GRID ref_grid, REF_DICT ref_dict,
                                    REF_BOOL *replace, REF_INT ldim);

REF_STATUS ref_phys_read_mapbc(REF_DICT ref_dict, const char *mapbc_filename);
REF_STATUS ref_phys_parse_tags(REF_DICT ref_dict, const char *tags);

REF_STATUS ref_phys_cc_fv_res(REF_GRID ref_grid, REF_INT nequ, REF_DBL *flux,
                              REF_DBL *res);
REF_STATUS ref_phys_cc_fv_embed(REF_GRID ref_grid, REF_INT nequ, REF_DBL *flux,
                                REF_DBL *res);

REF_STATUS ref_phys_spalding_yplus(REF_DBL uplus, REF_DBL *yplus);
REF_STATUS ref_phys_spalding_dyplus_duplus(REF_DBL uplus,
                                           REF_DBL *dyplus_duplus);
REF_STATUS ref_phys_spalding_uplus(REF_DBL yplus, REF_DBL *uplus);

REF_STATUS ref_phys_signed_distance(REF_GRID ref_grid, REF_DBL *field,
                                    REF_DBL *distance);

REF_STATUS ref_phys_wall_distance(REF_GRID ref_grid, REF_DICT ref_dict,
                                  REF_DBL *distance);

END_C_DECLORATION

#endif /* REF_PHYS_H */
