
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

#ifndef REF_FIXTURE_H
#define REF_FIXTURE_H

#include "ref_defs.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_fixture_tri_surf_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_tet_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_tet2_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_pyr_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_tri_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_tri2_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_tri_qua_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_twod_cubic_edge(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_pri_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_pri_tet_cap_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_pri_stack_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_hex_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);

REF_STATUS ref_fixture_hanging_hex_pri_grid(REF_GRID *ref_grid,
                                            REF_MPI ref_mpi);

REF_STATUS ref_fixture_hex_brick_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_hex_brick_args_grid(REF_GRID *ref_grid, REF_MPI ref_mpi,
                                           REF_DBL x0, REF_DBL x1, REF_DBL y0,
                                           REF_DBL y1, REF_DBL z0, REF_DBL z1,
                                           REF_INT l, REF_INT m, REF_INT n);

REF_STATUS ref_fixture_tet_brick_grid(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_fixture_tet_brick_args_grid(REF_GRID *ref_grid, REF_MPI ref_mpi,
                                           REF_DBL x0, REF_DBL x1, REF_DBL y0,
                                           REF_DBL y1, REF_DBL z0, REF_DBL z1,
                                           REF_INT l, REF_INT m, REF_INT n);
REF_STATUS ref_fixture_twod_brick_grid(REF_GRID *ref_grid, REF_MPI ref_mpi,
                                       REF_INT dim);

REF_STATUS ref_fixture_twod_square_circle(REF_GRID *ref_grid, REF_MPI ref_mpi);

END_C_DECLORATION

#endif /* REF_FIXTURE_H */
