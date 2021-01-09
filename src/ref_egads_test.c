
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

#include "ref_egads.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"
#include "ref_validation.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--recon", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY) {
    REF_GRID ref_grid;
    REF_INT node;
    REIS(4, argc, "required args: --recon grid.ext geom.egads");
    REIS(1, pos, "required args: --recon grid.ext geom.egads");
    printf("reconstruct geometry association\n");
    printf("grid source %s\n", argv[2]);
    printf("geometry source %s\n", argv[3]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "argv import");
    RSS(ref_egads_load(ref_grid_geom(ref_grid), argv[3]), "ld egads");
    RSS(ref_egads_recon(ref_grid), "geom recon");
    printf("verify topo and params\n");
    RSS(ref_geom_verify_topo(ref_grid), "geom topo conflict");
    RSS(ref_geom_verify_param(ref_grid), "test constrained params");
    printf("validate\n");
    RSS(ref_validation_all(ref_grid), "validate");
    printf("constrain\n");
    RSS(ref_export_tec_surf(ref_grid, "ref_geom_orig.tec"), "tec");
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RSS(ref_geom_constrain(ref_grid, node), "original params");
    }
    RSS(ref_geom_tec(ref_grid, "ref_geom_recon.tec"), "geom export");
    printf("validate\n");
    RSS(ref_validation_all(ref_grid), "validate");
    RSS(ref_export_by_extension(ref_grid, "ref_geom_recon.meshb"), "export");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    return 0;
  }

  if (ref_egads_allows_construction()) { /* single cylinder */
    REF_GEOM ref_geom;
    RSS(ref_geom_create(&ref_geom), "create geom");
    RSS(ref_egads_construct(ref_geom, "cylinder"), "create");
    RSS(ref_geom_free(ref_geom), "free geom");
  }

  if (ref_egads_allows_construction()) { /* single cylinder tess */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "cylinder"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* single cylinder enrich2 */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "cylinder"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_enrich2(ref_grid), "enrich2");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* single cylinder enrich3 */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "cylinder"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_enrich3(ref_grid), "enrich3");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* revolve enrich3 */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "revolve"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_enrich3(ref_grid), "enrich3");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* steinmetz */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "steinmetz"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_verify_param(ref_grid), "egads params");
    /* RSS(ref_geom_tec(ref_grid, "steinmetz.tec"),"geom"); */
    /* RSS(ref_export_by_extension(ref_grid, "steinmetz.meshb"), "meshb"); */
    /* RSS(ref_egads_save(ref_grid_geom(ref_grid), "steinmetz.egads"), "egd");*/
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (ref_egads_allows_construction()) { /* revolve */
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid, ref_mpi), "create grid");
    RSS(ref_egads_construct(ref_grid_geom(ref_grid), "revolve"), "create");
    RSS(ref_egads_tess(ref_grid, 0, NULL), "tess");
    RSS(ref_geom_verify_param(ref_grid), "egads params");
    /* RSS(ref_geom_tec(ref_grid, "revolve.tec"), "geom"); */
    /* RSS(ref_export_by_extension(ref_grid, "revolve.meshb"), "meshb"); */
    /* RSS(ref_egads_save(ref_grid_geom(ref_grid), "revolve.egads"), "egd"); */
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
