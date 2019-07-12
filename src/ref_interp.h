

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

#ifndef REF_INTERP_H
#define REF_INTERP_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_INTERP_STRUCT REF_INTERP_STRUCT;
typedef REF_INTERP_STRUCT *REF_INTERP;
END_C_DECLORATION

#include "ref_agents.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_mpi.h"
#include "ref_search.h"

BEGIN_C_DECLORATION
struct REF_INTERP_STRUCT {
  REF_MPI ref_mpi;
  REF_GRID from_grid;
  REF_GRID to_grid;
  REF_BOOL instrument;
  REF_BOOL continuously;
  REF_INT n_walk;
  REF_INT n_terminated;
  REF_INT walk_steps;
  REF_INT n_geom;
  REF_INT n_geom_fail;
  REF_INT n_tree;
  REF_INT tree_cells;
  REF_INT max;
  REF_BOOL *agent_hired;
  REF_INT *cell;
  REF_INT *part;
  REF_DBL *bary;
  REF_DBL inside;
  REF_DBL bound;
  REF_AGENTS ref_agents;
  REF_LIST visualize;
  REF_DBL search_fuzz;
  REF_DBL search_donor_scale;
  REF_SEARCH ref_search;
};

#define ref_interp_from_grid(ref_interp) ((ref_interp)->from_grid)
#define ref_interp_to_grid(ref_interp) ((ref_interp)->to_grid)
#define ref_interp_cell(ref_interp, node) ((ref_interp)->cell[(node)])
#define ref_interp_part(ref_interp, node) ((ref_interp)->part[(node)])
#define ref_interp_bary(ref_interp, j, node) \
  ((ref_interp)->bary[(j) + 4 * (node)])
#define ref_interp_max(ref_interp) ((ref_interp)->max)
#define ref_interp_continuously(ref_interp) ((ref_interp)->continuously)
#define ref_interp_search(ref_interp) ((ref_interp)->ref_search)
#define ref_interp_search_fuzz(ref_interp) ((ref_interp)->search_fuzz)
#define ref_interp_search_donor_scale(ref_interp) \
  ((ref_interp)->search_donor_scale)

REF_STATUS ref_interp_create(REF_INTERP *ref_interp, REF_GRID from_grid,
                             REF_GRID to_grid);
REF_STATUS ref_interp_create_identity(REF_INTERP *ref_interp,
                                      REF_GRID ref_grid);

REF_STATUS ref_interp_free(REF_INTERP ref_interp);

REF_STATUS ref_interp_pack(REF_INTERP ref_interp, REF_INT *n2o);

REF_STATUS ref_interp_remove(REF_INTERP ref_interp, REF_INT node);

REF_STATUS ref_interp_tattle(REF_INTERP ref_interp, REF_INT node);
REF_STATUS ref_interp_locate(REF_INTERP ref_interp);
REF_STATUS ref_interp_locate_subset(REF_INTERP ref_interp);
REF_STATUS ref_interp_locate_nearest(REF_INTERP ref_interp);
REF_STATUS ref_interp_locate_node(REF_INTERP ref_interp, REF_INT node);
REF_STATUS ref_interp_locate_between(REF_INTERP ref_interp, REF_INT node0,
                                     REF_INT node1, REF_INT new_node);

REF_STATUS ref_interp_scalar(REF_INTERP ref_interp, REF_INT leading_dim,
                             REF_DBL *from_scalar, REF_DBL *to_scalar);

REF_STATUS ref_interp_min_bary(REF_INTERP ref_interp, REF_DBL *min_bary);
REF_STATUS ref_interp_max_error(REF_INTERP ref_interp, REF_DBL *max_error);
REF_STATUS ref_interp_stats(REF_INTERP ref_interp);

REF_STATUS ref_interp_tec(REF_INTERP ref_interp, const char *filename);
REF_STATUS ref_interp_integrate(REF_GRID ref_grid, REF_DBL *canidate,
                                REF_DBL *truth, REF_INT norm_power,
                                REF_DBL *error);

REF_STATUS ref_interp_convergence_rate(REF_DBL f3, REF_DBL h3, REF_DBL f2,
                                       REF_DBL h2, REF_DBL f1, REF_DBL h1,
                                       REF_DBL *rate);

REF_STATUS ref_iterp_plt(REF_GRID ref_grid, const char *filename, REF_INT *ldim,
                         REF_DBL **scalar);

END_C_DECLORATION

#endif /* REF_INTERP_H */
