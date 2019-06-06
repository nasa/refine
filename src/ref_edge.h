
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

#ifndef REF_EDGE_H
#define REF_EDGE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_EDGE_STRUCT REF_EDGE_STRUCT;
typedef REF_EDGE_STRUCT *REF_EDGE;
END_C_DECLORATION

#include "ref_adj.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_EDGE_STRUCT {
  REF_INT n, max;
  REF_INT *e2n;
  REF_ADJ adj;
  REF_NODE node;
};

REF_STATUS ref_edge_create(REF_EDGE *ref_edge, REF_GRID ref_grid);
REF_STATUS ref_edge_free(REF_EDGE ref_edge);

#define ref_edge_n(ref_edge) ((ref_edge)->n)
#define ref_edge_max(ref_edge) ((ref_edge)->max)

#define ref_edge_e2n(ref_edge, node, edge) ((ref_edge)->e2n[node + 2 * edge])

#define ref_edge_adj(ref_edge) ((ref_edge)->adj)
#define ref_edge_node(ref_edge) ((ref_edge)->node)

#define each_ref_edge(ref_edge, edge) \
  for ((edge) = 0; (edge) < ref_edge_n(ref_edge); (edge)++)

#define each_edge_having_node(ref_edge, node, item, edge) \
  each_ref_adj_node_item_with_ref(ref_edge_adj(ref_edge), node, item, edge)

REF_STATUS ref_edge_uniq(REF_EDGE ref_edge, REF_INT node0, REF_INT node1);

REF_STATUS ref_edge_with(REF_EDGE ref_edge, REF_INT node0, REF_INT node1,
                         REF_INT *edge);

REF_STATUS ref_edge_part(REF_EDGE ref_edge, REF_INT edge, REF_INT *part);

REF_STATUS ref_edge_ghost_min_int(REF_EDGE ref_edge, REF_MPI ref_mpi,
                                  REF_INT *data);
REF_STATUS ref_edge_ghost_int(REF_EDGE ref_edge, REF_MPI ref_mpi,
                              REF_INT *data);
REF_STATUS ref_edge_ghost_glob(REF_EDGE ref_edge, REF_MPI ref_mpi,
                               REF_GLOB *data);
REF_STATUS ref_edge_ghost_dbl(REF_EDGE ref_edge, REF_MPI ref_mpi, REF_DBL *data,
                              REF_INT dim);

REF_STATUS ref_edge_tec_fill(REF_EDGE ref_edge, const char *filename);
REF_STATUS ref_edge_tec_int(REF_EDGE ref_edge, const char *filename,
                            REF_INT *data);
REF_STATUS ref_edge_tec_dbl(REF_EDGE ref_edge, const char *filename,
                            REF_DBL *data);

REF_STATUS ref_edge_tec_ratio(REF_EDGE ref_edge, const char *root_filename);

REF_STATUS ref_edge_rcm(REF_EDGE ref_edge, REF_INT **o2n, REF_INT **n2o);

END_C_DECLORATION

#endif /* REF_EDGE_H */
