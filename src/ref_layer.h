
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

#ifndef REF_LAYER_H
#define REF_LAYER_H

#include "ref_defs.h"
#include "ref_dict.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION
typedef struct REF_LAYER_STRUCT REF_LAYER_STRUCT;
typedef REF_LAYER_STRUCT *REF_LAYER;
END_C_DECLORATION

#include "ref_cell.h"
#include "ref_list.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

struct REF_LAYER_STRUCT {
  REF_LIST ref_list;
  REF_GRID ref_grid;
  REF_INT nnode_per_layer;
  REF_BOOL verbose;
};

REF_STATUS ref_layer_create(REF_LAYER *ref_layer, REF_MPI ref_mpi);
REF_STATUS ref_layer_free(REF_LAYER ref_layer);

#define ref_layer_list(ref_layer) ((ref_layer)->ref_list)
#define ref_layer_n(ref_layer) (ref_list_n(ref_layer_list(ref_layer)))
#define ref_layer_grid(ref_layer) ((ref_layer)->ref_grid)

REF_STATUS ref_layer_attach(REF_LAYER ref_layer, REF_GRID ref_grid,
                            REF_INT faceid);
REF_STATUS ref_layer_puff(REF_LAYER ref_layer, REF_GRID ref_grid);
REF_STATUS ref_layer_insert(REF_LAYER ref_layer, REF_GRID ref_grid);
REF_STATUS ref_layer_recon(REF_LAYER ref_layer, REF_GRID ref_grid);

REF_STATUS ref_layer_align_quad(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_LAYER_H */
