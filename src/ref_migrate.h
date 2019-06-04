
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

#ifndef REF_MIGRATE_H
#define REF_MIGRATE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_MIGRATE_STRUCT REF_MIGRATE_STRUCT;
typedef REF_MIGRATE_STRUCT *REF_MIGRATE;
typedef enum REF_MIGRATE_PARTIONERS { /* 0 */ REF_MIGRATE_RECOMMENDED,
                                      /* 1 */ REF_MIGRATE_SINGLE,
                                      /* 2 */ REF_MIGRATE_PARMETIS,
                                      /* 3 */ REF_MIGRATE_ZOLTAN_GRAPH,
                                      /* 4 */ REF_MIGRATE_ZOLTAN_RCB,
                                      /* 5 */ REF_MIGRATE_LAST
} REF_MIGRATE_PARTIONER;
END_C_DECLORATION

#include "ref_adj.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_MIGRATE_STRUCT {
  REF_GRID grid;
  REF_ADJ parent_local;
  REF_ADJ parent_part;
  REF_ADJ conn;
  REF_INT max;
  REF_GLOB *global;
  REF_DBL *xyz;
  REF_DBL *weight;
  REF_INT *age;
};

#define ref_migrate_grid(ref_migrate) ((ref_migrate)->grid)
#define ref_migrate_parent_local(ref_migrate) ((ref_migrate)->parent_local)
#define ref_migrate_parent_part(ref_migrate) ((ref_migrate)->parent_part)
#define ref_migrate_conn(ref_migrate) ((ref_migrate)->conn)

#define ref_migrate_max(ref_migrate) ((ref_migrate)->max)

#define ref_migrate_global(ref_migrate, node) ((ref_migrate)->global[(node)])
#define ref_migrate_valid(ref_migrate, node) \
  (REF_EMPTY != ref_migrate_global(ref_migrate, node))

#define ref_migrate_xyz(ref_migrate, ixyz, node) \
  ((ref_migrate)->xyz[(ixyz) + 3 * (node)])
#define ref_migrate_weight(ref_migrate, node) ((ref_migrate)->weight[(node)])
#define ref_migrate_age(ref_migrate, node) ((ref_migrate)->age[(node)])

#define each_ref_migrate_node(ref_migrate, node)                    \
  for ((node) = 0; (node) < ref_migrate_max(ref_migrate); (node)++) \
    if (ref_migrate_valid(ref_migrate, node))

REF_STATUS ref_migrate_create(REF_MIGRATE *ref_migrate, REF_GRID ref_grid);
REF_STATUS ref_migrate_free(REF_MIGRATE ref_migrate);

REF_STATUS ref_migrate_inspect(REF_MIGRATE ref_migrate);

REF_STATUS ref_migrate_2d_agglomeration_keep(REF_MIGRATE ref_migrate,
                                             REF_INT keep, REF_INT lose);
REF_STATUS ref_migrate_2d_agglomeration(REF_MIGRATE ref_migrate);

REF_STATUS ref_migrate_shufflin_cell(REF_NODE ref_node, REF_CELL ref_cell);
REF_STATUS ref_migrate_shufflin(REF_GRID ref_grid);

REF_STATUS ref_migrate_to_balance(REF_GRID ref_grid);

REF_ULONG ref_migrate_morton_id(REF_UINT x, REF_UINT y, REF_UINT z);

END_C_DECLORATION

#endif /* REF_MIGRATE_H */
