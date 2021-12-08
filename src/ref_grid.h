
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

#ifndef REF_GRID_H
#define REF_GRID_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_GRID_STRUCT REF_GRID_STRUCT;
typedef REF_GRID_STRUCT *REF_GRID;
END_C_DECLORATION

#include "ref_adapt.h"
#include "ref_cell.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_interp.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

typedef enum REF_GRID_COORDSYSS { /* 0 */ REF_GRID_XBYRZU,
                                  /* 1 */ REF_GRID_XBYUZL,
                                  /* 2 */ REF_GRID_XFYRZD,
                                  /* 3 */ REF_GRID_COORDSYS_LAST
} REF_GRID_COORDSYS;
typedef enum REF_GRID_UNITS { /* 0 */ REF_GRID_IN,
                              /* 1 */ REF_GRID_FT,
                              /* 2 */ REF_GRID_M,
                              /* 3 */ REF_GRID_CM,
                              /* 4 */ REF_GRID_UNIT_LAST } REF_GRID_UNIT;

struct REF_GRID_STRUCT {
  REF_MPI mpi;
  REF_NODE node;

  REF_CELL cell[REF_CELL_N_TYPE + 1];

  REF_GEOM geom;
  REF_GATHER gather;
  REF_ADAPT adapt;

  REF_INTERP interp;

  REF_MIGRATE_PARTIONER partitioner;
  REF_INT partitioner_seed;
  REF_BOOL partitioner_full;

  REF_INT meshb_version;
  REF_GRID_COORDSYS coordinate_system;
  REF_GRID_UNIT unit;
  REF_DBL references[7];

  REF_BOOL twod;
  REF_BOOL surf;
};

REF_STATUS ref_grid_create(REF_GRID *ref_grid, REF_MPI ref_mpi);
REF_STATUS ref_grid_free(REF_GRID ref_grid);

REF_STATUS ref_grid_deep_copy(REF_GRID *ref_grid, REF_GRID original);
REF_STATUS ref_grid_cache_background(REF_GRID ref_grid);
REF_STATUS ref_grid_stable_pack(REF_GRID ref_grid);
REF_STATUS ref_grid_pack(REF_GRID ref_grid);

#define ref_grid_mpi(ref_grid) ((ref_grid)->mpi)
#define ref_grid_once(ref_grid) ref_mpi_once(ref_grid_mpi(ref_grid))

#define ref_grid_node(ref_grid) ((ref_grid)->node)
#define ref_grid_cell(ref_grid, group) ((ref_grid)->cell[(group)])

#define ref_grid_edg(ref_grid) ref_grid_cell(ref_grid, REF_CELL_EDG)
#define ref_grid_ed2(ref_grid) ref_grid_cell(ref_grid, REF_CELL_ED2)
#define ref_grid_ed3(ref_grid) ref_grid_cell(ref_grid, REF_CELL_ED3)
#define ref_grid_tri(ref_grid) ref_grid_cell(ref_grid, REF_CELL_TRI)
#define ref_grid_tr2(ref_grid) ref_grid_cell(ref_grid, REF_CELL_TR2)
#define ref_grid_tr3(ref_grid) ref_grid_cell(ref_grid, REF_CELL_TR3)
#define ref_grid_qua(ref_grid) ref_grid_cell(ref_grid, REF_CELL_QUA)
#define ref_grid_tet(ref_grid) ref_grid_cell(ref_grid, REF_CELL_TET)
#define ref_grid_pyr(ref_grid) ref_grid_cell(ref_grid, REF_CELL_PYR)
#define ref_grid_pri(ref_grid) ref_grid_cell(ref_grid, REF_CELL_PRI)
#define ref_grid_hex(ref_grid) ref_grid_cell(ref_grid, REF_CELL_HEX)

#define ref_grid_geom(ref_grid) ((ref_grid)->geom)
#define ref_grid_gather(ref_grid) ((ref_grid)->gather)
#define ref_grid_adapt(ref_grid, param) (((ref_grid)->adapt)->param)
#define ref_grid_interp(ref_grid) ((ref_grid)->interp)
#define ref_grid_background(ref_grid)  \
  ((NULL == ref_grid_interp(ref_grid)) \
       ? NULL                          \
       : ref_interp_from_grid(ref_grid_interp(ref_grid)))

#define ref_grid_partitioner(ref_grid) ((ref_grid)->partitioner)
#define ref_grid_partitioner_seed(ref_grid) ((ref_grid)->partitioner_seed)
#define ref_grid_partitioner_full(ref_grid) ((ref_grid)->partitioner_full)

#define ref_grid_meshb_version(ref_grid) ((ref_grid)->meshb_version)
#define ref_grid_coordinate_system(ref_grid) ((ref_grid)->coordinate_system)
#define ref_grid_unit(ref_grid) ((ref_grid)->unit)
/* 0-2: reference length in each direction.
 *   3: reference area
 * 4-6: moment center (all used for forces/moments calculations) */
#define ref_grid_reference(ref_grid, ref_index) \
  ((ref_grid)->references[(ref_index)])

#define ref_grid_twod(ref_grid) ((ref_grid)->twod)
#define ref_grid_surf(ref_grid) ((ref_grid)->surf)

#define each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell)     \
  for ((group) = 7, (ref_cell) = ref_grid_cell(ref_grid, group); \
       (group) <= 11; (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_2d_ref_cell(ref_grid, group, ref_cell)                   \
  for ((group) = 3, (ref_cell) = ref_grid_cell(ref_grid, group); (group) <= 6; \
       (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_1d_ref_cell(ref_grid, group, ref_cell)                   \
  for ((group) = 0, (ref_cell) = ref_grid_cell(ref_grid, group); (group) <= 2; \
       (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_2d_3d_ref_cell(ref_grid, group, ref_cell)  \
  for ((group) = 3, (ref_cell) = ref_grid_cell(ref_grid, group); \
       (group) <= 11; (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_all_ref_cell(ref_grid, group, ref_cell)    \
  for ((group) = 0, (ref_cell) = ref_grid_cell(ref_grid, group); \
       (group) < REF_CELL_N_TYPE;                                \
       (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_face_ref_cell(ref_grid, group, ref_cell)                 \
  for ((group) = 3, (ref_cell) = ref_grid_cell(ref_grid, group); (group) <= 5; \
       (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

#define each_ref_grid_edge_ref_cell(ref_grid, group, ref_cell)                 \
  for ((group) = 0, (ref_cell) = ref_grid_cell(ref_grid, group); (group) <= 2; \
       (group)++, (ref_cell) = ref_grid_cell(ref_grid, group))

REF_STATUS ref_grid_parse_coordinate_system(REF_GRID ref_grid,
                                            const char *coordinate_system);
REF_STATUS ref_grid_parse_unit(REF_GRID ref_grid, const char *unit);

REF_STATUS ref_grid_inspect(REF_GRID ref_grid);
REF_STATUS ref_grid_tattle(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_grid_cell_with(REF_GRID ref_grid, REF_INT node_per,
                              REF_CELL *ref_cell);
REF_STATUS ref_grid_face_with(REF_GRID ref_grid, REF_INT node_per,
                              REF_CELL *ref_cell);

REF_STATUS ref_grid_cell_has_face(REF_GRID ref_grid, REF_INT *face_nodes,
                                  REF_BOOL *has_face);

REF_STATUS ref_grid_faceid_range(REF_GRID ref_grid, REF_INT *min_faceid,
                                 REF_INT *max_faceid);

REF_STATUS ref_grid_tri_qua_id_nodes(REF_GRID ref_grid, REF_INT cell_id,
                                     REF_INT *nnode, REF_INT *ncell,
                                     REF_INT **g2l, REF_INT **l2g);
REF_STATUS ref_grid_cell_id_nodes(REF_GRID ref_grid, REF_CELL ref_cell,
                                  REF_INT cell_tag, REF_INT *nnode,
                                  REF_INT *ncell, REF_INT **g2l, REF_INT **l2g);

REF_STATUS ref_grid_compact_cell_nodes(REF_GRID ref_grid, REF_CELL ref_cell,
                                       REF_GLOB *nnode, REF_LONG *ncell,
                                       REF_GLOB **l2c);
REF_STATUS ref_grid_compact_cell_id_nodes(REF_GRID ref_grid, REF_CELL ref_cell,
                                          REF_INT cell_id, REF_GLOB *nnode,
                                          REF_LONG *ncell, REF_GLOB **l2c);
REF_STATUS ref_grid_compact_surf_id_nodes(REF_GRID ref_grid, REF_INT cell_id,
                                          REF_GLOB *nnode, REF_LONG *ncell,
                                          REF_GLOB **l2c);

REF_STATUS ref_grid_inward_boundary_orientation(REF_GRID ref_grid);

REF_STATUS ref_grid_node_list_around(REF_GRID ref_grid, REF_INT node,
                                     REF_INT max_node, REF_INT *nnode,
                                     REF_INT *node_list);

REF_STATUS ref_grid_enclosing_tet(REF_GRID ref_grid, REF_DBL *xyz, REF_INT *tet,
                                  REF_DBL *bary);

REF_STATUS ref_grid_extrude_twod(REF_GRID *extruded, REF_GRID twod,
                                 REF_INT planes);
REF_STATUS ref_grid_orient_edg(REF_GRID ref_grid, REF_INT *nodes);

REF_STATUS ref_grid_drop_volume(REF_GRID ref_grid);

REF_STATUS ref_grid_ncell(REF_GRID ref_grid, REF_INT *ncell);
REF_STATUS ref_grid_contiguous_group_cell(REF_GRID ref_grid,
                                          REF_INT contiguous_cell,
                                          REF_INT *cell_group, REF_INT *cell);
REF_STATUS ref_grid_contiguous_cell_global(REF_GRID ref_grid, REF_LONG *global);

END_C_DECLORATION

#endif /* REF_GRID_H */
