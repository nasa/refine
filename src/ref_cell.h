
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

#ifndef REF_CELL_H
#define REF_CELL_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CELL_STRUCT REF_CELL_STRUCT;
typedef REF_CELL_STRUCT *REF_CELL;
typedef enum REF_CELL_TYPES { /*  0 */ REF_CELL_EDG,
                              /*  1 */ REF_CELL_ED2,
                              /*  2 */ REF_CELL_ED3,
                              /*  3 */ REF_CELL_TRI,
                              /*  4 */ REF_CELL_TR2,
                              /*  5 */ REF_CELL_TR3,
                              /*  6 */ REF_CELL_QUA,
                              /*  7 */ REF_CELL_TET,
                              /*  8 */ REF_CELL_PYR,
                              /*  9 */ REF_CELL_PRI,
                              /* 10 */ REF_CELL_HEX } REF_CELL_TYPE;
#define REF_CELL_N_TYPE (11)
#define REF_CELL_MAX_SIZE_PER (11)
#define REF_CELL_MAX_NODE_PER (10)
#define REF_CELL_MAX_FACE_PER (6)
END_C_DECLORATION

#include "ref_adj.h"
#include "ref_mpi.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

struct REF_CELL_STRUCT {
  REF_CELL_TYPE type;
  REF_BOOL last_node_is_an_id;
  REF_INT size_per, node_per, edge_per, face_per;
  REF_INT *e2n;
  REF_INT *f2n;
  REF_INT n, max;
  REF_INT blank;
  REF_INT *c2n;
  REF_ADJ ref_adj;
};

#define ref_cell_type(ref_cell) ((ref_cell)->type)
#define ref_cell_last_node_is_an_id(ref_cell) ((ref_cell)->last_node_is_an_id)
#define ref_cell_id_index(ref_cell) ((ref_cell)->node_per)

#define ref_cell_size_per(ref_cell) ((ref_cell)->size_per)
#define ref_cell_node_per(ref_cell) ((ref_cell)->node_per)
#define ref_cell_edge_per(ref_cell) ((ref_cell)->edge_per)
#define ref_cell_face_per(ref_cell) ((ref_cell)->face_per)

#define ref_cell_n(ref_cell) ((ref_cell)->n)
#define ref_cell_max(ref_cell) ((ref_cell)->max)
#define ref_cell_blank(ref_cell) ((ref_cell)->blank)
#define ref_cell_adj(ref_cell) ((ref_cell)->ref_adj)

#define ref_cell_valid(ref_cell, cell)          \
  ((cell) >= 0 && (cell) < ((ref_cell)->max) && \
   REF_EMPTY != (ref_cell)->c2n[ref_cell_size_per(ref_cell) * (cell)])

#define ref_cell_c2n(ref_cell, node, cell) \
  ((ref_cell)->c2n[(node) + ref_cell_size_per(ref_cell) * (cell)])

#define ref_cell_e2n_gen(ref_cell, node, edge) \
  ((ref_cell)->e2n[(node) + 2 * (edge)])

#define ref_cell_e2n(ref_cell, node, cell_edge, cell)          \
  ((ref_cell)->c2n[(ref_cell)->e2n[(node) + 2 * (cell_edge)] + \
                   ref_cell_size_per(ref_cell) * (cell)])

#define ref_cell_f2n_gen(ref_cell, node, face) \
  ((ref_cell)->f2n[(node) + 4 * (face)])

#define ref_cell_f2n(ref_cell, node, cell_face, cell)          \
  ((ref_cell)->c2n[(ref_cell)->f2n[(node) + 4 * (cell_face)] + \
                   ref_cell_size_per(ref_cell) * (cell)])

#define ref_cell_node_empty(ref_cell, node) \
  ref_adj_empty((ref_cell)->ref_adj, node)

#define ref_cell_first_with(ref_cell, node) \
  ref_adj_safe_ref(ref_cell_adj(ref_cell),  \
                   ref_adj_first(ref_cell_adj(ref_cell), (node)))

#define each_ref_cell_valid_cell(ref_cell, cell)              \
  for ((cell) = 0; (cell) < ref_cell_max(ref_cell); (cell)++) \
    if (ref_cell_valid(ref_cell, cell))

#define each_ref_cell_having_node(ref_cell, node, item, cell) \
  each_ref_adj_node_item_with_ref((ref_cell)->ref_adj, node, item, cell)

#define each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node,  \
                                   cell)                                     \
  each_ref_cell_having_node(                                                 \
      ref_cell, node0, item,                                                 \
      cell) for ((cell_node) = 0; (cell_node) < ref_cell_node_per(ref_cell); \
                 (cell_node)++) if ((node1) ==                               \
                                    ref_cell_c2n(ref_cell, cell_node, cell))

#define each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) \
  for ((cell) = 0; (cell) < ref_cell_max(ref_cell); (cell)++)      \
    if (REF_SUCCESS == ref_cell_nodes(ref_cell, cell, nodes))

#define each_ref_cell_cell_node(ref_cell, cell_node)               \
  for ((cell_node) = 0; (cell_node) < ref_cell_node_per(ref_cell); \
       (cell_node)++)

#define each_ref_cell_cell_edge(ref_cell, cell_edge)               \
  for ((cell_edge) = 0; (cell_edge) < ref_cell_edge_per(ref_cell); \
       (cell_edge)++)

#define each_ref_cell_cell_face(ref_cell, cell_face)               \
  for ((cell_face) = 0; (cell_face) < ref_cell_face_per(ref_cell); \
       (cell_face)++)

REF_STATUS ref_cell_create(REF_CELL *ref_cell, REF_CELL_TYPE type);
REF_STATUS ref_cell_free(REF_CELL ref_cell);

REF_STATUS ref_cell_deep_copy(REF_CELL *ref_cell, REF_CELL original);
REF_STATUS ref_cell_pack(REF_CELL ref_cell, REF_INT *o2n);

REF_STATUS ref_cell_meshb_keyword(REF_CELL ref_cell, REF_INT *keyword);

REF_STATUS ref_cell_inspect(REF_CELL ref_cell);
REF_STATUS ref_cell_tattle(REF_CELL ref_cell, REF_INT cell);

REF_STATUS ref_cell_add(REF_CELL ref_cell, REF_INT *nodes, REF_INT *cell);

REF_STATUS ref_cell_add_many_global(REF_CELL ref_cell, REF_NODE ref_node,
                                    REF_INT n, REF_GLOB *c2n, REF_INT *part,
                                    REF_INT exclude_part_id);

REF_STATUS ref_cell_remove(REF_CELL ref_cell, REF_INT cell);
REF_STATUS ref_cell_replace_whole(REF_CELL ref_cell, REF_INT cell,
                                  REF_INT *nodes);
REF_STATUS ref_cell_replace_node(REF_CELL ref_cell, REF_INT old_node,
                                 REF_INT new_node);
REF_STATUS ref_cell_compact(REF_CELL ref_cell, REF_INT **o2n, REF_INT **n2o);

REF_STATUS ref_cell_nodes(REF_CELL ref_cell, REF_INT cell, REF_INT *nodes);
REF_STATUS ref_cell_part_cell_node(REF_CELL ref_cell, REF_NODE ref_node,
                                   REF_INT cell, REF_INT *cell_node);
REF_STATUS ref_cell_part(REF_CELL ref_cell, REF_NODE ref_node, REF_INT cell,
                         REF_INT *part);
REF_STATUS ref_cell_all_local(REF_CELL ref_cell, REF_NODE ref_node,
                              REF_INT cell, REF_BOOL *all_local_nodes);
REF_STATUS ref_cell_local_gem(REF_CELL ref_cell, REF_NODE ref_node,
                              REF_INT node0, REF_INT node1, REF_BOOL *local);
REF_STATUS ref_cell_ncell(REF_CELL ref_cell, REF_NODE ref_node,
                          REF_LONG *ncell);

REF_STATUS ref_cell_orient_node0(REF_INT nnode, REF_INT node0, REF_INT *nodes);

REF_STATUS ref_cell_has_side(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                             REF_BOOL *has_side);

REF_STATUS ref_cell_side_has_id(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                                REF_INT id, REF_BOOL *has_id);
REF_STATUS ref_cell_id_range(REF_CELL ref_cell, REF_MPI ref_mpi,
                             REF_INT *min_id, REF_INT *max_id);

REF_STATUS ref_cell_with_face(REF_CELL ref_cell, REF_INT *face_nodes,
                              REF_INT *cell0, REF_INT *cell1);

REF_STATUS ref_cell_ntri_with_tet_nodes(REF_CELL ref_cell, REF_INT *nodes,
                                        REF_INT *ntri);

REF_STATUS ref_cell_with(REF_CELL ref_cell, REF_INT *nodes, REF_INT *cell);

REF_STATUS ref_cell_degree_with2(REF_CELL ref_cell, REF_INT node0,
                                 REF_INT node1, REF_INT *degree);
REF_STATUS ref_cell_list_with2(REF_CELL ref_cell, REF_INT node0, REF_INT node1,
                               REF_INT max_cell, REF_INT *ncell,
                               REF_INT *cell_list);

REF_STATUS ref_cell_node_list_around(REF_CELL ref_cell, REF_INT node,
                                     REF_INT max_node, REF_INT *nnode,
                                     REF_INT *node_list);

REF_STATUS ref_cell_id_list_around(REF_CELL ref_cell, REF_INT node,
                                   REF_INT max_faceid, REF_INT *nfaceid,
                                   REF_INT *faceids);
REF_STATUS ref_cell_id_list_around_both(REF_CELL ref_cell_a,
                                        REF_CELL ref_cell_b, REF_INT node,
                                        REF_INT max_faceid, REF_INT *nfaceid,
                                        REF_INT *faceids);

REF_STATUS ref_cell_gen_edge_face(REF_CELL ref_cell, REF_INT edge,
                                  REF_INT *face0, REF_INT *face1);

REF_STATUS ref_cell_ghost_long(REF_CELL ref_cell, REF_NODE ref_node,
                               REF_LONG *data);

REF_STATUS ref_cell_global(REF_CELL ref_cell, REF_NODE ref_node,
                           REF_LONG **global);

REF_STATUS ref_cell_tec_fill(REF_CELL ref_cell, const char *filename);

REF_STATUS ref_cell_shape(REF_CELL ref_cell, REF_DBL *bary, REF_DBL *shape);

END_C_DECLORATION

#endif /* REF_CELL_H */

/*
  tecplot "brick"
      7---6
     /|  /|
    4-+-5 |
    | | | |
    | 3-+-2
    |/  |/
    0---1
 */

#define REF_CELL_TEC_BRICK_TET(brick, nodes) \
  {                                          \
    (brick)[0] = (nodes)[0];                 \
    (brick)[1] = (nodes)[1];                 \
    (brick)[2] = (nodes)[2];                 \
    (brick)[3] = (nodes)[2];                 \
    (brick)[4] = (nodes)[3];                 \
    (brick)[5] = (nodes)[3];                 \
    (brick)[6] = (nodes)[3];                 \
    (brick)[7] = (nodes)[3];                 \
  }

#define REF_CELL_TEC_BRICK_PYR(brick, nodes) \
  {                                          \
    (brick)[0] = (nodes)[0];                 \
    (brick)[1] = (nodes)[1];                 \
    (brick)[2] = (nodes)[2];                 \
    (brick)[3] = (nodes)[2];                 \
    (brick)[4] = (nodes)[3];                 \
    (brick)[5] = (nodes)[4];                 \
    (brick)[6] = (nodes)[4];                 \
    (brick)[7] = (nodes)[3];                 \
  }

#define REF_CELL_TEC_BRICK_PRI(brick, nodes) \
  {                                          \
    (brick)[0] = (nodes)[0];                 \
    (brick)[1] = (nodes)[1];                 \
    (brick)[2] = (nodes)[2];                 \
    (brick)[3] = (nodes)[2];                 \
    (brick)[4] = (nodes)[3];                 \
    (brick)[5] = (nodes)[4];                 \
    (brick)[6] = (nodes)[5];                 \
    (brick)[7] = (nodes)[5];                 \
  }

#define REF_CELL_TEC_BRICK_HEX(brick, nodes) \
  {                                          \
    (brick)[0] = (nodes)[0];                 \
    (brick)[1] = (nodes)[1];                 \
    (brick)[2] = (nodes)[2];                 \
    (brick)[3] = (nodes)[3];                 \
    (brick)[4] = (nodes)[4];                 \
    (brick)[5] = (nodes)[5];                 \
    (brick)[6] = (nodes)[6];                 \
    (brick)[7] = (nodes)[7];                 \
  }

/* http://www.simcenter.msstate.edu/docs/solidmesh/ugridconnectivity.html
in c numbering
boundary elements right hand normal point into element, see f2n

                              inode3------5------inode2
                                 / \              . /
                                /   \          .   /
                               /     \      .     /
                              /       \  .       /
                             /        .\        /
                            2      1    4      3
                           /    .        \    /
                          /  .            \  /
                         /.                \/
                      inode0------0------inode1

                     inode3-------7----inode4
                         |    .            | \
                         |       .         |  \
                         |          .      |   \
                         |             5   |    6
                         |                .|     \
                         |                 | .    \
                         |                 |    .  \
                         2                 4       inode2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                         |               . |     /
                         |            .    |    3
                         |         1       |   /
                         |      .          |  /
                         |   .             | /
                         |.                |/
                       inode0------0-----inode1

                                                  inode5
                                                  . /|
                                               .   / |
                                            .     /  |
                                         .       /   |
                                      .         /    |
                                   7           8     5
                                .             /      |
                             .               /       |
                          .                 /        |
                       inode3-----6------inode4    inode2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                         |               . |     /
                         2            .    |    3
                         |         1       4   /
                         |      .          |  /
                         |   .             | /
                         |.                |/
                       inode0-----0------inode1

                               inode7-----11-----inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             9    .           10    6
                            /     7           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4-8----------inode5     |
                         |      inode3.....|...5..inode2
                         |       .         |       /
                         |      .          |      /
                         |     .           |     /
                         2    1            4    3
                         |   .             |   /
                         |  .              |  /
                         | .               | /
                         |.                |/
                       inode0------0-----inode1


                       inode0-----m2----inode1

                       inode0---m2--m3--inode1

                              inode2
                                 / \
                                /   \
                               /     \
                              /       \
                            m5        m4
                            /           \
                           /             \
                          /               \
                         /                 \
                      inode0----m3-----inode1

                              inode2
                                 / \
                                /   \
                              m7    m6
                              /       \
                             /         \
                            /    m9     \
                          m8            m5
                          /               \
                         /                 \
                      inode0--m3----m4---inode1

*/
