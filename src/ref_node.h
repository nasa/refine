
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

#ifndef REF_NODE_H
#define REF_NODE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_NODE_STRUCT REF_NODE_STRUCT;
typedef REF_NODE_STRUCT *REF_NODE;
END_C_DECLORATION

#include "ref_mpi.h"

BEGIN_C_DECLORATION

struct REF_NODE_STRUCT {
  REF_INT n, max;
  REF_INT blank;
  REF_GLOB *global;
  REF_GLOB *sorted_global;
  REF_GLOB *sorted_local;
  REF_INT *part;
  REF_INT *age;
  REF_DBL *real;
  REF_INT naux;
  REF_DBL *aux;
  REF_MPI ref_mpi;
  REF_INT n_unused, max_unused;
  REF_GLOB *unused_global;
  REF_GLOB old_n_global, new_n_global;
  REF_DBL twod_mid_plane;
  REF_DBL min_volume;
  REF_DBL min_uv_area;
  REF_INT tet_quality;
  REF_INT tri_quality;
};

#define REF_NODE_REAL_PER (15) /* x,y,z, m[6], log_m[6] */

#define REF_NODE_EPIC_QUALITY (1)
#define REF_NODE_JAC_QUALITY (2)

#define ref_node_n(ref_node) ((ref_node)->n)
#define ref_node_max(ref_node) ((ref_node)->max)

#define ref_node_n_global(ref_node) ((ref_node)->old_n_global)

#define ref_node_valid(ref_node, node)               \
  ((node) > -1 && (node) < ref_node_max(ref_node) && \
   (ref_node)->global[(node)] >= 0)

#define ref_node_global(ref_node, node)               \
  (((node) > -1 && (node) < ref_node_max(ref_node) && \
    (ref_node)->global[(node)] >= 0)                  \
       ? (ref_node)->global[(node)]                   \
       : REF_EMPTY)

#define each_ref_node_valid_node(ref_node, node)              \
  for ((node) = 0; (node) < ref_node_max(ref_node); (node)++) \
    if (ref_node_valid(ref_node, node))

#define ref_node_xyz(ref_node, ixyz, node) \
  ((ref_node)->real[(ixyz) + REF_NODE_REAL_PER * (node)])
#define ref_node_xyz_ptr(ref_node, node) \
  (&((ref_node)->real[REF_NODE_REAL_PER * (node)]))

#define ref_node_real(ref_node, ireal, node) \
  ((ref_node)->real[(ireal) + REF_NODE_REAL_PER * (node)])

#define ref_node_owned(ref_node, node) \
  (ref_mpi_rank(ref_node_mpi(ref_node)) == ref_node_part(ref_node, node))

#define ref_node_part(ref_node, node) ((ref_node)->part[(node)])
#define ref_node_age(ref_node, node) ((ref_node)->age[(node)])

#define ref_node_naux(ref_node) ((ref_node)->naux)
#define ref_node_aux(ref_node, iaux, node) \
  ((ref_node)->aux[(iaux) + ref_node_naux(ref_node) * (node)])

#define ref_node_mpi(ref_node) ((ref_node)->ref_mpi)

#define ref_node_n_unused(ref_node) ((ref_node)->n_unused)
#define ref_node_max_unused(ref_node) ((ref_node)->max_unused)

#define ref_node_twod_mid_plane(ref_node) ((ref_node)->twod_mid_plane)
#define ref_node_min_volume(ref_node) ((ref_node)->min_volume)
#define ref_node_min_uv_area(ref_node) ((ref_node)->min_uv_area)

REF_STATUS ref_node_create(REF_NODE *ref_node, REF_MPI ref_mpi);
REF_STATUS ref_node_free(REF_NODE ref_node);

REF_STATUS ref_node_deep_copy(REF_NODE *ref_node_ptr, REF_NODE original);
REF_STATUS ref_node_pack(REF_NODE ref_node, REF_INT *o2n, REF_INT *n2o);

REF_STATUS ref_node_inspect(REF_NODE ref_node);
REF_STATUS ref_node_location(REF_NODE ref_node, REF_INT node);
REF_STATUS ref_node_tattle_global(REF_NODE ref_node, REF_INT global);

REF_STATUS ref_node_local(REF_NODE ref_node, REF_GLOB global, REF_INT *node);

REF_STATUS ref_node_initialize_n_global(REF_NODE ref_node, REF_GLOB n_global);
REF_STATUS ref_node_next_global(REF_NODE ref_node, REF_GLOB *global);

REF_STATUS ref_node_synchronize_globals(REF_NODE ref_node);
REF_STATUS ref_node_shift_new_globals(REF_NODE ref_node);
REF_STATUS ref_node_eliminate_unused_globals(REF_NODE ref_node);
REF_STATUS ref_node_collect_ghost_age(REF_NODE ref_node);

REF_STATUS ref_node_add(REF_NODE ref_node, REF_GLOB global, REF_INT *node);
REF_STATUS ref_node_add_many(REF_NODE ref_node, REF_INT n, REF_GLOB *global);

REF_STATUS ref_node_remove(REF_NODE ref_node, REF_INT node);
REF_STATUS ref_node_remove_invalidates_sorted(REF_NODE ref_node, REF_INT node);
REF_STATUS ref_node_remove_without_global(REF_NODE ref_node, REF_INT node);
REF_STATUS ref_node_remove_without_global_invalidates_sorted(REF_NODE ref_node,
                                                             REF_INT node);
REF_STATUS ref_node_rebuild_sorted_global(REF_NODE ref_node);
REF_STATUS ref_node_implicit_global_from_local(REF_NODE ref_node);

REF_STATUS ref_node_compact(REF_NODE ref_node, REF_INT **o2n, REF_INT **n2o);

REF_STATUS ref_node_ghost_real(REF_NODE ref_node);
REF_STATUS ref_node_ghost_int(REF_NODE ref_node, REF_INT *vector, REF_INT ldim);
REF_STATUS ref_node_ghost_glob(REF_NODE ref_node, REF_GLOB *vector, REF_INT ldim);
REF_STATUS ref_node_ghost_dbl(REF_NODE ref_node, REF_DBL *vector, REF_INT ldim);
REF_STATUS ref_node_localize_ghost_int(REF_NODE ref_node, REF_INT *scalar);

REF_STATUS ref_node_edge_twod(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                              REF_BOOL *twod);

REF_STATUS ref_node_node_twod(REF_NODE ref_node, REF_INT node, REF_BOOL *twod);

REF_STATUS ref_node_metric_form(REF_NODE ref_node, REF_INT node, REF_DBL m11,
                                REF_DBL m12, REF_DBL m13, REF_DBL m22,
                                REF_DBL m23, REF_DBL m33);
REF_STATUS ref_node_metric_set(REF_NODE ref_node, REF_INT node, REF_DBL *m);
REF_STATUS ref_node_metric_get(REF_NODE ref_node, REF_INT node, REF_DBL *m);
REF_STATUS ref_node_metric_set_log(REF_NODE ref_node, REF_INT node,
                                   REF_DBL *log_m);
REF_STATUS ref_node_metric_get_log(REF_NODE ref_node, REF_INT node,
                                   REF_DBL *log_m);

REF_STATUS ref_node_ratio(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                          REF_DBL *ratio);
REF_STATUS ref_node_dratio_dnode0(REF_NODE ref_node, REF_INT node0,
                                  REF_INT node1, REF_DBL *ratio,
                                  REF_DBL *dratio_dnode0);
REF_STATUS ref_node_ratio_node0(REF_NODE ref_node, REF_INT node0, REF_INT node1,
                                REF_DBL *ratio_node0);

REF_STATUS ref_node_tri_normal(REF_NODE ref_node, REF_INT *nodes,
                               REF_DBL *normal);
REF_STATUS ref_node_tri_centroid(REF_NODE ref_node, REF_INT *nodes,
                                 REF_DBL *centroid);

REF_STATUS ref_node_tri_y_projection(REF_NODE ref_node, REF_INT *nodes,
                                     REF_DBL *y_projection);
REF_STATUS ref_node_tri_twod_orientation(REF_NODE ref_node, REF_INT *nodes,
                                         REF_BOOL *valid);
REF_STATUS ref_node_tri_node_angle(REF_NODE ref_node, REF_INT *nodes,
                                   REF_INT node, REF_DBL *angle);
REF_STATUS ref_node_tri_area(REF_NODE ref_node, REF_INT *nodes, REF_DBL *area);
REF_STATUS ref_node_tri_darea_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                     REF_DBL *area, REF_DBL *darea_dnode0);

REF_STATUS ref_node_tri_quality(REF_NODE ref_node, REF_INT *nodes,
                                REF_DBL *quality);
REF_STATUS ref_node_tri_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                        REF_DBL *quality,
                                        REF_DBL *dquality_dnode0);

REF_STATUS ref_node_xyz_vol(REF_DBL *xyzs[4], REF_DBL *volume);
REF_STATUS ref_node_tet_vol(REF_NODE ref_node, REF_INT *nodes, REF_DBL *volume);
REF_STATUS ref_node_tet_dvol_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                    REF_DBL *vol, REF_DBL *dvol_dnode0);

REF_STATUS ref_node_tet_quality(REF_NODE ref_node, REF_INT *nodes,
                                REF_DBL *quality);
REF_STATUS ref_node_tet_dquality_dnode0(REF_NODE ref_node, REF_INT *nodes,
                                        REF_DBL *quality,
                                        REF_DBL *dquality_dnode0);

REF_STATUS ref_node_twod_clone(REF_NODE ref_node, REF_INT original,
                               REF_INT *clone);
REF_STATUS ref_node_interpolate_edge(REF_NODE ref_node, REF_INT node0,
                                     REF_INT node1, REF_DBL node1_weight,
                                     REF_INT new_node);
REF_STATUS ref_node_interpolate_face(REF_NODE ref_node, REF_INT node0,
                                     REF_INT node1, REF_INT node2,
                                     REF_INT new_node);
REF_STATUS ref_node_resize_aux(REF_NODE ref_node);

REF_STATUS ref_node_bary3(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                          REF_DBL *bary);
REF_STATUS ref_node_bary3d(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                           REF_DBL *bary);
REF_STATUS ref_node_bary4(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                          REF_DBL *bary);
REF_STATUS ref_node_clip_bary4(REF_DBL *orig_bary, REF_DBL *clip_bary);

REF_STATUS ref_node_tri_projection(REF_NODE ref_node, REF_INT *nodes,
                                   REF_DBL *xyz, REF_DBL *projection);
REF_STATUS ref_node_dist_to_edge(REF_NODE ref_node, REF_INT *nodes,
                                 REF_DBL *xyz, REF_DBL *distance);
REF_STATUS ref_node_dist_to_tri(REF_NODE ref_node, REF_INT *nodes, REF_DBL *xyz,
                                REF_DBL *distance);

REF_STATUS ref_node_xyz_grad(REF_DBL *xyzs[4], REF_DBL *scalar,
                             REF_DBL *gradient);
REF_STATUS ref_node_tet_grad_nodes(REF_NODE ref_node, REF_INT *nodes,
                                   REF_DBL *scalar, REF_DBL *gradient);

REF_STATUS ref_node_nearest_xyz(REF_NODE ref_node, REF_DBL *xyz,
                                REF_INT *closest_node, REF_DBL *distance);

REF_STATUS ref_node_selection(REF_NODE ref_node, REF_DBL *elements,
                              REF_INT position, REF_DBL *value);

REF_STATUS ref_node_push_unused(REF_NODE ref_node, REF_GLOB unused_global);
REF_STATUS ref_node_pop_unused(REF_NODE ref_node, REF_GLOB *new_global);
REF_STATUS ref_node_shift_unused(REF_NODE ref_node, REF_GLOB equal_and_above,
                                 REF_GLOB shift);
REF_STATUS ref_node_sort_unused(REF_NODE ref_node);
REF_STATUS ref_node_erase_unused(REF_NODE ref_node);
REF_STATUS ref_node_allgather_unused(REF_NODE ref_node);

END_C_DECLORATION

#endif /* REF_NODE_H */
