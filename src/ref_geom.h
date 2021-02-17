
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

#ifndef REF_GEOM_H
#define REF_GEOM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_GEOM_STRUCT REF_GEOM_STRUCT;
typedef REF_GEOM_STRUCT *REF_GEOM;

#define REF_GEOM_NODE (0)
#define REF_GEOM_EDGE (1)
#define REF_GEOM_FACE (2)
#define REF_GEOM_SOLID (3)

#define REF_GEOM_DESCR_SIZE (6)
#define REF_GEOM_DESCR_TYPE (0)
#define REF_GEOM_DESCR_ID (1)
#define REF_GEOM_DESCR_GREF (2)
#define REF_GEOM_DESCR_JUMP (3)
#define REF_GEOM_DESCR_DEGEN (4)
#define REF_GEOM_DESCR_NODE (5)

END_C_DECLORATION

#include "ref_adj.h"
#include "ref_facelift.h"
#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_GEOM_STRUCT {
  REF_INT n, max;
  REF_INT blank;
  REF_INT *descr;
  REF_DBL *param;
  REF_DBL *uv_area_sign;
  REF_DBL *initial_cell_height;
  REF_DBL *face_min_length;
  REF_DBL *face_seg_per_rad;
  REF_DBL segments_per_radian_of_curvature;
  REF_DBL segments_per_bounding_box_diagonal;
  REF_DBL tolerance_protection;
  REF_DBL gap_protection;
  REF_ADJ ref_adj;
  REF_INT nnode, nedge, nface;
  REF_BOOL effective;
  REF_BOOL effective_curvature;
  REF_BOOL manifold;
  REF_BOOL contex_owned;
  void *context;
  void *model;
  void *body;
  void *faces;
  void *edges;
  void *nodes;
  REF_SIZE cad_data_size;
  REF_BYTE *cad_data;
  void *meshlink;
  void *meshlink_projection;
  REF_FACELIFT ref_facelift;
};

#define ref_geom_n(ref_geom) ((ref_geom)->n)
#define ref_geom_max(ref_geom) ((ref_geom)->max)
#define ref_geom_blank(ref_geom) ((ref_geom)->blank)
#define ref_geom_adj(ref_geom) ((ref_geom)->ref_adj)
#define ref_geom_effective(ref_geom) ((ref_geom)->effective)
#define ref_geom_effective_curvature(ref_geom) ((ref_geom)->effective_curvature)
#define ref_geom_manifold(ref_geom) ((ref_geom)->manifold)
#define ref_geom_cad_data(ref_geom) ((ref_geom)->cad_data)
#define ref_geom_cad_data_size(ref_geom) ((ref_geom)->cad_data_size)
#define ref_geom_facelift(ref_geom) ((ref_geom)->ref_facelift)

#define ref_geom_model_loaded(ref_geom) (NULL != (void *)((ref_geom)->model))
#define ref_geom_meshlinked(ref_geom) (NULL != (void *)((ref_geom)->meshlink))

#define ref_geom_descr(ref_geom, attribute, geom) \
  ((ref_geom)->descr[(attribute) + REF_GEOM_DESCR_SIZE * (geom)])

#define ref_geom_type(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_TYPE, (geom)))
#define ref_geom_id(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_ID, (geom)))
#define ref_geom_gref(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_GREF, (geom)))
#define ref_geom_jump(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_JUMP, (geom)))
#define ref_geom_degen(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_DEGEN, (geom)))
#define ref_geom_node(ref_geom, geom) \
  (ref_geom_descr((ref_geom), REF_GEOM_DESCR_NODE, (geom)))

#define ref_geom_param(ref_geom, dimension, geom) \
  ((ref_geom)->param[(dimension) + 2 * (geom)])

#define ref_geom_face_initial_cell_height(ref_geom, face) \
  (((face) >= 0 && (face) < (ref_geom)->nface &&          \
    NULL != (ref_geom)->initial_cell_height)              \
       ? (ref_geom)->initial_cell_height[(face)]          \
       : -1.0)
#define ref_geom_face_min_length(ref_geom, face) \
  (((face) >= 0 && (face) < (ref_geom)->nface && \
    NULL != (ref_geom)->face_min_length)         \
       ? (ref_geom)->face_min_length[(face)]     \
       : -1.0)
#define ref_geom_face_segments_per_radian_of_curvature(ref_geom, face) \
  (((face) >= 0 && (face) < (ref_geom)->nface &&                       \
    NULL != (ref_geom)->face_seg_per_rad)                              \
       ? (ref_geom)->face_seg_per_rad[(face)]                          \
       : -999.0)
#define ref_geom_segments_per_radian_of_curvature(ref_geom) \
  ((ref_geom)->segments_per_radian_of_curvature)
#define ref_geom_segments_per_bounding_box_diagonal(ref_geom) \
  ((ref_geom)->segments_per_bounding_box_diagonal)
#define ref_geom_curvature_unlimited(ref_geom) \
  (0.1 > ref_geom_segments_per_radian_of_curvature(ref_geom))
#define ref_geom_tolerance_protection(ref_geom) \
  ((ref_geom)->tolerance_protection)
#define ref_geom_gap_protection(ref_geom) ((ref_geom)->gap_protection)

#define each_ref_type(refx_geom, type) for ((type) = 0; (type) < 3; (type)++)
#define each_ref_descr(ref_geom, item) \
  for ((item) = 0; (item) < REF_GEOM_DESCR_SIZE; (item)++)

#define each_ref_geom(ref_geom, geom)                         \
  for ((geom) = 0; (geom) < ref_geom_max(ref_geom); (geom)++) \
    if (REF_EMPTY != ref_geom_type(ref_geom, geom))

#define each_ref_geom_of(ref_geom, type, geom)                \
  for ((geom) = 0; (geom) < ref_geom_max(ref_geom); (geom)++) \
    if ((type) == ref_geom_type(ref_geom, geom))

#define each_ref_geom_node(ref_geom, geom)                    \
  for ((geom) = 0; (geom) < ref_geom_max(ref_geom); (geom)++) \
    if (REF_GEOM_NODE == ref_geom_type(ref_geom, geom))

#define each_ref_geom_edge(ref_geom, geom)                    \
  for ((geom) = 0; (geom) < ref_geom_max(ref_geom); (geom)++) \
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom))

#define each_ref_geom_face(ref_geom, geom)                    \
  for ((geom) = 0; (geom) < ref_geom_max(ref_geom); (geom)++) \
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom))

#define each_ref_geom_having_node(ref_geom, node, item, geom) \
  each_ref_adj_node_item_with_ref((ref_geom)->ref_adj, node, item, geom)

#define each_ref_geom_node_id(ref_geom, id) \
  for ((id) = 1; (id) <= (ref_geom)->nnode; (id)++)

#define each_ref_geom_edge_id(ref_geom, id) \
  for ((id) = 1; (id) <= (ref_geom)->nedge; (id)++)

#define each_ref_geom_face_id(ref_geom, id) \
  for ((id) = 1; (id) <= (ref_geom)->nface; (id)++)

REF_STATUS ref_geom_create(REF_GEOM *ref_geom);
REF_STATUS ref_geom_initialize(REF_GEOM ref_geom);
REF_STATUS ref_geom_free(REF_GEOM ref_geom);

REF_STATUS ref_geom_deep_copy(REF_GEOM *ref_geom, REF_GEOM original);
REF_STATUS ref_geom_share_context(REF_GEOM ref_geom_recipient,
                                  REF_GEOM ref_geom_donor);
REF_STATUS ref_geom_pack(REF_GEOM ref_geom, REF_INT *o2n);

REF_STATUS ref_geom_infer_nedge_nface(REF_GRID ref_grid);

REF_STATUS ref_geom_uv_area(REF_GEOM ref_geom, REF_INT *nodes,
                            REF_DBL *uv_area);
REF_STATUS ref_geom_uv_area_sign(REF_GRID ref_grid, REF_INT id, REF_DBL *sign);
REF_STATUS ref_geom_uv_area_report(REF_GRID ref_grid);

REF_STATUS ref_geom_inspect(REF_GEOM ref_geom);
REF_STATUS ref_geom_tattle(REF_GEOM ref_geom, REF_INT node);

REF_STATUS ref_geom_supported(REF_GEOM ref_geom, REF_INT node,
                              REF_BOOL *has_support);
REF_STATUS ref_geom_tri_supported(REF_GEOM ref_geom, REF_INT *nodes,
                                  REF_BOOL *has_support);
REF_STATUS ref_geom_id_supported(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                                 REF_INT id, REF_BOOL *has_support);

REF_STATUS ref_geom_add(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param);
REF_STATUS ref_geom_add_with_descr(REF_GEOM ref_geom, REF_INT *descr,
                                   REF_DBL *param);

REF_STATUS ref_geom_remove(REF_GEOM ref_geom, REF_INT geom);
REF_STATUS ref_geom_remove_all(REF_GEOM ref_geom, REF_INT node);
REF_STATUS ref_geom_remove_without_cell(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_geom_is_a(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_BOOL *it_is);
REF_STATUS ref_geom_unique_id(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                              REF_BOOL *id);

REF_STATUS ref_geom_find(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_INT id, REF_INT *geom);

REF_STATUS ref_geom_tuv(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param);
REF_STATUS ref_geom_cell_tuv_supported(REF_GEOM ref_geom, REF_INT *nodes,
                                       REF_INT type, REF_BOOL *supported);
REF_STATUS ref_geom_cell_tuv(REF_GEOM ref_geom, REF_INT node, REF_INT *nodes,
                             REF_INT type, REF_DBL *param, REF_INT *sens);

REF_STATUS ref_geom_between_face_area(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_INT new_node,
                                      const char *msg);
REF_STATUS ref_geom_add_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_DBL node1_weight, REF_INT new_node);
REF_STATUS ref_geom_support_between(REF_GRID ref_grid, REF_INT node0,
                                    REF_INT node1, REF_BOOL *needs_support);

REF_STATUS ref_geom_tri_uv_bounding_box(REF_GRID ref_grid, REF_INT node,
                                        REF_DBL *uv_min, REF_DBL *uv_max);
REF_STATUS ref_geom_tri_uv_bounding_box2(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_DBL *uv_min,
                                         REF_DBL *uv_max);

REF_STATUS ref_geom_constrain_all(REF_GRID ref_grid);
REF_STATUS ref_geom_constrain(REF_GRID ref_grid, REF_INT node);

REF_STATUS ref_geom_radian_request(REF_GEOM ref_geom, REF_INT geom,
                                   REF_DBL *delta_radian);

REF_STATUS ref_geom_face_rsn(REF_GEOM ref_geom, REF_INT faceid, REF_DBL *uv,
                             REF_DBL *r, REF_DBL *s, REF_DBL *n);
REF_STATUS ref_geom_uv_rsn(REF_DBL *uv, REF_DBL *r, REF_DBL *s, REF_DBL *n,
                           REF_DBL *drsduv);
REF_STATUS ref_geom_tri_centroid(REF_GRID ref_grid, REF_INT *nodes,
                                 REF_DBL *uv);
REF_STATUS ref_geom_tri_norm_deviation(REF_GRID ref_grid, REF_INT *nodes,
                                       REF_DBL *dot_product);
REF_STATUS ref_geom_crease(REF_GRID ref_grid, REF_INT node,
                           REF_DBL *dot_product);

REF_STATUS ref_geom_max_gap(REF_GRID ref_grid, REF_DBL *max_gap);
REF_STATUS ref_geom_verify_param(REF_GRID ref_grid);
REF_STATUS ref_geom_verify_topo(REF_GRID ref_grid);

REF_STATUS ref_geom_usable(REF_GEOM ref_geom, REF_INT geom, REF_BOOL *usable);

REF_STATUS ref_geom_reliability(REF_GEOM ref_geom, REF_INT geom, REF_DBL *slop);

REF_STATUS ref_geom_feedback(REF_GRID ref_grid);

REF_STATUS ref_geom_has_jump(REF_GEOM ref_grid, REF_INT node,
                             REF_BOOL *has_jump);

REF_STATUS ref_geom_aflr_volume(REF_GRID ref_grid, const char *project,
                                const char *options);
REF_STATUS ref_geom_tetgen_volume(REF_GRID ref_grid, const char *project,
                                  const char *options);

REF_STATUS ref_geom_edge_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file);
REF_STATUS ref_geom_face_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file);
REF_STATUS ref_geom_norm_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file);
REF_STATUS ref_geom_tec(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_geom_curve_tec(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_geom_tec_para_shard(REF_GRID ref_grid,
                                   const char *root_filename);

REF_STATUS ref_geom_ghost(REF_GEOM ref_geom, REF_NODE ref_node);

REF_STATUS ref_geom_faceid_range(REF_GRID ref_grid, REF_INT *min_faceid,
                                 REF_INT *max_faceid);
REF_STATUS ref_geom_edgeid_range(REF_GRID ref_grid, REF_INT *min_edgeid,
                                 REF_INT *max_edgeid);

REF_STATUS ref_geom_report_tri_area_normdev(REF_GRID ref_grid);

REF_STATUS ref_geom_bary3(REF_GEOM ref_geom, REF_INT *nodes, REF_DBL *uv,
                          REF_DBL *bary);
REF_STATUS ref_geom_tri_uv_bounding_sphere3(REF_GEOM ref_geom, REF_INT *nodes,
                                            REF_DBL *center, REF_DBL *radius);
REF_STATUS ref_geom_bary2(REF_GEOM ref_geom, REF_INT *nodes, REF_DBL t,
                          REF_DBL *bary);
REF_STATUS ref_geom_edg_t_bounding_sphere2(REF_GEOM ref_geom, REF_INT *nodes,
                                           REF_DBL *center, REF_DBL *radius);

REF_STATUS ref_geom_enrich2(REF_GRID ref_grid);
REF_STATUS ref_geom_enrich3(REF_GRID ref_grid);

END_C_DECLORATION

#endif /* REF_GEOM_H */
