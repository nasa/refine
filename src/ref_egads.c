
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

#include "ref_cloud.h"
#include "ref_malloc.h"
#include "ref_math.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_EGADS
#include "egads.h"
#endif

REF_STATUS ref_egads_open(REF_GEOM ref_geom) {
#ifdef HAVE_EGADS
  ego context;
  REIS(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  RAS(EG_setOutLevel(context, 0) >= 0, "make silent");
  ref_geom->context = (void *)context;
#else
  ref_geom->context = NULL;
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_egads_close(REF_GEOM ref_geom) {
#ifdef HAVE_EGADS
  if (NULL != ref_geom->faces) EG_free((ego *)(ref_geom->faces));
  if (NULL != ref_geom->edges) EG_free((ego *)(ref_geom->edges));
  if (NULL != ref_geom->nodes) EG_free((ego *)(ref_geom->nodes));
  if (NULL != ref_geom->context)
    REIS(EGADS_SUCCESS, EG_close((ego)(ref_geom->context)), "EG close");
#endif

  ref_geom->context = NULL;
  ref_geom->solid = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_egads_load(REF_GEOM ref_geom, const char *filename) {
#ifdef HAVE_EGADS
  ego context;
  ego model = NULL;
  ego geom, *bodies, *children;
  int oclass, mtype, nbody, *senses, nchild;
  ego solid, *faces, *edges, *nodes;
  int nface, nedge, nnode;
  REF_INT face;
  ego ref;
  double uv_box[4];

  context = (ego)(ref_geom->context);

#ifdef HAVE_EGADS_LITE
  {
    /* entry point NOT in egads.h */
    int EG_importModel(egObject * context, const size_t nbytes,
                       const char stream[], egObject **model);

    SUPRESS_UNUSED_COMPILER_WARNING(filename);

    RAS(0 < ref_geom_cad_data_size(ref_geom), "zero size cad_data");
    RNS(ref_geom_cad_data(ref_geom), "cad_data NULL");
    REIS(EGADS_SUCCESS,
         EG_importModel(context, (size_t)ref_geom_cad_data_size(ref_geom),
                        (char *)ref_geom_cad_data(ref_geom), &model),
         "EG load");
  }
#else
  RNS(filename, "filename NULL for EGADS(full) load");
  REIS(EGADS_SUCCESS, EG_loadModel(context, 0, filename, &model), "EG load");

  {
    /* entry point NOT in egads.h */
    int EG_exportModel(ego mobject, size_t * nbytes, char *stream[]);

    REF_SIZE cad_data_size;
    REF_BYTE *cad_data;

    REIS(EGADS_SUCCESS, EG_exportModel(model, &cad_data_size, &cad_data),
         "EG stream");
    ref_geom_cad_data_size(ref_geom) = cad_data_size;
    /* safe non-NULL free, if already allocated, to prevent memory leaks */
    ref_free(ref_geom->cad_data);
    ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                      ref_geom_cad_data_size(ref_geom), REF_BYTE);
    memcpy(ref_geom_cad_data(ref_geom), cad_data,
           ref_geom_cad_data_size(ref_geom));
    EG_free(cad_data);
  }
#endif

  REIS(EGADS_SUCCESS,
       EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                      &senses),
       "EG topo bodies");
  REIS(1, nbody, "expected 1 body");
  solid = bodies[0];
  REIS(EGADS_SUCCESS,
       EG_getTopology(solid, &geom, &oclass, &mtype, NULL, &nchild, &children,
                      &senses),
       "EG topo body type");
  RAB(SOLIDBODY == mtype || FACEBODY == mtype, "expected SOLIDBODY or FACEBODY",
      { printf("mtype %d\n", mtype); });
  ref_geom->solid = (void *)solid;
  ref_geom->manifold = FACEBODY != mtype;

  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, NODE, &nnode, &nodes),
       "EG node topo");
  ref_geom->nnode = nnode;
  ref_geom->nodes = (void *)nodes;
  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EDGE, &nedge, &edges),
       "EG edge topo");
  ref_geom->nedge = nedge;
  ref_geom->edges = (void *)edges;
  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, FACE, &nface, &faces),
       "EG face topo");
  ref_geom->nface = nface;
  ref_geom->faces = (void *)faces;

  /* use face mtype SFORWARD, SREVERSE to set uv_area_sign */
  /* If it is SFORWARD (1) then the Face's Normal is in the same direction as
     the surface (u cross v), which points outward of the solid. If it is
     SREVERSE (-1), then the natural surface normal points inward and the
     Face points consistently out of the solid. */
  ref_malloc_init(ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 0.0);
  for (face = 0; face < nface; face++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->faces))[face], &ref, &oclass, &mtype,
                        uv_box, &nchild, &children, &senses),
         "topo");
    switch (mtype) {
      /* refine assumes normal point into the domain (solid), flip sign */
      case SFORWARD:
        (ref_geom->uv_area_sign)[face] = -1.0;
        break;
      case SREVERSE:
        (ref_geom->uv_area_sign)[face] = 1.0;
        break;
      default:
        printf("mtype %d\n", mtype);
        RSS(REF_IMPLEMENT, "unknown face type, expected SFORWARD or SREVERSE");
    }
  }

  ref_malloc_init(ref_geom->initial_cell_height, ref_geom->nface, REF_DBL,
                  -1.0);
  for (face = 0; face < nface; face++) {
    int len, atype;
    const double *preals;
    const int *pints;
    const char *string;
    if (EGADS_SUCCESS == EG_attributeRet(((ego *)(ref_geom->faces))[face],
                                         "initial_cell_height", &atype, &len,
                                         &pints, &preals, &string)) {
      if (ATTRREAL == atype && len == 1) {
        ref_geom->initial_cell_height[face] = preals[0];
      }
    }
  }

  ref_malloc_init(ref_geom->face_seg_per_rad, ref_geom->nface, REF_DBL, -999.0);
  for (face = 0; face < nface; face++) {
    int len, atype;
    const double *preals;
    const int *pints;
    const char *string;
    if (EGADS_SUCCESS == EG_attributeRet(((ego *)(ref_geom->faces))[face],
                                         "seg_per_rad", &atype, &len, &pints,
                                         &preals, &string)) {
      if (ATTRREAL == atype && len == 1) {
        ref_geom->face_seg_per_rad[face] = preals[0];
      }
    }
  }

#else
  printf("nothing for %s, No EGADS linked for %s\n", __func__, filename);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif

  return REF_SUCCESS;
}

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_face_surface_type(REF_GEOM ref_geom, REF_INT faceid,
                                              int *surface_type) {
  ego esurf, *eloops, eref;
  int oclass, mtype, nloop, *senses, *pinfo;
  double data[18], *preal;

  RNS(ref_geom->faces, "faces not loaded");
  if (faceid < 1 || faceid > ref_geom->nface) return REF_INVALID;

  REIS(EGADS_SUCCESS,
       EG_getTopology(((ego *)(ref_geom->faces))[faceid - 1], &esurf, &oclass,
                      &mtype, data, &nloop, &eloops, &senses),
       "topo");
  REIS(EGADS_SUCCESS,
       EG_getGeometry(esurf, &oclass, surface_type, &eref, &pinfo, &preal),
       "geom");
  EG_free(pinfo);
  EG_free(preal);
  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_tess_fill_vertex(REF_GRID ref_grid, ego tess,
                                             REF_GLOB *n_global) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);

  int node, new_node, pty, pin;
  double verts[3];
  int tess_status, nvert;
  ego geom;

  REIS(EGADS_SUCCESS, EG_statusTessBody(tess, &geom, &tess_status, &nvert),
       "EG tess");
  REIS(1, tess_status, "tess not closed");

  for (node = 0; node < nvert; node++) {
    REIS(EGADS_SUCCESS, EG_getGlobal(tess, node + 1, &pty, &pin, verts),
         "global node info");
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    REIS(node, new_node, "node index");
    ref_node_xyz(ref_node, 0, node) = verts[0];
    ref_node_xyz(ref_node, 1, node) = verts[1];
    ref_node_xyz(ref_node, 2, node) = verts[2];
    /* pty: point type (-) Face local index, (0) Node, (+) Edge local index */
    if (0 == pty) {
      RSS(ref_geom_add(ref_geom, node, REF_GEOM_NODE, pin, NULL), "node");
    }
  }

  *n_global = (REF_GLOB)nvert;

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_tess_fill_tri(REF_GRID ref_grid, ego tess) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  int face, tlen, plen;
  REF_INT node, tri, new_cell;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;
  REF_DBL param[2];

  for (face = 0; face < (ref_geom->nface); face++) {
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    for (node = 0; node < plen; node++) {
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, node + 1, &(nodes[0])), "l2g0");
      nodes[0] -= 1;
      param[0] = uv[0 + 2 * node];
      param[1] = uv[1 + 2 * node];
      RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_FACE, face + 1, param),
          "face uv");
    }
    for (tri = 0; tri < tlen; tri++) {
      /* triangle orientation flipped, per refine convention */
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[0 + 3 * tri], &(nodes[1])),
           "l2g0");
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[1 + 3 * tri], &(nodes[0])),
           "l2g1");
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[2 + 3 * tri], &(nodes[2])),
           "l2g2");
      nodes[0] -= 1;
      nodes[1] -= 1;
      nodes[2] -= 1;
      nodes[3] = face + 1;
      RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell), "new tri");
    }
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_tess_fill_edg(REF_GRID ref_grid, ego tess) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL param[2];
  int node, edge, plen;
  const double *points, *t;
  REF_INT new_cell;

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    int egads_status;
    REF_BOOL degenerate;
    degenerate = REF_FALSE;
    REIS(EGADS_SUCCESS, EG_getTessEdge(tess, edge + 1, &plen, &points, &t),
         "tess query edge");
    for (node = 0; node < plen; node++) {
      egads_status = EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0]));
      if (EGADS_DEGEN == egads_status) {
        degenerate = REF_TRUE;
      } else {
        REIS(EGADS_SUCCESS, egads_status, "l2g0");
        nodes[0] -= 1;
        param[0] = t[node];
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, edge + 1, param),
            "edge t");
      }
    }
    if (!degenerate)
      for (node = 0; node < (plen - 1); node++) {
        /* assue edge index is 1-bias */
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0])),
             "l2g0");
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 2, &(nodes[1])),
             "l2g1");
        nodes[0] -= 1;
        nodes[1] -= 1;
        nodes[2] = edge + 1;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell), "new edge");
      }
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_merge_tparams(REF_CLOUD object_tp_augment,
                                          REF_INT id, REF_DBL *new_params) {
  REF_DBL params[3];
  REF_INT i, item;

  if (ref_cloud_has_global(object_tp_augment, (REF_GLOB)id)) {
    RSS(ref_cloud_item(object_tp_augment, (REF_GLOB)id, &item),
        "find existing entry");
    each_ref_cloud_aux(object_tp_augment, i) {
      params[i] = ref_cloud_aux(object_tp_augment, i, item);
      params[i] = MIN(params[i], new_params[i]);
    }
    RSS(ref_cloud_store(object_tp_augment, (REF_GLOB)id, params),
        "cache merged .tParams");
  } else {
    RSS(ref_cloud_store(object_tp_augment, (REF_GLOB)id, new_params),
        "cache new .tParams");
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_face_width(REF_GEOM ref_geom, REF_INT faceid,
                                       REF_CLOUD edge_tp_augment) {
  ego faceobj;
  double diag, box[6];

  double len, xyz0[18], xyz1[18], dx[3], param[2], s;
  REF_INT it, tsamples = 3;
  REF_INT ineligible_cad_node0, ineligible_cad_node1, cad_node0, cad_node1;
  ego edgeobj0, edgeobj1;
  ego *cadnodes0, *cadnodes1;
  ego *edges0, *edges1;
  ego *loops, ref;
  int nloop;
  int ncadnode0, ncadnode1;
  int nedge0, nedge1;
  int *senses, mtype, oclass;
  double trange0[2], trange1[2], data[18];
  REF_INT edge0, edge1;
  REF_INT loop0, loop1;
  REF_DBL width, aspect_ratio, adjusted;
  double params[3];
  int edgeid;

  RAS(0 < faceid && faceid <= (ref_geom->nface), "invalid faceid");
  faceobj = ((ego *)(ref_geom->faces))[faceid - 1];
  REIS(EGADS_SUCCESS, EG_getBoundingBox(faceobj, box), "EG bounding box");
  diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
              (box[1] - box[4]) * (box[1] - box[4]) +
              (box[2] - box[5]) * (box[2] - box[5]));

  REIS(EGADS_SUCCESS,
       EG_getTopology(faceobj, &ref, &oclass, &mtype, data, &nloop, &loops,
                      &senses),
       "topo");
  for (loop0 = 0; loop0 < nloop; loop0++) {
    /* loop through all Edges associated with this Loop */
    REIS(EGADS_SUCCESS,
         EG_getTopology(loops[loop0], &ref, &oclass, &mtype, data, &nedge0,
                        &edges0, &senses),
         "topo");
    for (edge0 = 0; edge0 < nedge0; edge0++) {
      edgeobj0 = edges0[edge0];
      REIS(EGADS_SUCCESS,
           EG_getTopology(edgeobj0, &ref, &oclass, &mtype, trange0, &ncadnode0,
                          &cadnodes0, &senses),
           "EG topo edge0");
      if (mtype == DEGENERATE) continue; /* skip DEGENERATE */
      RAS(0 < ncadnode0 && ncadnode0 < 3, "edge children");
      ineligible_cad_node0 = EG_indexBodyTopo(ref_geom->solid, cadnodes0[0]);
      if (2 == ncadnode0) {
        ineligible_cad_node1 = EG_indexBodyTopo(ref_geom->solid, cadnodes0[1]);
      } else {
        ineligible_cad_node1 = ineligible_cad_node0; /* ONENODE edge */
      }
      width = diag;
      for (loop1 = 0; loop1 < nloop; loop1++) {
        /* loop through all Edges associated with this Loop */
        REIS(EGADS_SUCCESS,
             EG_getTopology(loops[loop1], &ref, &oclass, &mtype, data, &nedge1,
                            &edges1, &senses),
             "topo");
        for (edge1 = 0; edge1 < nedge1; edge1++) {
          edgeobj1 = edges1[edge1];
          REIS(EGADS_SUCCESS,
               EG_getTopology(edgeobj1, &ref, &oclass, &mtype, trange1,
                              &ncadnode1, &cadnodes1, &senses),
               "EG topo edge0");
          if (mtype == DEGENERATE) continue; /* skip DEGENERATE */
          RAS(0 < ncadnode1 && ncadnode1 < 3, "edge children");
          cad_node0 = EG_indexBodyTopo(ref_geom->solid, cadnodes1[0]);
          if (2 == ncadnode1) {
            cad_node1 = EG_indexBodyTopo(ref_geom->solid, cadnodes1[1]);
          } else {
            cad_node1 = cad_node0; /* ONENODE edge */
          }
          if (cad_node0 == ineligible_cad_node0 ||
              cad_node0 == ineligible_cad_node1 ||
              cad_node1 == ineligible_cad_node0 ||
              cad_node1 == ineligible_cad_node1)
            continue;
          for (it = 0; it < tsamples; it++) {
            s = (REF_DBL)(it + 1) / (REF_DBL)(tsamples + 1);
            param[0] = s * trange0[1] + (1.0 - s) * trange0[0];
            REIS(EGADS_SUCCESS, EG_evaluate(edgeobj0, param, xyz0),
                 "sample edge");
            /* inverse projections */
            REIS(EGADS_SUCCESS, EG_invEvaluate(edgeobj1, xyz0, param, xyz1),
                 "inv eval other edge");
            dx[0] = xyz1[0] - xyz0[0];
            dx[1] = xyz1[1] - xyz0[1];
            dx[2] = xyz1[2] - xyz0[2];
            len = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
            width = MIN(width, len);
          }
        }
      }
      edgeid = EG_indexBodyTopo(ref_geom->solid, edgeobj0);
      if (ref_math_divisible(diag, width)) {
        aspect_ratio = diag / width;
        adjusted = MIN(MAX(1.0, aspect_ratio - 10.0), 10.0) * width;
        params[0] = adjusted;
        params[1] = 0.1 * params[0];
        params[2] = 15.0;
        RSS(ref_egads_merge_tparams(edge_tp_augment, edgeid, params),
            "update tparams");
      }
    }
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_adjust_tparams(REF_GEOM ref_geom, ego tess,
                                           REF_CLOUD face_tp_augment,
                                           REF_CLOUD edge_tp_augment,
                                           REF_INT auto_tparams) {
  ego faceobj;
  int face, tlen, plen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;

  double params[3], diag, box[6];

  for (face = 0; face < (ref_geom->nface); face++) {
    faceobj = ((ego *)(ref_geom->faces))[face];
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    if (0 == plen || 0 == tlen) {
      printf("face %d has %d nodes and %d triangles\n", face + 1, plen, tlen);
      if (auto_tparams & REF_EGADS_MISSING_TPARAM) {
        REIS(EGADS_SUCCESS, EG_getBoundingBox(faceobj, box), "EG bounding box");
        diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                    (box[1] - box[4]) * (box[1] - box[4]) +
                    (box[2] - box[5]) * (box[2] - box[5]));
        params[0] = 0.1 * diag;
        params[1] = 0.01 * diag;
        params[2] = 15.0;
        RSS(ref_egads_merge_tparams(face_tp_augment, face + 1, params),
            "update tparams");
      }
    } else {
      REF_INT tri, side, n0, n1, i;
      const REF_DBL *xyz0, *xyz1, *uv0, *uv1;
      REF_DBL uvm[2], xyz[3], dside[3], dmid[3];
      REF_DBL length, offset, max_chord, max_chord_length, max_chord_offset;
      max_chord = 0.0;
      max_chord_length = 0.0;
      max_chord_offset = 0.0;
      for (tri = 0; tri < tlen; tri++) {
        for (side = 0; side < 3; side++) {
          n0 = side;
          n1 = side + 1;
          if (n1 > 2) n1 -= 3;
          n0 = tris[n0 + 3 * tri] - 1;
          n1 = tris[n1 + 3 * tri] - 1;
          xyz0 = &(points[3 * n0]);
          xyz1 = &(points[3 * n1]);
          uv0 = &(uv[2 * n0]);
          uv1 = &(uv[2 * n1]);
          for (i = 0; i < 2; i++) uvm[i] = 0.5 * (uv0[i] + uv1[i]);
          RSS(ref_geom_eval_at(ref_geom, REF_GEOM_FACE, face + 1, uvm, xyz,
                               NULL),
              "eval mid uv");
          for (i = 0; i < 3; i++) dside[i] = xyz1[i] - xyz0[i];
          for (i = 0; i < 3; i++) dmid[i] = xyz[i] - 0.5 * (xyz1[i] + xyz0[i]);
          length = sqrt(dside[0] * dside[0] + dside[1] * dside[1] +
                        dside[2] * dside[2]);
          offset =
              sqrt(dmid[0] * dmid[0] + dmid[1] * dmid[1] + dmid[2] * dmid[2]);
          if (ref_math_divisible(offset, length)) {
            if (max_chord < offset / length) {
              max_chord = offset / length;
              max_chord_length = length;
              max_chord_offset = offset;
            }
          }
        }
      }
      REIS(EGADS_SUCCESS, EG_getBoundingBox(faceobj, box), "EG bounding box");
      diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                  (box[1] - box[4]) * (box[1] - box[4]) +
                  (box[2] - box[5]) * (box[2] - box[5]));
      if (max_chord > 0.2)
        printf("face %d rel chord %f abs len %f abs chord %f diag %f\n",
               face + 1, max_chord, max_chord_length, max_chord_offset, diag);
      if (max_chord > 0.2 && (auto_tparams & REF_EGADS_CHORD_TPARAM)) {
        params[0] = 0.1 * diag;
        params[1] = 0.001 * diag;
        params[2] = 15.0;
        RSS(ref_egads_merge_tparams(face_tp_augment, face + 1, params),
            "update tparams");
      }
    }

    /* face width parameter to all edges */
    if (auto_tparams & REF_EGADS_WIDTH_TPARAM) {
      RSS(ref_egads_face_width(ref_geom, face + 1, edge_tp_augment),
          "face width");
    }
  }
  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_update_tparams_attributes(
    REF_GEOM ref_geom, REF_CLOUD face_tp_original, REF_CLOUD edge_tp_original,
    REF_CLOUD face_tp_augment, REF_CLOUD edge_tp_augment, REF_BOOL *rebuild) {
  ego faceobj;
  int face;
  ego edgeobj;
  int edge;
  int i, item;
  int len, atype;
  const double *preals;
  const int *pints;
  const char *string;
  double params[3];
  REF_BOOL needs_update;
  REF_DBL tol = 1.0e-7;

  *rebuild = REF_FALSE;
  for (face = 0; face < (ref_geom->nface); face++) {
    faceobj = ((ego *)(ref_geom->faces))[face];
    /* respect manually set .tParams */
    if (ref_cloud_has_global(face_tp_original, (REF_GLOB)(face + 1))) continue;
    /* merge update with existing attributes */
    if (ref_cloud_has_global(face_tp_augment, (REF_GLOB)(face + 1))) {
      needs_update = REF_FALSE;
      RSS(ref_cloud_item(face_tp_augment, (REF_GLOB)(face + 1), &item),
          "find existing entry");
      each_ref_cloud_aux(face_tp_augment, i) {
        params[i] = ref_cloud_aux(face_tp_augment, i, item);
      }
      if (EGADS_SUCCESS == EG_attributeRet(faceobj, ".tParams", &atype, &len,
                                           &pints, &preals, &string)) {
        /* examine exisiting .tParams and detect change within tol */
        if (ATTRREAL != atype || len != 3) THROW("malformed .tParams");
        each_ref_cloud_aux(face_tp_augment, i) {
          needs_update =
              needs_update || (ABS(params[i] - preals[i]) > tol * params[i]);
        }
      } else {
        /* create new .tParams attribute */
        needs_update = REF_TRUE;
      }
      if (needs_update) {
#ifdef HAVE_EGADS_LITE
        RSS(REF_IMPLEMENT, "full EGADS required to adjust .tParams");
#else
        REIS(EGADS_SUCCESS,
             EG_attributeAdd(faceobj, ".tParams", ATTRREAL, 3, NULL, params,
                             NULL),
             "set new .tParams");
        printf("select face %d\nattribute .tParams  %f;%f;%f\n", face + 1,
               params[0], params[1], params[2]);
#endif
        *rebuild = REF_TRUE;
      }
    }
  }
  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    edgeobj = ((ego *)(ref_geom->edges))[edge];
    /* respect manually set .tParams */
    if (ref_cloud_has_global(edge_tp_original, (REF_GLOB)(edge + 1))) continue;
    /* merge update with existing attributes */
    if (ref_cloud_has_global(edge_tp_augment, (REF_GLOB)(edge + 1))) {
      needs_update = REF_FALSE;
      RSS(ref_cloud_item(edge_tp_augment, (REF_GLOB)(edge + 1), &item),
          "find existing entry");
      each_ref_cloud_aux(edge_tp_augment, i) {
        params[i] = ref_cloud_aux(edge_tp_augment, i, item);
      }
      if (EGADS_SUCCESS == EG_attributeRet(edgeobj, ".tParams", &atype, &len,
                                           &pints, &preals, &string)) {
        /* examine exisiting .tParams and detect change within tol */
        if (ATTRREAL != atype || len != 3) THROW("malformed .tParams");
        each_ref_cloud_aux(edge_tp_augment, i) {
          needs_update =
              needs_update || (ABS(params[i] - preals[i]) > tol * params[i]);
        }
      } else {
        /* create new .tParams attribute */
        needs_update = REF_TRUE;
      }
      if (needs_update) {
#ifdef HAVE_EGADS_LITE
        RSS(REF_IMPLEMENT, "full EGADS required to adjust .tParams");
#else
        REIS(EGADS_SUCCESS,
             EG_attributeAdd(edgeobj, ".tParams", ATTRREAL, 3, NULL, params,
                             NULL),
             "set new .tParams");
        printf("select edge %d\nattribute .tParams  %f;%f;%f\n", edge + 1,
               params[0], params[1], params[2]);
#endif
        *rebuild = REF_TRUE;
      }
    }
  }
  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_cache_tparams(REF_GEOM ref_geom,
                                          REF_CLOUD face_tp_original,
                                          REF_CLOUD edge_tp_original) {
  ego faceobj, edgeobj;
  int face, edge;
  int i;
  int len, atype;
  const double *preals;
  const int *pints;
  const char *string;
  double params[3];

  for (face = 0; face < (ref_geom->nface); face++) {
    faceobj = ((ego *)(ref_geom->faces))[face];
    if (EGADS_SUCCESS == EG_attributeRet(faceobj, ".tParams", &atype, &len,
                                         &pints, &preals, &string)) {
      if (ATTRREAL == atype && len == 3) {
        for (i = 0; i < 3; i++) params[i] = preals[i];
        RSS(ref_cloud_store(face_tp_original, (REF_GLOB)(face + 1), params),
            "cache original .tParams");
      } else {
        printf("  wrong format .tParams atype %d len %d\n", atype, len);
      }
    }
  }

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    edgeobj = ((ego *)(ref_geom->edges))[edge];
    if (EGADS_SUCCESS == EG_attributeRet(edgeobj, ".tParams", &atype, &len,
                                         &pints, &preals, &string)) {
      if (ATTRREAL == atype && len == 3) {
        for (i = 0; i < 3; i++) params[i] = preals[i];
        RSS(ref_cloud_store(edge_tp_original, (REF_GLOB)(edge + 1), params),
            "cache original .tParams");
      } else {
        printf("  wrong format .tParams atype %d len %d\n", atype, len);
      }
    }
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_tess_create(REF_GEOM ref_geom, ego *tess,
                                        REF_INT auto_tparams) {
  ego solid, geom;
  int tess_status, nvert;
  double params[3], diag, box[6];
  REF_BOOL rebuild;
  REF_INT tries;
  REF_CLOUD face_tp_original, face_tp_augment;
  REF_CLOUD edge_tp_original, edge_tp_augment;
  solid = (ego)(ref_geom->solid);
  /* maximum length of an EDGE segment or triangle side (in physical space) */
  /* curvature-based value that looks locally at the deviation between
     the centroid of the discrete object and the underlying geometry */
  /* maximum interior dihedral angle (in degrees) */

  RSS(ref_cloud_create(&face_tp_original, 3), "create tparams cache");
  RSS(ref_cloud_create(&edge_tp_original, 3), "create tparams cache");
  RSS(ref_cloud_create(&face_tp_augment, 3), "create tparams augment");
  RSS(ref_cloud_create(&edge_tp_augment, 3), "create tparams augment");
  RSS(ref_egads_cache_tparams(ref_geom, face_tp_original, edge_tp_original),
      "tparams cache");

  REIS(EGADS_SUCCESS, EG_getBoundingBox(solid, box), "EG bounding box");
  diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
              (box[1] - box[4]) * (box[1] - box[4]) +
              (box[2] - box[5]) * (box[2] - box[5]));

  params[0] = 0.25 * diag;
  params[1] = 0.001 * diag;
  params[2] = 15.0;

  rebuild = REF_TRUE;
  tries = 0;
  while (rebuild && tries < 5) {
    tries++;
    REIS(EGADS_SUCCESS, EG_makeTessBody(solid, params, tess), "EG tess");
    REIS(EGADS_SUCCESS, EG_statusTessBody(*tess, &geom, &tess_status, &nvert),
         "EG tess");
    REIS(1, tess_status, "tess not closed");

    RSS(ref_egads_adjust_tparams(ref_geom, *tess, face_tp_augment,
                                 edge_tp_augment, auto_tparams),
        "adjust params");
    RSS(ref_egads_update_tparams_attributes(ref_geom, face_tp_original,
                                            edge_tp_original, face_tp_augment,
                                            edge_tp_augment, &rebuild),
        "adjust params");

    if (rebuild)
      printf("rebuild EGADS tessilation after .tParams adjustment, try %d\n",
             tries);
  }

  RSS(ref_cloud_free(face_tp_augment), "free tparams augment");
  RSS(ref_cloud_free(edge_tp_original), "free tparams cache");
  RSS(ref_cloud_free(face_tp_original), "free tparams cache");

  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_tess(REF_GRID ref_grid, REF_INT auto_tparams) {
#ifdef HAVE_EGADS
  ego tess;
  REF_GLOB n_global;

  if (ref_mpi_once(ref_grid_mpi(ref_grid))) {
    RSS(ref_egads_tess_create(ref_grid_geom(ref_grid), &tess, auto_tparams),
        "create tess object");

    RSS(ref_egads_tess_fill_vertex(ref_grid, tess, &n_global),
        "fill tess vertex");
    RSS(ref_egads_tess_fill_tri(ref_grid, tess), "fill tess triangles");
    RSS(ref_egads_tess_fill_edg(ref_grid, tess), "fill tess edges");
  }

  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &n_global, 1, REF_GLOB_TYPE),
      "bcast glob");
  RSS(ref_node_initialize_n_global(ref_grid_node(ref_grid), n_global),
      "init glob");

  RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps");
  ref_grid_surf(ref_grid) = REF_TRUE;

#else
  printf("returning empty grid from %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  SUPRESS_UNUSED_COMPILER_WARNING(auto_tparams);
#endif

  return REF_SUCCESS;
}

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_edge_faces(REF_GRID ref_grid,
                                       REF_INT **edge_face_arg) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT *e2f, *nface;
  REF_INT face, edge;

  ego esurf, *eloops;
  int oclass, mtype, nloop, *senses;
  double data[18];
  ego ecurve, *eedges;
  int iloop, iedge, nedge;

  ref_malloc_init(*edge_face_arg, 2 * (ref_geom->nedge), REF_INT, REF_EMPTY);
  e2f = *edge_face_arg;
  ref_malloc_init(nface, (ref_geom->nedge), REF_INT, 0);

  for (face = 0; face < (ref_geom->nface); face++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->faces))[face], &esurf, &oclass,
                        &mtype, data, &nloop, &eloops, &senses),
         "topo");
    for (iloop = 0; iloop < nloop; iloop++) {
      /* loop through all Edges associated with this Loop */
      REIS(EGADS_SUCCESS,
           EG_getTopology(eloops[iloop], &ecurve, &oclass, &mtype, data, &nedge,
                          &eedges, &senses),
           "topo");
      for (iedge = 0; iedge < nedge; iedge++) {
        edge = EG_indexBodyTopo((ego)(ref_geom->solid), eedges[iedge]) - 1;
        RAS(2 > nface[edge], "edge has more than 2 faces");
        e2f[nface[edge] + 2 * edge] = face + 1;
        nface[edge]++;
      }
    }
  }

  ref_free(nface);
  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_node_faces(REF_GRID ref_grid,
                                       REF_ADJ *ref_adj_arg) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj;
  REF_INT *e2f, id, toponode;
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double trange[2];
  RSS(ref_adj_create(ref_adj_arg), "create ref_adj");
  ref_adj = *ref_adj_arg;
  RSS(ref_egads_edge_faces(ref_grid, &e2f), "edge2face");
  for (id = 1; id <= ref_geom->nedge; id++) {
    object = ((ego *)(ref_geom->edges))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");

    if (0 < nchild) {
      toponode = EG_indexBodyTopo(ref_geom->solid, pchldrn[0]);
      if (0 < e2f[0 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[0 + 2 * (id - 1)]),
            "add");
      }
      if (0 < e2f[1 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[1 + 2 * (id - 1)]),
            "add");
      }
    }

    if (1 < nchild) {
      toponode = EG_indexBodyTopo(ref_geom->solid, pchldrn[1]);
      if (0 < e2f[0 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[0 + 2 * (id - 1)]),
            "add");
      }
      if (0 < e2f[1 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[1 + 2 * (id - 1)]),
            "add");
      }
    }
  }
  ref_free(e2f);
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_mark_jump_degen(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT node, geom, edge, face, cad_node;
  REF_INT nfound, node_geom, edge_geom, face_geom;
  ego eref;
  int oclass, mtype, *senses;
  double trange[2];
  double uv[2];
  ego *echilds;
  int nchild;
  int *e2f;
  REF_DBL du, dv;
  REF_DBL xyz[3], dxyz_duv[15];
  REF_INT geom_node_id, degen;

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->edges))[edge], &eref, &oclass,
                        &mtype, trange, &nchild, &echilds, &senses),
         "edge topo");
    if (mtype == ONENODE) {
      REIS(1, nchild, "ONENODE should have one node");
      cad_node = EG_indexBodyTopo(ref_geom->solid, echilds[0]);
      if (ref_grid_once(ref_grid)) {
        printf("edge id %d is ONENODE at geom node %d\n", edge + 1, cad_node);
      }
      node = REF_EMPTY;
      each_ref_geom_node(ref_geom, geom) {
        if (cad_node == ref_geom_id(ref_geom, geom)) {
          node = ref_geom_node(ref_geom, geom);
        }
      }
      nfound = 0;
      each_ref_geom_edge(ref_geom, geom) {
        if (node == ref_geom_node(ref_geom, geom) &&
            edge + 1 == ref_geom_id(ref_geom, geom)) {
          ref_geom_jump(ref_geom, geom) = edge + 1;
          REIS(0, nfound, "edge geom already found");
          nfound++;
        }
      }
    }
  }

  RSS(ref_egads_edge_faces(ref_grid, &e2f), "edge2face");

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->edges))[edge], &eref, &oclass,
                        &mtype, trange, &nchild, &echilds, &senses),
         "edge topo");
    if (mtype == DEGENERATE) {
      geom_node_id = EG_indexBodyTopo((ego)(ref_geom->solid), echilds[0]);
      face = e2f[0 + 2 * edge] - 1;
      REIS(REF_EMPTY, e2f[1 + 2 * edge], "DEGENERATE edge has two faces");

      face_geom = REF_EMPTY;
      each_ref_geom_node(ref_geom, node_geom) {
        if (geom_node_id == ref_geom_id(ref_geom, node_geom)) {
          node = ref_geom_node(ref_geom, node_geom);
          RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, face + 1,
                            &face_geom),
              "face for degen edge at node not found");
        }
      }

      /* in parallel, may not have this node */
      if (REF_EMPTY != face_geom) {
        uv[0] = ref_geom_param(ref_geom, 0, face_geom);
        uv[1] = ref_geom_param(ref_geom, 1, face_geom);
        RSS(ref_geom_eval_at(ref_geom, REF_GEOM_FACE, face + 1, uv, xyz,
                             dxyz_duv),
            "eval at");
        du = sqrt(ref_math_dot(&(dxyz_duv[0]), &(dxyz_duv[0])));
        dv = sqrt(ref_math_dot(&(dxyz_duv[3]), &(dxyz_duv[3])));
        if (du > dv) {
          degen = (edge + 1); /* positive edge, trust uv[0], larger dxyz/du */
        } else {
          degen = -(edge + 1); /* negative edge, trust uv[1], larger dxyz/dv */
        }
        ref_geom_degen(ref_geom, face_geom) = degen;
      }

      if (ref_grid_once(ref_grid)) {
        printf("edge id %d is degen for face id %d\n", edge + 1, face + 1);
      }
    }
  }

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    if (e2f[0 + 2 * edge] == e2f[1 + 2 * edge]) {
      nfound = 0;
      each_ref_geom_edge(ref_geom, edge_geom) {
        if (edge + 1 == ref_geom_id(ref_geom, edge_geom)) {
          each_ref_geom_face(ref_geom, face_geom) {
            if (e2f[0 + 2 * edge] == ref_geom_id(ref_geom, face_geom) &&
                ref_geom_node(ref_geom, edge_geom) ==
                    ref_geom_node(ref_geom, face_geom)) {
              ref_geom_jump(ref_geom, face_geom) = edge + 1;
              nfound++;
            }
          }
        }
      }
      if (ref_grid_once(ref_grid)) {
        printf("edge id %d is used twice by face id %d\n", edge + 1,
               e2f[0 + 2 * edge]);
      }
    }
  }

  ref_free(e2f);
#else
  printf("unable to %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_recon_nodes(REF_GRID ref_grid,
                                        REF_INT **cad_nodes) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_ADJ n2f;
  REF_INT id, item, best_node, node, faceid;
  REF_DBL best_dist, dist;
  REF_INT max_faceids = 50;
  REF_INT *grid_faceids, grid_nfaceids;
  REF_INT *cad_faceids, cad_nfaceids;
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double xyz[3];
  REF_BOOL show_xyz = REF_FALSE;
  REF_INT i, j;
  REF_BOOL found, all_found;
  RSS(ref_egads_node_faces(ref_grid, &n2f), "build n2f");
  ref_malloc(grid_faceids, max_faceids, REF_INT);
  ref_malloc(cad_faceids, max_faceids, REF_INT);
  printf("searching for %d topo nodes\n", ref_geom->nnode);
  ref_malloc(*cad_nodes, ref_geom->nnode, REF_INT);
  for (id = 1; id <= ref_geom->nnode; id++) {
    object = ((ego *)(ref_geom->nodes))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, xyz, &nchild, &pchldrn,
                        &psens),
         "EG topo node");
    cad_nfaceids = 0;
    each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
      cad_faceids[cad_nfaceids] = faceid;
      cad_nfaceids++;
    }
    best_node = REF_EMPTY;
    best_dist = 1.0e20;
    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node, max_faceids,
                                  &grid_nfaceids, grid_faceids),
          "count faceids");
      if (grid_nfaceids != cad_nfaceids) continue;
      all_found = REF_TRUE;
      for (i = 0; i < cad_nfaceids; i++) {
        found = REF_FALSE;
        for (j = 0; j < grid_nfaceids; j++) {
          found = found || (cad_faceids[i] == grid_faceids[j]);
        }
        all_found = all_found && found;
      }
      if (!all_found) continue;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      if (dist < best_dist) {
        best_node = node;
        best_dist = dist;
      }
    }
    if (REF_EMPTY == best_node) {
      printf("cad node %d expects faces", id);
      each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
        printf(" %d", faceid);
      }
      printf("\n");
      printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      THROW("unable to find a grid vertex that matches cad topo");
    }
    printf(" topo node id %3d node %6d dist %.4e fid", id, best_node,
           best_dist);
    RSS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node, max_faceids,
                                &grid_nfaceids, grid_faceids),
        "count faceids");
    for (i = 0; i < grid_nfaceids; i++) printf(" %d", grid_faceids[i]);
    printf(" expects");
    each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
      printf(" %d", faceid);
    }
    printf("\n");
    if (show_xyz) {
      printf(" d %23.15e %23.15e %23.15e\n",
             ref_node_xyz(ref_node, 0, best_node),
             ref_node_xyz(ref_node, 1, best_node),
             ref_node_xyz(ref_node, 2, best_node));
      printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
    }
    (*cad_nodes)[id - 1] = best_node;
    RSS(ref_geom_add(ref_geom, best_node, REF_GEOM_NODE, id, NULL), "node");
  }
  ref_free(cad_faceids);
  ref_free(grid_faceids);
  ref_adj_free(n2f);
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_recon(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double xyz[9], trange[2];
  REF_INT node, id, best_node;
  REF_DBL best_dist, dist, best_param;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *cad_nodes;
  double t;
  double param[2];
  REF_INT i, cell, edge_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pass, updates;
  REF_BOOL show_xyz = REF_FALSE;
  REF_BOOL debug = REF_FALSE;
#define REF_GEOM_MAX_FACEIDS (50)
  REF_INT nfaceid, faceids[REF_GEOM_MAX_FACEIDS];
  REF_INT nfaceid0, faceids0[REF_GEOM_MAX_FACEIDS];
  REF_INT nfaceid1, faceids1[REF_GEOM_MAX_FACEIDS];
  REF_INT *e2f;
  REF_INT degree, max_node = 50, *node_list;

  /* to allow recon after meshb load times */
  RSS(ref_geom_initialize(ref_geom), "clear out previous assoc");
  RSS(ref_cell_free(ref_grid_edg(ref_grid)), "clear out edge");
  RSS(ref_cell_create(&ref_grid_edg(ref_grid), REF_CELL_EDG), "edg");

  RSS(ref_egads_recon_nodes(ref_grid, &cad_nodes), "recover nodes");

  RSS(ref_egads_edge_faces(ref_grid, &e2f), "compute edge faces");

  ref_malloc(node_list, max_node, REF_INT);

  for (id = 1; id <= ref_geom->nedge; id++) {
    object = ((ego *)(ref_geom->edges))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");
    if (mtype == DEGENERATE) {
      printf(" topo edge id %3d degen\n", id);
    } else {
      int toponode0, toponode1;
      REF_INT node0, node1;
      double closest[9];
      REF_INT next_node, current_node;
      REF_INT geom;
      REIS(TWONODE, mtype, "ONENODE edge not implemented");
      REIS(2, nchild, "expect two topo node for edge");
      toponode0 = EG_indexBodyTopo(ref_geom->solid, pchldrn[0]);
      toponode1 = EG_indexBodyTopo(ref_geom->solid, pchldrn[1]);
      node0 = cad_nodes[toponode0 - 1];
      node1 = cad_nodes[toponode1 - 1];
      printf(" topo edge id %3d fid %d %d trange [%f,%f]\n", id,
             e2f[0 + 2 * (id - 1)], e2f[1 + 2 * (id - 1)], trange[0],
             trange[1]);
      REIS(EGADS_SUCCESS, EG_evaluate(object, &(trange[0]), xyz), "EG eval");
      node = node0;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      printf("  node0 id %2d index %3d t %f dist %e fid", toponode0, node0,
             trange[0], dist);
      RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node0,
                                  REF_GEOM_MAX_FACEIDS, &nfaceid0, faceids0),
          REF_INCREASE_LIMIT, "count faceids");
      for (i = 0; i < nfaceid0; i++) printf(" %d", faceids0[i]);
      printf("\n");
      if (show_xyz) {
        printf(" d %23.15e %23.15e %23.15e\n", ref_node_xyz(ref_node, 0, node),
               ref_node_xyz(ref_node, 1, node),
               ref_node_xyz(ref_node, 2, node));
        printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      }
      REIS(EGADS_SUCCESS, EG_evaluate(object, &(trange[1]), xyz), "EG eval");
      node = node1;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      printf("  node1 id %2d index %3d t %f dist %e fid", toponode1, node1,
             trange[1], dist);
      RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node1,
                                  REF_GEOM_MAX_FACEIDS, &nfaceid1, faceids1),
          REF_INCREASE_LIMIT, "count faceids");
      for (i = 0; i < nfaceid1; i++) printf(" %d", faceids1[i]);
      printf("\n");
      if (show_xyz) {
        printf(" d %23.15e %23.15e %23.15e\n", ref_node_xyz(ref_node, 0, node),
               ref_node_xyz(ref_node, 1, node),
               ref_node_xyz(ref_node, 2, node));
        printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      }
      current_node = node0;
      t = trange[0];
      {
        best_node = node0;
        param[0] = trange[0];
        RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
                                  ref_node_xyz_ptr(ref_node, best_node), param),
            "inv wrapper");
        REIS(EGADS_SUCCESS, EG_evaluate(object, &(param[0]), closest),
             "EG eval");
        dist = sqrt(pow(closest[0] - ref_node_xyz(ref_node, 0, best_node), 2) +
                    pow(closest[1] - ref_node_xyz(ref_node, 1, best_node), 2) +
                    pow(closest[2] - ref_node_xyz(ref_node, 2, best_node), 2));
        printf("   best_node %5d t %f best_dist %e fid", best_node, param[0],
               dist);
        RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node,
                                    REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
            REF_INCREASE_LIMIT, "count faceids");
        for (i = 0; i < nfaceid; i++) printf(" %d", faceids[i]);
        printf("\n");
        if (param[0] < trange[0] || trange[1] < param[0])
          THROW("edge t out of range");
        if (show_xyz) {
          printf(" d %23.15e %23.15e %23.15e\n",
                 ref_node_xyz(ref_node, 0, best_node),
                 ref_node_xyz(ref_node, 1, best_node),
                 ref_node_xyz(ref_node, 2, best_node));
          printf(" c %23.15e %23.15e %23.15e\n", closest[0], closest[1],
                 closest[2]);
        }
      }
      while (current_node != node1) {
        RSS(ref_cell_node_list_around(ref_grid_tri(ref_grid), current_node,
                                      max_node, &degree, node_list),
            "next node");
        best_node = REF_EMPTY;
        best_dist = 1.0e20;
        best_param = 1.0e20;
        for (i = 0; i < degree; i++) {
          REF_INT fid;
          REF_BOOL have_faceid0, have_faceid1;
          next_node = node_list[i];
          if (REF_SUCCESS ==
              ref_geom_find(ref_geom, next_node, REF_GEOM_EDGE, id, &geom))
            continue; /* this canidate is already part of the edge */
          RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), next_node,
                                      REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
              REF_INCREASE_LIMIT, "count faceids");
          if (nfaceid < 2) continue; /* should be between faces */
          have_faceid0 = REF_FALSE;
          for (fid = 0; fid < nfaceid; fid++)
            if (e2f[0 + 2 * (id - 1)] == faceids[fid]) have_faceid0 = REF_TRUE;
          have_faceid1 = REF_FALSE;
          for (fid = 0; fid < nfaceid; fid++)
            if (e2f[1 + 2 * (id - 1)] == faceids[fid]) have_faceid1 = REF_TRUE;
          if (!have_faceid0 || !have_faceid1)
            continue; /* must have expected faceids */
          param[0] = t;
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
                                    ref_node_xyz_ptr(ref_node, next_node),
                                    param),
              "inv wrapper");
          REIS(EGADS_SUCCESS, EG_evaluate(object, &(param[0]), closest),
               "EG eval");
          dist =
              sqrt(pow(closest[0] - ref_node_xyz(ref_node, 0, next_node), 2) +
                   pow(closest[1] - ref_node_xyz(ref_node, 1, next_node), 2) +
                   pow(closest[2] - ref_node_xyz(ref_node, 2, next_node), 2));
          if (dist < best_dist) {
            best_node = next_node;
            best_dist = dist;
            best_param = param[0];
          }
        }
        if (REF_EMPTY == best_node) {
          RSS(REF_FAILURE, "count not find next node");
        }
        param[0] = best_param;
        printf("   best_node %5d t %f best_dist %e fid", best_node, param[0],
               best_dist);
        RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node,
                                    REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
            REF_INCREASE_LIMIT, "count faceids");
        for (i = 0; i < nfaceid; i++) printf(" %d", faceids[i]);
        printf("\n");
        if (param[0] < trange[0] || trange[1] < param[0])
          THROW("edge t out of range");
        if (show_xyz) {
          printf(" d %23.15e %23.15e %23.15e\n",
                 ref_node_xyz(ref_node, 0, best_node),
                 ref_node_xyz(ref_node, 1, best_node),
                 ref_node_xyz(ref_node, 2, best_node));
          printf(" c %23.15e %23.15e %23.15e\n", closest[0], closest[1],
                 closest[2]);
        }
        param[1] = t;
        nodes[0] = best_node;
        nodes[1] = current_node;
        nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &cell), "add e");
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, id, &(param[0])),
            "add geom edge 0");
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_EDGE, id, &(param[1])),
            "add geom edge 1");
        current_node = best_node;
        t = param[0];
      }
    }
    if (debug) ref_geom_tec(ref_grid, "ref_geom_each_edge.tec");
  }
  ref_free(cad_nodes);
  ref_free(node_list);
  ref_free(e2f);

  printf("seting face UV under edges\n");
  each_ref_cell_valid_cell_with_nodes(ref_grid_edg(ref_grid), cell,
                                      edge_nodes) {
    REF_INT edgeid, faceid;
    int sense;
    REF_INT ncell, tri_list[2], tri_nodes[REF_CELL_MAX_SIZE_PER];
    edgeid = edge_nodes[2];
    RSS(ref_cell_list_with2(ref_grid_tri(ref_grid), edge_nodes[0],
                            edge_nodes[1], 2, &ncell, tri_list),
        "tri list for edge");
    REIS(2, ncell, "expect two tri for edge");
    for (i = 0; i < ncell; i++) {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri_list[i], tri_nodes),
          "nodes");
      faceid = tri_nodes[3];
      sense = 0; /* when edge is used twice in loop, this is 1 or -1 */
      RSS(ref_geom_tuv(ref_geom, edge_nodes[0], REF_GEOM_EDGE, edgeid, &t),
          "edge t0");
      REIS(
          EGADS_SUCCESS,
          EG_getEdgeUV(((ego *)(ref_geom->faces))[faceid - 1],
                       ((ego *)(ref_geom->edges))[edgeid - 1], sense, t, param),
          "eval edge face uv");
      RSS(ref_geom_add(ref_geom, edge_nodes[0], REF_GEOM_FACE, faceid, param),
          "add geom face for edge");
      RSS(ref_geom_tuv(ref_geom, edge_nodes[1], REF_GEOM_EDGE, edgeid, &t),
          "edge t0");
      REIS(
          EGADS_SUCCESS,
          EG_getEdgeUV(((ego *)(ref_geom->faces))[faceid - 1],
                       ((ego *)(ref_geom->edges))[edgeid - 1], sense, t, param),
          "eval edge face uv");
      RSS(ref_geom_add(ref_geom, edge_nodes[1], REF_GEOM_FACE, faceid, param),
          "add geom face for edge");
    }
  }

  pass = 0;
  updates = REF_EMPTY;
  while (0 != updates) {
    REF_INT geom0, geom1, geom2;
    REF_INT faceid;
    int surface_type;
    double closest[3];
    REF_BOOL inv_eval_wrapper = REF_TRUE;
    updates = 0;
    pass++;
    each_ref_cell_valid_cell_with_nodes(ref_grid_tri(ref_grid), cell, nodes) {
      faceid = nodes[3];
      RSS(ref_egads_face_surface_type(ref_geom, faceid, &surface_type), "styp");
      RXS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_FACE, faceid, &geom0),
          REF_NOT_FOUND, "find0");
      RXS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_FACE, faceid, &geom1),
          REF_NOT_FOUND, "find1");
      RXS(ref_geom_find(ref_geom, nodes[2], REF_GEOM_FACE, faceid, &geom2),
          REF_NOT_FOUND, "find2");
      if (REF_EMPTY != geom0 && REF_EMPTY == geom1) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[0], REF_GEOM_FACE, faceid, param),
            "geom0");
        if (inv_eval_wrapper && surface_type == PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[1]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[1]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_FACE, faceid, param),
            "add face");
      }
      if (REF_EMPTY != geom1 && REF_EMPTY == geom2) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[1], REF_GEOM_FACE, faceid, param),
            "geom1");
        if (inv_eval_wrapper && surface_type == PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[2]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[2]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[2], REF_GEOM_FACE, faceid, param),
            "add face");
      }
      if (REF_EMPTY != geom2 && REF_EMPTY == geom0) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[2], REF_GEOM_FACE, faceid, param),
            "geom2");
        if (inv_eval_wrapper && surface_type == PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[0]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[0]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_FACE, faceid, param),
            "add face");
      }
    }
    printf(" pass %3d updates %d\n", pass, updates);
  }
  return REF_SUCCESS;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  printf("No EGADS linked for %s\n", __func__);
  return REF_IMPLEMENT;
#endif
}
