
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
  REIS(SOLIDBODY, mtype, "expected SOLIDBODY");
  ref_geom->solid = (void *)solid;

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
static REF_STATUS ref_egads_tess_fill_vertex(REF_GRID ref_grid, ego tess) {
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

  RSS(ref_node_initialize_n_global(ref_node, nvert), "init glob");

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
static REF_STATUS ref_egads_tess_adjust(REF_GEOM ref_geom, ego tess,
                                        REF_BOOL *rebuild) {
  int face, tlen, plen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;

  *rebuild = REF_FALSE;
  for (face = 0; face < (ref_geom->nface); face++) {
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    if (0 == plen || 0 == tlen) {
      printf("face %d has %d nodes and %d triangles\n", face + 1, plen, tlen);
    }
  }
  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_tess_create(REF_GEOM ref_geom, ego *tess) {
  ego solid, geom;
  int tess_status, nvert;
  double params[3], diag, box[6];
  REF_BOOL rebuild;
  solid = (ego)(ref_geom->solid);

  /* maximum length of an EDGE segment or triangle side (in physical space) */
  /* curvature-based value that looks locally at the deviation between
     the centroid of the discrete object and the underlying geometry */
  /* maximum interior dihedral angle (in degrees) */

  REIS(EGADS_SUCCESS, EG_getBoundingBox(solid, box), "EG bounding box");
  diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
              (box[1] - box[4]) * (box[1] - box[4]) +
              (box[2] - box[5]) * (box[2] - box[5]));

  params[0] = 0.25 * diag;
  params[1] = 0.001 * diag;
  params[2] = 15.0;

  rebuild = REF_TRUE;
  while (rebuild) {
    REIS(EGADS_SUCCESS, EG_makeTessBody(solid, params, tess), "EG tess");
    REIS(EGADS_SUCCESS, EG_statusTessBody(*tess, &geom, &tess_status, &nvert),
         "EG tess");
    REIS(1, tess_status, "tess not closed");

    RSS(ref_egads_tess_adjust(ref_geom, *tess, &rebuild), "adjust params");
  }

  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_tess(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  ego tess;

  RSS(ref_egads_tess_create(ref_grid_geom(ref_grid), &tess),
      "create tess object");

  RSS(ref_egads_tess_fill_vertex(ref_grid, tess), "fill tess vertex");
  RSS(ref_egads_tess_fill_tri(ref_grid, tess), "fill tess triangles");
  RSS(ref_egads_tess_fill_edg(ref_grid, tess), "fill tess edges");

  RSS(ref_egads_mark_jump_degen(ref_grid), "T and UV jumps");
  ref_grid_surf(ref_grid) = REF_TRUE;

#else
  printf("returning empty grid from %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
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
