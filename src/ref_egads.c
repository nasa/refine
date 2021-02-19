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
#include "ref_dict.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"

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
  ref_geom->body = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_egads_out_level(REF_GEOM ref_geom, REF_INT out_level) {
#ifdef HAVE_EGADS
  ego context;
  context = (ego)(ref_geom->context);
  RAS(EG_setOutLevel(context, (int)out_level) >= 0, "set verbosity");
#else
  printf("nothing for %s, No EGADS linked\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(out_level);
#endif
  return REF_SUCCESS;
}

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_cache_body_objects(REF_GEOM ref_geom) {
  ego body = (ego)(ref_geom->body);
  ego *faces, *edges, *nodes;
  int nface, nedge, nnode;
  REF_INT face;
  int oclass, mtype, *senses, nchild;
  ego *children;
  ego ref;
  double uv_box[4];

#ifdef HAVE_EGADS_EFFECTIVE
  if (ref_geom_effective(ref_geom)) {
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes),
         "EG node topo");
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EEDGE, &nedge, &edges),
         "EG edge topo");
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EFACE, &nface, &faces),
         "EG face topo");
  } else {
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes),
         "EG node topo");
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges),
         "EG edge topo");
    REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, FACE, &nface, &faces),
         "EG face topo");
  }
#else
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes),
       "EG node topo");
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges),
       "EG edge topo");
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, FACE, &nface, &faces),
       "EG face topo");
#endif

  ref_geom->nnode = nnode;
  ref_geom->nodes = (void *)nodes;
  ref_geom->nedge = nedge;
  ref_geom->edges = (void *)edges;
  ref_geom->nface = nface;
  ref_geom->faces = (void *)faces;

  /* use face mtype SFORWARD, SREVERSE to set uv_area_sign */
  /* If it is SFORWARD (1) then the Face's Normal is in the same direction as
     the surface (u cross v), which points outward of the body. If it is
     SREVERSE (-1), then the natural surface normal points inward and the
     Face points consistently out of the body. */
  ref_malloc_init(ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 0.0);
  for (face = 0; face < nface; face++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->faces))[face], &ref, &oclass, &mtype,
                        uv_box, &nchild, &children, &senses),
         "topo");
    switch (mtype) {
      /* refine assumes normal point into the domain (body), flip sign */
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

  ref_malloc_init(ref_geom->face_min_length, ref_geom->nface, REF_DBL, -1.0);
  for (face = 0; face < nface; face++) {
    int len, atype;
    const double *preals;
    const int *pints;
    const char *string;
    if (EGADS_SUCCESS == EG_attributeRet(((ego *)(ref_geom->faces))[face],
                                         "min_length", &atype, &len, &pints,
                                         &preals, &string)) {
      if (ATTRREAL == atype && len == 1) {
        ref_geom->face_min_length[face] = preals[0];
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
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_load(REF_GEOM ref_geom, const char *filename) {
#ifdef HAVE_EGADS
  ego context;
  ego model = NULL;
  ego geom, *bodies, *children;
  int oclass, nego, mtype, nbody, *senses, nchild;
  ego body;

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
  if (NULL == filename) THROW("filename NULL for EGADS(full) load");
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

  ref_geom->model = (void *)model;

  REIS(EGADS_SUCCESS,
       EG_getTopology(model, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                      &senses),
       "EG topo bodies");
  REIS(1, nbody, "expected 1 body");
  body = bodies[0];
#ifdef HAVE_EGADS_EFFECTIVE
  {
    int ibody;
    int bodyclass, bodytype;
    ego owner, prev, next;
    for (ibody = 0; ibody < nego; ibody++) {
      REIS(EGADS_SUCCESS,
           EG_getInfo(bodies[ibody], &bodyclass, &bodytype, &owner, &prev,
                      &next),
           "info");
      if (EBODY == bodyclass) {
        body = bodies[ibody];
        ref_geom_effective(ref_geom) = REF_TRUE;
      }
    }
  }

#endif

  ref_geom->body = (void *)body;

  REIS(EGADS_SUCCESS,
       EG_getTopology(body, &geom, &oclass, &mtype, NULL, &nchild, &children,
                      &senses),
       "EG topo body type");
  RAB(SOLIDBODY == mtype || FACEBODY == mtype || SHEETBODY == mtype,
      "expected SOLIDBODY or FACEBODY or SHEETBODY",
      { printf("mtype %d\n", mtype); });
  ref_geom->manifold = (SOLIDBODY == mtype);

  RSS(ref_egads_cache_body_objects(ref_geom), "cache egads objects");

#else
  printf("nothing for %s, No EGADS linked for %s\n", __func__, filename);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_egads_save(REF_GEOM ref_geom, const char *filename) {
#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE)
  remove(filename); /* ignore failure */
  REIS(EGADS_SUCCESS, EG_saveModel((ego)(ref_geom->model), filename),
       "EG save");
#else
  printf("nothing for %s, No EGADS(full) linked for %s\n", __func__, filename);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif

  return REF_SUCCESS;
}

REF_BOOL ref_egads_allows_construction(void) {
#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE)
  return REF_TRUE;
#else
  return REF_FALSE;
#endif
}

REF_BOOL ref_egads_allows_effective(void) {
#if defined(HAVE_EGADS) && defined(HAVE_EGADS_EFFECTIVE) && \
    !defined(HAVE_EGADS_LITE)
  return REF_TRUE;
#else
  return REF_FALSE;
#endif
}

REF_STATUS ref_egads_construct(REF_GEOM ref_geom, const char *description) {
#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE)
  ego body = NULL;

  RAS(ref_egads_allows_construction(), "construction not allowed");

  if (0 == strcmp("cylinder", description)) {
    int stype = CYLINDER;
    double data[7] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    REIS(EGADS_SUCCESS,
         EG_makeSolidBody((ego)(ref_geom->context), stype, data, &body),
         "make solid body");
  }
  if (0 == strcmp("sphere", description)) {
    int stype = SPHERE;
    double data[4] = {0.0, 0.0, 0.0, 1.0};
    REIS(EGADS_SUCCESS,
         EG_makeSolidBody((ego)(ref_geom->context), stype, data, &body),
         "make solid body");
  }
  if (0 == strcmp("boxbox", description)) {
    int stype = BOX;
    ego box1, box2;
    {
      double data[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 0.5};
      REIS(EGADS_SUCCESS,
           EG_makeSolidBody((ego)(ref_geom->context), stype, data, &box1),
           "make solid body");
    }
    {
      double data[6] = {0.0, 0.0, 0.5, 1.0, 1.0, 0.5};
      REIS(EGADS_SUCCESS,
           EG_makeSolidBody((ego)(ref_geom->context), stype, data, &box2),
           "make solid body");
    }
    {
      ego boxbox, geom, *bodies;
      int oclass, nego, nbody, *senses;
      REIS(EGADS_SUCCESS, EG_generalBoolean(box1, box2, FUSION, 0.0, &boxbox),
           "make solid body");
      REIS(EGADS_SUCCESS,
           EG_getTopology(boxbox, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
      REIS(1, nbody, "expected 1 body");
      REIS(EGADS_SUCCESS, EG_copyObject(bodies[0], NULL, &body), "copy body");
      REIS(0, EG_deleteObject(boxbox), "delete temp model");
      REIS(0, EG_deleteObject(box1), "delete box1");
      REIS(0, EG_deleteObject(box2), "delete box2");
    }
  }
  if (0 == strcmp("steinmetz", description)) {
    int stype = CYLINDER;
    ego cyl1, cyl2;
    {
      double data[7] = {-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
      REIS(EGADS_SUCCESS,
           EG_makeSolidBody((ego)(ref_geom->context), stype, data, &cyl1),
           "make solid body");
    }
    {
      double data[7] = {0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
      REIS(EGADS_SUCCESS,
           EG_makeSolidBody((ego)(ref_geom->context), stype, data, &cyl2),
           "make solid body");
    }
    {
      ego stein, geom, *bodies;
      int oclass, nego, nbody, *senses;
      REIS(EGADS_SUCCESS,
           EG_generalBoolean(cyl1, cyl2, INTERSECTION, 0.0, &stein),
           "make solid body");
      REIS(EGADS_SUCCESS,
           EG_getTopology(stein, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
      REIS(1, nbody, "expected 1 body");
      REIS(EGADS_SUCCESS, EG_copyObject(bodies[0], NULL, &body), "copy body");
      REIS(0, EG_deleteObject(stein), "delete temp model");
      REIS(0, EG_deleteObject(cyl1), "delete cyl1");
      REIS(0, EG_deleteObject(cyl2), "delete cyl2");
    }
  }
  if (0 == strcmp("square", description)) {
    ego nodes[4], curves[4], edges[4];
    { /* nodes */
      double xyz[3];
      xyz[0] = 0;
      xyz[1] = 0;
      xyz[2] = 0;
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, NODE, 0, xyz, 0,
                           NULL, NULL, &nodes[0]),
           "make node");
      xyz[0] = 1;
      xyz[1] = 0;
      xyz[2] = 0;
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, NODE, 0, xyz, 0,
                           NULL, NULL, &nodes[1]),
           "make node");
      xyz[0] = 1;
      xyz[1] = 1;
      xyz[2] = 0;
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, NODE, 0, xyz, 0,
                           NULL, NULL, &nodes[2]),
           "make node");
      xyz[0] = 0;
      xyz[1] = 1;
      xyz[2] = 0;
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, NODE, 0, xyz, 0,
                           NULL, NULL, &nodes[3]),
           "make node");
    }

    { /* curves */
      ego ref, *pchldrn;
      int oclass, mtype, nchild, *psens;
      double data[6];
      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[0], &ref, &oclass, &mtype, &(data[0]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[1], &ref, &oclass, &mtype, &(data[3]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      data[3] -= data[0];
      data[4] -= data[1];
      data[5] -= data[2];
      REIS(EGADS_SUCCESS,
           EG_makeGeometry((ego)(ref_geom->context), CURVE, LINE, NULL, NULL,
                           data, &curves[0]),
           "line");

      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[1], &ref, &oclass, &mtype, &(data[0]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[2], &ref, &oclass, &mtype, &(data[3]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      data[3] -= data[0];
      data[4] -= data[1];
      data[5] -= data[2];
      REIS(EGADS_SUCCESS,
           EG_makeGeometry((ego)(ref_geom->context), CURVE, LINE, NULL, NULL,
                           data, &curves[1]),
           "line");

      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[2], &ref, &oclass, &mtype, &(data[0]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[3], &ref, &oclass, &mtype, &(data[3]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      data[3] -= data[0];
      data[4] -= data[1];
      data[5] -= data[2];
      REIS(EGADS_SUCCESS,
           EG_makeGeometry((ego)(ref_geom->context), CURVE, LINE, NULL, NULL,
                           data, &curves[2]),
           "line");

      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[3], &ref, &oclass, &mtype, &(data[0]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS,
           EG_getTopology(nodes[0], &ref, &oclass, &mtype, &(data[3]), &nchild,
                          &pchldrn, &psens),
           "node xyz");
      data[3] -= data[0];
      data[4] -= data[1];
      data[5] -= data[2];
      REIS(EGADS_SUCCESS,
           EG_makeGeometry((ego)(ref_geom->context), CURVE, LINE, NULL, NULL,
                           data, &curves[3]),
           "line");
    }

    { /* edges */
      ego ref, *pchldrn, objs[2];
      int oclass, mtype, nchild, *psens;
      double xyz[3], dum[3], range[2];
      objs[0] = nodes[0];
      objs[1] = nodes[1];
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[0], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[0], xyz, &range[0], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[1], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[0], xyz, &range[1], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), curves[0], EDGE, TWONODE,
                           range, 2, objs, NULL, &edges[0]),
           "make edge");

      objs[0] = nodes[1];
      objs[1] = nodes[2];
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[0], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[1], xyz, &range[0], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[1], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[1], xyz, &range[1], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), curves[1], EDGE, TWONODE,
                           range, 2, objs, NULL, &edges[1]),
           "make edge");

      objs[0] = nodes[2];
      objs[1] = nodes[3];
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[0], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[2], xyz, &range[0], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[1], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[2], xyz, &range[1], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), curves[2], EDGE, TWONODE,
                           range, 2, objs, NULL, &edges[2]),
           "make edge");

      objs[0] = nodes[3];
      objs[1] = nodes[0];
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[0], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[3], xyz, &range[0], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_getTopology(objs[1], &ref, &oclass, &mtype, xyz, &nchild,
                          &pchldrn, &psens),
           "node xyz");
      REIS(EGADS_SUCCESS, EG_invEvaluate(curves[3], xyz, &range[1], dum),
           "trange");
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), curves[3], EDGE, TWONODE,
                           range, 2, objs, NULL, &edges[3]),
           "make edge");
    }

    { /* loop and face */
      int senses[4] = {1, 1, 1, 1};
      ego loop, face;
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, LOOP, CLOSED, NULL,
                           4, edges, senses, &loop),
           "make loop");
      REIS(EGADS_SUCCESS, EG_makeFace(loop, SREVERSE, NULL, &face), "face");
      REIS(EGADS_SUCCESS,
           EG_makeTopology((ego)(ref_geom->context), NULL, BODY, FACEBODY, NULL,
                           1, &face, NULL, &body),
           "facebody");
      REIS(0, EG_deleteObject((ego)(ref_geom->context)),
           "delete construction objects");
    }
  }
  if (0 == strcmp("revolve", description)) {
    ego face;
    double data[6];
    RSS(ref_egads_construct(ref_geom, "square"), "create");
    REIS(EGADS_SUCCESS, EG_copyObject((ego)(ref_geom->body), NULL, &face),
         "copy body");
    REIS(0, EG_deleteObject((ego)(ref_geom->model)), "delete body model");

    data[0] = 0.0;
    data[1] = 0.0;
    data[2] = 0.0;
    data[3] = 1.0;
    data[4] = 0.0;
    data[5] = 0.0;
    REIS(EGADS_SUCCESS, EG_rotate(face, 360.0, data, &body), "revolve");
    REIS(0, EG_deleteObject((ego)(ref_geom->context)),
         "delete construction objects");
  }
  RNB(body, "unknown description", { printf(">%s<\n", description); });
  {
    ego model;
    REIS(EGADS_SUCCESS,
         EG_makeTopology((ego)(ref_geom->context), NULL, MODEL, 0, NULL, 1,
                         &body, NULL, &model),
         "make Topo Model");
    ref_geom->model = (void *)model;
  }
  ref_geom->body = (void *)body;
  {
    ego geom, *children;
    int oclass, mtype, *senses, nchild;
    REIS(EGADS_SUCCESS,
         EG_getTopology(body, &geom, &oclass, &mtype, NULL, &nchild, &children,
                        &senses),
         "EG topo body type");
    ref_geom->manifold = (SOLIDBODY == mtype);
  }

  RSS(ref_egads_cache_body_objects(ref_geom), "cache egads objects");

#else
  printf("nothing for %s, Full EGADS not linked for %s\n", __func__,
         description);
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
static REF_STATUS ref_egads_edge_faces(REF_GEOM ref_geom,
                                       REF_INT **edge_face_arg) {
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
        edge = EG_indexBodyTopo((ego)(ref_geom->body), eedges[iedge]) - 1;
        RAB(2 > nface[edge], "edge has more than 2 faces",
            printf("face ids %d %d %d edge id %d\n", e2f[0 + 2 * edge],
                   e2f[1 + 2 * edge], face + 1, edge + 1));
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
static REF_STATUS ref_egads_node_faces(REF_GEOM ref_geom,
                                       REF_ADJ *ref_adj_arg) {
  REF_ADJ ref_adj;
  REF_INT *e2f, id, toponode;
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double trange[2];
  RSS(ref_adj_create(ref_adj_arg), "create ref_adj");
  ref_adj = *ref_adj_arg;
  RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");
  for (id = 1; id <= ref_geom->nedge; id++) {
    object = ((ego *)(ref_geom->edges))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");

    if (0 < nchild) {
      toponode = EG_indexBodyTopo(ref_geom->body, pchldrn[0]);
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
      toponode = EG_indexBodyTopo(ref_geom->body, pchldrn[1]);
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
static REF_STATUS ref_egads_tess_fill_tri(REF_GRID ref_grid, ego tess,
                                          REF_GLOB n_global) {
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
      RAB(0 < nodes[0] && nodes[0] <= n_global, "tri global out of range", {
        printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[0]);
      });
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
      RAB(0 < nodes[1] && nodes[1] <= n_global, "tri n1 global out of range", {
        printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[0]);
      });
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[1 + 3 * tri], &(nodes[0])),
           "l2g1");
      RAB(0 < nodes[0] && nodes[0] <= n_global, "tri n0 global out of range", {
        printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[0]);
      });
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[2 + 3 * tri], &(nodes[2])),
           "l2g2");
      RAB(0 < nodes[2] && nodes[2] <= n_global, "tri n2 global out of range", {
        printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[2]);
      });
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
static REF_STATUS ref_egads_tess_fill_edg(REF_GRID ref_grid, ego tess,
                                          REF_GLOB n_global) {
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
        RAB(0 < nodes[0] && nodes[0] <= n_global, "edg global out of range", {
          printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[0]);
        });
        nodes[0] -= 1;
        param[0] = t[node];
        RSB(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, edge + 1, param),
            "edge t", {
              printf("edge %d of %d plen %d node %d global %d\n", edge + 1,
                     ref_geom->nedge, plen, node + 1, nodes[0] + 1);
            });
      }
    }
    if (!degenerate)
      for (node = 0; node < (plen - 1); node++) {
        /* assume edge index is 1-bias */
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0])),
             "l2g0");
        RAB(0 < nodes[0] && nodes[0] <= n_global, "edg0 global out of range", {
          printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[0]);
        });
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 2, &(nodes[1])),
             "l2g1");
        RAB(0 < nodes[1] && nodes[1] <= n_global, "edg1 global out of range", {
          printf("nvert " REF_GLOB_FMT " global %d\n", n_global, nodes[1]);
        });
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
                                       REF_CLOUD edge_tp_augment, REF_INT *e2f,
                                       REF_LIST face_locked) {
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
  REF_BOOL contains0, contains1;

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
      ineligible_cad_node0 = EG_indexBodyTopo(ref_geom->body, cadnodes0[0]);
      if (2 == ncadnode0) {
        ineligible_cad_node1 = EG_indexBodyTopo(ref_geom->body, cadnodes0[1]);
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
          cad_node0 = EG_indexBodyTopo(ref_geom->body, cadnodes1[0]);
          if (2 == ncadnode1) {
            cad_node1 = EG_indexBodyTopo(ref_geom->body, cadnodes1[1]);
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
      edgeid = EG_indexBodyTopo(ref_geom->body, edgeobj0);
      if (ref_math_divisible(diag, width)) {
        aspect_ratio = diag / width;
        adjusted = MIN(MAX(1.0, aspect_ratio - 10.0), 10.0) * width;
        params[0] = 2.0 * adjusted;
        params[1] = 0.3 * params[0];
        params[2] = 20.0;
        RSS(ref_list_contains(face_locked, e2f[0 + 2 * (edgeid - 1)],
                              &contains0),
            "lock face0");
        RSS(ref_list_contains(face_locked, e2f[1 + 2 * (edgeid - 1)],
                              &contains1),
            "lock face1");
        if (!contains0 && !contains1)
          RSS(ref_egads_merge_tparams(edge_tp_augment, edgeid, params),
              "update tparams");
      }
    }
  }

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_adjust_tparams_missing_face(
    REF_GEOM ref_geom, ego tess, REF_CLOUD face_tp_augment,
    REF_CLOUD edge_tp_augment, REF_DBL seg_per_diag, REF_LIST face_locked) {
  ego faceobj;
  int face, tlen, plen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;

  double params[3], diag, box[6];
  REF_INT edge, *e2f;
  REF_BOOL contains;

  RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");

  for (face = 0; face < (ref_geom->nface); face++) {
    faceobj = ((ego *)(ref_geom->faces))[face];
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    if (0 == plen || 0 == tlen) {
      printf("face %d has %d nodes and %d triangles\n", face + 1, plen, tlen);
      REIS(EGADS_SUCCESS, EG_getBoundingBox(faceobj, box), "EG bounding box");
      diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                  (box[1] - box[4]) * (box[1] - box[4]) +
                  (box[2] - box[5]) * (box[2] - box[5]));
      params[0] = seg_per_diag * diag;
      params[1] = 0.1 * params[0];
      params[2] = 15.0;
      RSS(ref_list_contains(face_locked, face + 1, &contains), "lock face");
      if (!contains) RSS(ref_list_push(face_locked, face + 1), "lock face");
      RSS(ref_egads_merge_tparams(face_tp_augment, face + 1, params),
          "update tparams");
      for (edge = 0; edge < (ref_geom->nedge); edge++) {
        if (face + 1 == e2f[0 + 2 * edge] || face + 1 == e2f[0 + 2 * edge]) {
          RSS(ref_egads_merge_tparams(edge_tp_augment, edge + 1, params),
              "update tparams");
        }
      }
    }
  }

  ref_free(e2f);

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_adjust_tparams_dup_edge(REF_GEOM ref_geom, ego tess,
                                                    REF_CLOUD edge_tp_augment,
                                                    REF_DBL seg_per_diag) {
  ego edgeobj;
  int edge, plen;
  const double *points, *t;
  ego ref;
  int oclass, mtype;
  double trange[2];
  int ncadnode;
  ego *cadnodes;
  int *senses;

  REF_CELL ref_cell;
  REF_INT node, new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, ncell, cell_list[2];
  double params[3], diag, box[6];

  RSS(ref_cell_create(&ref_cell, REF_CELL_EDG), "temp edg create");

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    edgeobj = ((ego *)(ref_geom->edges))[edge];
    REIS(EGADS_SUCCESS,
         EG_getTopology(edgeobj, &ref, &oclass, &mtype, trange, &ncadnode,
                        &cadnodes, &senses),
         "EG topo edge0");
    if (mtype == DEGENERATE) continue; /* skip DEGENERATE */
    REIS(EGADS_SUCCESS, EG_getTessEdge(tess, edge + 1, &plen, &points, &t),
         "tess edge nodes");
    for (node = 0; node < (plen - 1); node++) {
      /* assume edge index is 1-bias */
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0])), "l2g0");
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, -(edge + 1), node + 2, &(nodes[1])), "l2g1");
      nodes[0] -= 1;
      nodes[1] -= 1;
      nodes[2] = edge + 1;
      RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new edge");
    }
  }

  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_list_with2(ref_cell, ref_cell_c2n(ref_cell, 0, cell),
                            ref_cell_c2n(ref_cell, 1, cell), 2, &ncell,
                            cell_list),
        "edge list for edge");
    if (2 == ncell) {
      REF_INT id;
      printf("edge ids %d and %d share a segment\n",
             ref_cell_c2n(ref_cell, 2, cell_list[0]),
             ref_cell_c2n(ref_cell, 2, cell_list[1]));
      id = ref_cell_c2n(ref_cell, 2, cell_list[0]);
      edge = id - 1;
      edgeobj = ((ego *)(ref_geom->edges))[edge];
      REIS(EGADS_SUCCESS, EG_getBoundingBox(edgeobj, box), "EG bounding box");
      diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                  (box[1] - box[4]) * (box[1] - box[4]) +
                  (box[2] - box[5]) * (box[2] - box[5]));
      params[0] = seg_per_diag * diag;
      params[1] = 0.1 * params[0];
      params[2] = 15.0;
      RSS(ref_egads_merge_tparams(edge_tp_augment, edge + 1, params),
          "update tparams");
      id = ref_cell_c2n(ref_cell, 2, cell_list[1]);
      edge = id - 1;
      edgeobj = ((ego *)(ref_geom->edges))[edge];
      REIS(EGADS_SUCCESS, EG_getBoundingBox(edgeobj, box), "EG bounding box");
      diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                  (box[1] - box[4]) * (box[1] - box[4]) +
                  (box[2] - box[5]) * (box[2] - box[5]));
      params[0] = seg_per_diag * diag;
      params[1] = 0.1 * params[0];
      params[2] = 15.0;
      RSS(ref_egads_merge_tparams(edge_tp_augment, edge + 1, params),
          "update tparams");
    }
  }

  RSS(ref_cell_free(ref_cell), "free temp edg");

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_adjust_tparams_single_edge(
    REF_GEOM ref_geom, ego tess, REF_CLOUD edge_tp_augment,
    REF_INT auto_tparams, REF_LIST face_locked) {
  ego edgeobj;
  int edge, plen;
  const double *points, *t;
  ego ref;
  int oclass, mtype;
  double trange[2];
  int ncadnode;
  ego *cadnodes;
  int *senses;

  REF_INT *e2f;

  REF_INT i;
  REF_DBL t_mid, xyz[3], dside[3], dmid[3];
  const REF_DBL *xyz0, *xyz1;
  REF_DBL length, offset, chord;

  double params[3];
  REF_DBL chord_limit = 0.49; /* semicircle has an chord of 0.5 */
  REF_BOOL contains0, contains1;

  RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    edgeobj = ((ego *)(ref_geom->edges))[edge];
    REIS(EGADS_SUCCESS,
         EG_getTopology(edgeobj, &ref, &oclass, &mtype, trange, &ncadnode,
                        &cadnodes, &senses),
         "EG topo edge0");
    if (mtype == DEGENERATE) continue; /* skip DEGENERATE */
    REIS(EGADS_SUCCESS, EG_getTessEdge(tess, edge + 1, &plen, &points, &t),
         "tess edge nodes");
    if (2 == plen) {
      xyz0 = &(points[3 * 0]);
      xyz1 = &(points[3 * 1]);
      t_mid = 0.5 * (t[0] + t[1]);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, edge + 1, &t_mid, xyz,
                            NULL),
          "eval mid uv");
      for (i = 0; i < 3; i++) dside[i] = xyz1[i] - xyz0[i];
      for (i = 0; i < 3; i++) dmid[i] = xyz[i] - 0.5 * (xyz1[i] + xyz0[i]);
      length =
          sqrt(dside[0] * dside[0] + dside[1] * dside[1] + dside[2] * dside[2]);
      offset = sqrt(dmid[0] * dmid[0] + dmid[1] * dmid[1] + dmid[2] * dmid[2]);
      chord = 0.0;
      if (ref_math_divisible(offset, length)) {
        chord = offset / length;
      }
      if (chord > chord_limit)
        printf("edge id %d single segment with chord %.2e\n", edge + 1, chord);
      if (chord > chord_limit &&
          (auto_tparams & REF_EGADS_SINGLE_EDGE_TPARAM)) {
        params[0] = 0.4 * length;
        params[1] = 0.1 * length;
        params[2] = 15.0;
        RSS(ref_list_contains(face_locked, e2f[0 + 2 * edge], &contains0),
            "lock face0");
        RSS(ref_list_contains(face_locked, e2f[1 + 2 * edge], &contains1),
            "lock face1");
        if (!contains0 && !contains1) {
          RSS(ref_egads_merge_tparams(edge_tp_augment, edge + 1, params),
              "update tparams");
        }
      }
    }
  }

  ref_free(e2f);

  return REF_SUCCESS;
}
#endif

#ifdef HAVE_EGADS
static REF_STATUS ref_egads_adjust_tparams_chord(REF_GEOM ref_geom, ego tess,
                                                 REF_CLOUD face_tp_augment,
                                                 REF_CLOUD edge_tp_augment,
                                                 REF_INT auto_tparams,
                                                 REF_LIST face_locked) {
  ego faceobj;
  int face, tlen, plen;
  const double *points, *uv;
  const int *ptype, *pindex, *tris, *tric;

  double params[3], diag, box[6];
  REF_INT edge, *e2f;

  REF_INT tri, side, n0, n1, i;
  const REF_DBL *xyz0, *xyz1, *uv0, *uv1;
  REF_DBL uvm[2], xyz[3], dside[3], dmid[3];
  REF_DBL length, offset, max_chord, max_chord_length, max_chord_offset;
  REF_BOOL contains, contains0, contains1;

  RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");

  for (face = 0; face < (ref_geom->nface); face++) {
    faceobj = ((ego *)(ref_geom->faces))[face];
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    RAB(0 < plen && 0 < tlen, "missing face",
        { printf("face id %d plen %d tlen %d\n", face + 1, plen, tlen); });

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
        RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, face + 1, uvm, xyz,
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
    if (max_chord > 0.2) {
      REIS(EGADS_SUCCESS, EG_getBoundingBox(faceobj, box), "EG bounding box");
      diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                  (box[1] - box[4]) * (box[1] - box[4]) +
                  (box[2] - box[5]) * (box[2] - box[5]));
      printf("face %d rel chord %f abs len %f abs chord %f diag %f\n", face + 1,
             max_chord, max_chord_length, max_chord_offset, diag);
      if (auto_tparams & REF_EGADS_CHORD_TPARAM) {
        params[0] = 0.25 * diag;
        params[1] = 0.025 * diag;
        params[2] = 15.0;
        RSS(ref_list_contains(face_locked, face + 1, &contains), "lock face");
        if (!contains)
          RSS(ref_egads_merge_tparams(face_tp_augment, face + 1, params),
              "update tparams");
        for (edge = 0; edge < (ref_geom->nedge); edge++) {
          if (face + 1 == e2f[0 + 2 * edge] || face + 1 == e2f[0 + 2 * edge]) {
            RSS(ref_list_contains(face_locked, e2f[0 + 2 * edge], &contains0),
                "lock face0");
            RSS(ref_list_contains(face_locked, e2f[1 + 2 * edge], &contains1),
                "lock face1");
            if (!contains0 && !contains1)
              RSS(ref_egads_merge_tparams(edge_tp_augment, edge + 1, params),
                  "update tparams");
          }
        }
      }
    }

    /* face width parameter to all edges */
    if (auto_tparams & REF_EGADS_WIDTH_TPARAM) {
      RSS(ref_egads_face_width(ref_geom, face + 1, edge_tp_augment, e2f,
                               face_locked),
          "face width");
    }
  }

  ref_free(e2f);

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
        /* examine existing .tParams and detect change within tol */
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
        /* examine existing .tParams and detect change within tol */
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
                                        REF_INT auto_tparams,
                                        REF_DBL *global_params) {
  ego body, geom;
  int tess_status, nvert;
  double params[3], diag, box[6];
  REF_BOOL rebuild;
  REF_INT tries;
  REF_CLOUD face_tp_original, face_tp_augment;
  REF_CLOUD edge_tp_original, edge_tp_augment;
  REF_LIST face_locked;
  REF_DBL seg_per_diag;
  REF_INT face;

  body = (ego)(ref_geom->body);
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

  RSS(ref_list_create(&face_locked), "create locked face list");

  REIS(EGADS_SUCCESS, EG_getBoundingBox(body, box), "EG bounding box");
  diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
              (box[1] - box[4]) * (box[1] - box[4]) +
              (box[2] - box[5]) * (box[2] - box[5]));

  if (NULL != global_params) {
    params[0] = global_params[0];
    params[1] = global_params[1];
    params[2] = global_params[2];
  } else {
    params[0] = 0.025 * diag;
    params[1] = 0.0075 * diag;
    params[2] = 20.0;
  }

  rebuild = REF_TRUE;
  tries = 0;
  seg_per_diag = 2.0;
  while (rebuild) {
    seg_per_diag *= 0.5;
    tries++;
    RAS(tries < 10, "exhausted tries");
    printf("makeTessBody global param ( %e %e %f)\n", params[0], params[1],
           params[2]);
    REIS(EGADS_SUCCESS, EG_makeTessBody(body, params, tess), "EG tess");
    REIS(EGADS_SUCCESS, EG_statusTessBody(*tess, &geom, &tess_status, &nvert),
         "EG tess");
    REIS(1, tess_status, "tess not closed");

    RSS(ref_egads_adjust_tparams_missing_face(ref_geom, *tess, face_tp_augment,
                                              edge_tp_augment, seg_per_diag,
                                              face_locked),
        "adjust missing face params");
    RSS(ref_egads_adjust_tparams_dup_edge(ref_geom, *tess, edge_tp_augment,
                                          seg_per_diag),
        "adjust dup edge params");
    rebuild = REF_FALSE;
    RSS(ref_egads_update_tparams_attributes(ref_geom, face_tp_original,
                                            edge_tp_original, face_tp_augment,
                                            edge_tp_augment, &rebuild),
        "update tparams mandatory");

    RSS(ref_list_inspect(face_locked), "show");

    if (rebuild) {
      printf(
          "rebuild EGADS tessellation after missing face .tParams adjustment, "
          "try %d\n",
          tries);
      REIS(0, EG_deleteObject(*tess), "delete previous try at tess");
    }
  }

  for (face = 0; face < (ref_geom->nface); face++) {
    int tlen, plen;
    const double *points, *uv;
    const int *ptype, *pindex, *tris, *tric;
    REIS(EGADS_SUCCESS,
         EG_getTessFace(*tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    RAB(0 < plen && 0 < tlen, "missing face",
        { printf("face id %d plen %d tlen %d\n", face + 1, plen, tlen); });
  }

  rebuild = REF_TRUE;
  tries = 0;
  while (rebuild) {
    RSS(ref_egads_adjust_tparams_single_edge(ref_geom, *tess, edge_tp_augment,
                                             auto_tparams, face_locked),
        "adjust single edge params");
    RSS(ref_egads_adjust_tparams_chord(ref_geom, *tess, face_tp_augment,
                                       edge_tp_augment, auto_tparams,
                                       face_locked),
        "adjust chord params");
    rebuild = REF_FALSE;
    RSS(ref_egads_update_tparams_attributes(ref_geom, face_tp_original,
                                            edge_tp_original, face_tp_augment,
                                            edge_tp_augment, &rebuild),
        "update auto tparams");

    if (rebuild) {
      tries++;
      RAS(tries < 5, "exhausted tries");
      printf(
          "rebuild EGADS tessellation after chord .tParams adjustment, "
          "try %d\n",
          tries);
      REIS(0, EG_deleteObject(*tess), "delete previous try at tess");
      REIS(EGADS_SUCCESS, EG_makeTessBody(body, params, tess), "EG tess");
      REIS(EGADS_SUCCESS, EG_statusTessBody(*tess, &geom, &tess_status, &nvert),
           "EG tess");
      REIS(1, tess_status, "tess not closed");
    }
  }

  RSS(ref_list_free(face_locked), "free face list");

  RSS(ref_cloud_free(edge_tp_augment), "free tparams augment");
  RSS(ref_cloud_free(face_tp_augment), "free tparams augment");
  RSS(ref_cloud_free(edge_tp_original), "free tparams cache");
  RSS(ref_cloud_free(face_tp_original), "free tparams cache");

  for (face = 0; face < (ref_geom->nface); face++) {
    int tlen, plen;
    const double *points, *uv;
    const int *ptype, *pindex, *tris, *tric;
    REIS(EGADS_SUCCESS,
         EG_getTessFace(*tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");
    RAB(0 < plen && 0 < tlen, "missing face",
        { printf("face id %d plen %d tlen %d\n", face + 1, plen, tlen); });
  }

  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_tess(REF_GRID ref_grid, REF_INT auto_tparams,
                          REF_DBL *global_params) {
#ifdef HAVE_EGADS
  ego tess;
  REF_GLOB n_global;

  if (ref_mpi_once(ref_grid_mpi(ref_grid))) {
    RSS(ref_egads_tess_create(ref_grid_geom(ref_grid), &tess, auto_tparams,
                              global_params),
        "create tess object");

    RSS(ref_egads_tess_fill_vertex(ref_grid, tess, &n_global),
        "fill tess vertex");
    RSS(ref_egads_tess_fill_tri(ref_grid, tess, n_global),
        "fill tess triangles");
    RSS(ref_egads_tess_fill_edg(ref_grid, tess, n_global), "fill tess edges");
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
  SUPRESS_UNUSED_COMPILER_WARNING(global_params);
#endif

  return REF_SUCCESS;
}

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
      cad_node = EG_indexBodyTopo(ref_geom->body, echilds[0]);
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

  RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->edges))[edge], &eref, &oclass,
                        &mtype, trange, &nchild, &echilds, &senses),
         "edge topo");
    if (mtype == DEGENERATE) {
      geom_node_id = EG_indexBodyTopo((ego)(ref_geom->body), echilds[0]);
      face = e2f[0 + 2 * edge] - 1;
      REIS(REF_EMPTY, e2f[1 + 2 * edge], "DEGENERATE edge has two faces");

      face_geom = REF_EMPTY;
      each_ref_geom_node(ref_geom, node_geom) {
        if (geom_node_id == ref_geom_id(ref_geom, node_geom)) {
          node = ref_geom_node(ref_geom, node_geom);
          RSB(ref_geom_find(ref_geom, node, REF_GEOM_FACE, face + 1,
                            &face_geom),
              "face for degen edge at node not found",
              { printf("edgeid %d faceid %d\n", edge + 1, face + 1); });
        }
      }

      /* in parallel, may not have this node */
      if (REF_EMPTY != face_geom) {
        uv[0] = ref_geom_param(ref_geom, 0, face_geom);
        uv[1] = ref_geom_param(ref_geom, 1, face_geom);
        RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, face + 1, uv, xyz,
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
  RSS(ref_egads_node_faces(ref_geom, &n2f), "build n2f");
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

  RSS(ref_egads_edge_faces(ref_geom, &e2f), "compute edge faces");

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
      toponode0 = EG_indexBodyTopo(ref_geom->body, pchldrn[0]);
      toponode1 = EG_indexBodyTopo(ref_geom->body, pchldrn[1]);
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
        RSS(ref_egads_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
                                   ref_node_xyz_ptr(ref_node, best_node),
                                   param),
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
            continue; /* this candidate is already part of the edge */
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
          RSS(ref_egads_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
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

  printf("setting face UV under edges\n");
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
          RSS(ref_egads_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
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
          RSS(ref_egads_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
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
          RSS(ref_egads_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
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

REF_STATUS ref_egads_diagonal(REF_GEOM ref_geom, REF_INT geom, REF_DBL *diag) {
#ifdef HAVE_EGADS
  ego object;
  double box[6];

  object = (ego)(ref_geom->body);
  RNS(object, "EGADS body object is NULL. Has the geometry been loaded?");

  if (geom < 0) {
    object = (ego)(ref_geom->body);
  } else {
    switch (ref_geom_type(ref_geom, geom)) {
      case REF_GEOM_EDGE:
        object = ((ego *)(ref_geom->edges))[ref_geom_id(ref_geom, geom) - 1];
        break;
      case REF_GEOM_FACE:
        object = ((ego *)(ref_geom->faces))[ref_geom_id(ref_geom, geom) - 1];
        break;
      default:
        *diag = 0.0; /* for node */
        return REF_SUCCESS;
    }
  }

  REIS(EGADS_SUCCESS, EG_getBoundingBox(object, box), "EG bounding box");
  *diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
               (box[1] - box[4]) * (box[1] - box[4]) +
               (box[2] - box[5]) * (box[2] - box[5]));

#else
  printf("returning 1.0 from %s, No EGADS\n", __func__);
  *diag = 1.0;
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(geom);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_egads_tolerance(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                               REF_DBL *tolerance) {
#ifdef HAVE_EGADS
  ego object, *objects;
  double tol;
  int egads_status;

  *tolerance = 0.0;

  object = (ego)NULL;
  switch (type) {
    case REF_GEOM_NODE:
      if (id < 1 || id > ref_geom->nnode) return REF_INVALID;
      objects = (ego *)(ref_geom->nodes);
      object = objects[id - 1];
      break;
    case REF_GEOM_EDGE:
      if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
      objects = (ego *)(ref_geom->edges);
      object = objects[id - 1];
      break;
    case REF_GEOM_FACE:
      if (id < 1 || id > ref_geom->nface) return REF_INVALID;
      objects = (ego *)(ref_geom->faces);
      object = objects[id - 1];
      break;
    case REF_GEOM_SOLID:
      object = (ego)(ref_geom->body);
      break;
    default:
      printf("ref_geom type %d unknown\n", type);
      RSS(REF_IMPLEMENT, "unknown surface type");
  }

  egads_status = EG_getTolerance(object, &tol);
  if (EGADS_SUCCESS == egads_status) {
    *tolerance = tol;
  }
  if (EGADS_NOTTOPO == egads_status) {
    /* EFFECTIVE */
    *tolerance = 0.0;
    return REF_SUCCESS;
  }
  REIS(EGADS_SUCCESS, egads_status, "EG tolerance");

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
  *tolerance = -1.0;
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_egads_edge_curvature(REF_GEOM ref_geom, REF_INT geom, REF_DBL *k,
                                    REF_DBL *normal) {
#ifdef HAVE_EGADS
  double curvature[4];
  ego *edges;
  ego object;
  int edgeid;
  double t;
  int egads_status;
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  REIS(REF_GEOM_EDGE, ref_geom_type(ref_geom, geom), "expected edge geom");
  RNS(ref_geom->edges, "edges not loaded");
  edgeid = ref_geom_id(ref_geom, geom);
  edges = (ego *)(ref_geom->edges);
  object = edges[edgeid - 1];
  RNS(object, "EGADS object is NULL. Has the geometry been loaded?");

  t = ref_geom_param(ref_geom, 0, geom); /* ignores periodic */

#if !defined(HAVE_EGADS_LITE) && defined(HAVE_EGADS_EFFECTIVE)
  if (ref_geom_effective(ref_geom) && !ref_geom_effective_curvature(ref_geom)) {
    ego underlying_object;
    double underlying_t;
    REIS(EGADS_SUCCESS,
         EG_effectiveMap(object, &t, &underlying_object, &underlying_t),
         "map effective to brep");
    object = underlying_object;
    t = underlying_t;
  }
#endif

  egads_status = EG_curvature(object, &t, curvature);
  if (EGADS_DEGEN == egads_status) {
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double t_range[2];
    double t_offset;
    REF_DBL shift = 1.0e-2;
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, t_range, &nchild,
                        &pchldrn, &psens),
         "EG topo face");
    t_offset = (1.0 - shift) * t + shift * 0.5 * (t_range[0] + t_range[1]);
    egads_status = EG_curvature(object, &t_offset, curvature);
  }
  if (EGADS_SUCCESS == egads_status) {
    *k = curvature[0];
    normal[0] = curvature[1];
    normal[1] = curvature[2];
    normal[2] = curvature[3];
    return REF_SUCCESS;
  }
  if (EGADS_NOTGEOM == egads_status) {
    /* EFFECTIVE */
    *k = 0;
    normal[0] = 1;
    normal[1] = 0;
    normal[2] = 0;
    return REF_SUCCESS;
  }
  printf("EG_curvature %d (-24 is DEGEN) edgeid %d t %e\n", egads_status,
         edgeid, t);
  *k = 0;
  normal[0] = 1;
  normal[1] = 0;
  normal[2] = 0;
  return REF_FAILURE;
#else
  printf("curvature 0: No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(geom);
  *k = 0.0;
  normal[0] = 1.0;
  normal[1] = 0.0;
  normal[2] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_face_curvature(REF_GEOM ref_geom, REF_INT geom,
                                    REF_DBL *kr, REF_DBL *r, REF_DBL *ks,
                                    REF_DBL *s) {
#ifdef HAVE_EGADS
  REF_DBL uv[2];
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  REIS(REF_GEOM_FACE, ref_geom_type(ref_geom, geom), "expected face geom");

  uv[0] = ref_geom_param(ref_geom, 0, geom); /* ignores periodic */
  uv[1] = ref_geom_param(ref_geom, 1, geom);
  RAISE(ref_egads_face_curvature_at(ref_geom, ref_geom_id(ref_geom, geom),
                                    ref_geom_degen(ref_geom, geom), uv, kr, r,
                                    ks, s));
  return REF_SUCCESS;
#else
  printf("curvature 0, 0: No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(geom);
  *kr = 0.0;
  r[0] = 1.0;
  r[1] = 0.0;
  r[2] = 0.0;
  *ks = 0.0;
  s[0] = 0.0;
  s[1] = 1.0;
  s[2] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_face_curvature_at(REF_GEOM ref_geom, REF_INT faceid,
                                       REF_INT degen, REF_DBL *uv, REF_DBL *kr,
                                       REF_DBL *r, REF_DBL *ks, REF_DBL *s) {
#ifdef HAVE_EGADS
  double curvature[8];
  ego *faces;
  ego object;
  int egads_status;
  RNS(ref_geom->faces, "faces not loaded");
  faces = (ego *)(ref_geom->faces);
  RAS(0 < faceid && faceid <= ref_geom->nface, "face out of range");
  object = faces[faceid - 1];
  RNS(object, "EGADS object is NULL. Has the geometry been loaded?");

#if !defined(HAVE_EGADS_LITE) && defined(HAVE_EGADS_EFFECTIVE)
  if (ref_geom_effective(ref_geom) && !ref_geom_effective_curvature(ref_geom)) {
    ego underlying_object;
    double underlying_uv[2];
    REIS(EGADS_SUCCESS,
         EG_effectiveMap(object, uv, &underlying_object, underlying_uv),
         "map effective to brep");
    object = underlying_object;
    uv[0] = underlying_uv[0];
    uv[1] = underlying_uv[1];
  }
#endif

  egads_status = EG_curvature(object, uv, curvature);
  /* classic marked degen where u or v collapses to a point, move tangent */
  if (0 != degen || EGADS_DEGEN == egads_status) {
    REF_DBL du, dv;
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double uv_range[4];
    double params[2];
    double eval[18];
    REF_DBL shift = 1.0e-2;
    params[0] = uv[0];
    params[1] = uv[1];
    REIS(EGADS_SUCCESS, EG_evaluate(object, params, eval), "eval derivs");
    du = sqrt(ref_math_dot(&(eval[3]), &(eval[3])));
    dv = sqrt(ref_math_dot(&(eval[6]), &(eval[6])));
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, uv_range, &nchild,
                        &pchldrn, &psens),
         "EG topo face");
    /* move toward the center along the uv with largest dervative norm */
    if (du > dv) {
      params[0] =
          (1.0 - shift) * params[0] + shift * 0.5 * (uv_range[0] + uv_range[1]);
    } else {
      params[1] =
          (1.0 - shift) * params[1] + shift * 0.5 * (uv_range[2] + uv_range[3]);
    }
    egads_status = EG_curvature(object, params, curvature);
  }
  /* line or internal degen, move normal */
  if (EGADS_DEGEN == egads_status) {
    REF_DBL du, dv;
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double uv_range[4];
    double params[2];
    double eval[18];
    REF_DBL shift = 1.0e-2;
    params[0] = uv[0];
    params[1] = uv[1];
    REIS(EGADS_SUCCESS, EG_evaluate(object, params, eval), "eval derivs");
    du = sqrt(ref_math_dot(&(eval[3]), &(eval[3])));
    dv = sqrt(ref_math_dot(&(eval[6]), &(eval[6])));
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, uv_range, &nchild,
                        &pchldrn, &psens),
         "EG topo face");
    /* move toward the center along the uv with largest dervative norm */
    if (du <= dv) {
      params[0] =
          (1.0 - shift) * params[0] + shift * 0.5 * (uv_range[0] + uv_range[1]);
    } else {
      params[1] =
          (1.0 - shift) * params[1] + shift * 0.5 * (uv_range[2] + uv_range[3]);
    }
    egads_status = EG_curvature(object, params, curvature);
  }
  if (EGADS_SUCCESS == egads_status) {
    *kr = curvature[0];
    r[0] = curvature[1];
    r[1] = curvature[2];
    r[2] = curvature[3];
    *ks = curvature[4];
    s[0] = curvature[5];
    s[1] = curvature[6];
    s[2] = curvature[7];
    return REF_SUCCESS;
  }
  if (EGADS_NOTGEOM == egads_status) {
    /* EFFECTIVE */
    *kr = 0.0;
    r[0] = 1.0;
    r[1] = 0.0;
    r[2] = 0.0;
    *ks = 0.0;
    s[0] = 0.0;
    s[1] = 1.0;
    s[2] = 0.0;
    return REF_SUCCESS;
  }
  printf("EG_curvature %d (-24 is DEGEN) faceid %d u %f v %f\n", egads_status,
         faceid, uv[0], uv[1]);
  *kr = 0.0;
  r[0] = 1.0;
  r[1] = 0.0;
  r[2] = 0.0;
  *ks = 0.0;
  s[0] = 0.0;
  s[1] = 1.0;
  s[2] = 0.0;
  return REF_FAILURE;
#else
  printf("curvature 0, 0: No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(faceid);
  SUPRESS_UNUSED_COMPILER_WARNING(degen);
  SUPRESS_UNUSED_COMPILER_WARNING(uv);
  *kr = 0.0;
  r[0] = 1.0;
  r[1] = 0.0;
  r[2] = 0.0;
  *ks = 0.0;
  s[0] = 0.0;
  s[1] = 1.0;
  s[2] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_edge_trange(REF_GEOM ref_geom, REF_INT id,
                                 REF_DBL *trange) {
#ifdef HAVE_EGADS
  ego *edges;
  ego object;
  int periodic;
  int status;
  RNS(ref_geom->edges, "edges not loaded");
  edges = (ego *)(ref_geom->edges);
  RAS(1 <= id && id <= ref_geom->nedge, "edge id out of range");
  object = edges[id - 1];
  status = EG_getRange(object, trange, &periodic);
  /* returns -2 EGADS_NULLOBJ for EGADSlite of hemisphere
     REIB(EGADS_SUCCESS, EG_getRange(edge_ego, trange, &periodic),
     "edge trange", {
     printf("for edge %d (%p) face %d\n", edgeid, (void *)edge_ego,
     ref_geom_id(ref_geom, geom));
     });
  */
  /* use EG_getTopology as an alternate to EG_getRange */
  if (EGADS_NULLOBJ == status) {
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "topo for edge range");
  } else {
    REIS(EGADS_SUCCESS, status, "edge range");
  }
  return REF_SUCCESS;
#else
  printf("No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
  trange[0] = 0.0;
  trange[1] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_edge_crease(REF_GEOM ref_geom, REF_INT edgeid,
                                 REF_DBL *min_angle, REF_DBL *max_angle) {
#ifdef HAVE_EGADS
  REF_INT i, n = 11;
  REF_INT face, edge;
  REF_INT faces[2], nface;
  REF_DBL trange[2];

  ego esurf, *eloops;
  int oclass, mtype, nloop, *senses;
  double data[18];
  ego ecurve, *eedges;
  int iloop, iedge, nedge;

  RSS(ref_egads_edge_trange(ref_geom, edgeid, trange), "trange");

  nface = 0;
  faces[0] = REF_EMPTY;
  faces[1] = REF_EMPTY;
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
        edge = EG_indexBodyTopo((ego)(ref_geom->body), eedges[iedge]) - 1;
        if (edgeid == edge + 1) {
          RAB(2 > nface, "edge has more than 2 faces",
              printf("face ids %d %d %d edge id %d\n", faces[0], faces[1],
                     face + 1, edge + 1));
          faces[nface] = face + 1;
          nface++;
        }
      }
    }
  }
  REIS(2, nface, "expected two faces for edge");

  for (i = 0; i < n; i++) {
    REF_DBL ss, t;
    REF_DBL uv[2];
    REF_DBL kr, r[3], ks, s[3];
    REF_DBL n0[3], n1[3];
    REF_DBL angle;
    ss = i / (double)(n - 1);
    t = ss * trange[1] + (ss - 1.0) * trange[0];
    RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faces[0], 0, t, uv), "uv0");
    RSS(ref_egads_face_curvature_at(ref_geom, faces[0], 0, uv, &kr, r, &ks, s),
        "curve0");
    ref_math_cross_product(r, s, n0);
    RSS(ref_math_normalize(n0), "verify");
    RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faces[1], 0, t, uv), "uv0");
    RSS(ref_egads_face_curvature_at(ref_geom, faces[1], 0, uv, &kr, r, &ks, s),
        "curve0");
    ref_math_cross_product(r, s, n1);
    RSS(ref_math_normalize(n1), "verify");
    angle = ref_math_in_degrees(acos(ref_math_dot(n0, n1)));
    if (0 == i) {
      *min_angle = angle;
      *max_angle = angle;
    } else {
      *min_angle = MIN(*min_angle, angle);
      *max_angle = MIN(*max_angle, angle);
    }
  }
  return REF_SUCCESS;
#else
  printf("No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(edgeid);
  *min_angle = 0.0;
  *max_angle = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_edge_face_uv(REF_GEOM ref_geom, REF_INT edgeid,
                                  REF_INT faceid, REF_INT sense, REF_DBL t,
                                  REF_DBL *uv) {
#ifdef HAVE_EGADS
  ego *faces, *edges;
  ego face_ego, edge_ego;

  RNS(ref_geom->edges, "edges not loaded");
  edges = (ego *)(ref_geom->edges);
  RAS(1 <= edgeid && edgeid <= ref_geom->nedge, "edge id out of range");
  edge_ego = edges[edgeid - 1];

  RNS(ref_geom->faces, "faces not loaded");
  faces = (ego *)(ref_geom->faces);
  RAS(1 <= faceid && faceid <= ref_geom->nface, "face id out of range");
  face_ego = faces[faceid - 1];

  REIB(EGADS_SUCCESS, EG_getEdgeUV(face_ego, edge_ego, sense, t, uv),
       "eval edge face uv", {
         REF_DBL trange[2];
         printf("faceid %d edgeid %d sense %d t %.18e\n", faceid, edgeid, sense,
                t);
         printf("ref_egads_edge_trange status %d\n",
                ref_egads_edge_trange(ref_geom, edgeid, trange));
         printf("edgeid %d trange %.18e %.18e\n", edgeid, trange[0], trange[1]);
       });

  return REF_SUCCESS;
#else
  printf("No EGADS linked for %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(edgeid);
  SUPRESS_UNUSED_COMPILER_WARNING(faceid);
  SUPRESS_UNUSED_COMPILER_WARNING(sense);
  SUPRESS_UNUSED_COMPILER_WARNING(t);
  uv[0] = 0.0;
  uv[1] = 0.0;
  return REF_IMPLEMENT;
#endif
}

/*
  [x_t,y_t,z_t] edge [6]
  [x_tt,y_tt,z_tt]
  [x_u,y_u,z_u] [x_v,y_v,z_v] face [15]
  [x_uu,y_uu,z_uu] [x_uv,y_uv,z_uv] [x_vv,y_vv,z_vv]
*/
REF_STATUS ref_egads_eval(REF_GEOM ref_geom, REF_INT geom, REF_DBL *xyz,
                          REF_DBL *dxyz_dtuv) {
  REF_INT type, id, i;
  REF_DBL params[2];
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  params[0] = 0.0;
  params[1] = 0.0;
  type = ref_geom_type(ref_geom, geom);
  id = ref_geom_id(ref_geom, geom);

  for (i = 0; i < type; i++) {
    params[i] = ref_geom_param(ref_geom, i, geom);
  }
  RSS(ref_egads_eval_at(ref_geom, type, id, params, xyz, dxyz_dtuv), "eval at");
  return REF_SUCCESS;
}

REF_STATUS ref_egads_eval_at(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                             REF_DBL *params, REF_DBL *xyz,
                             REF_DBL *dxyz_dtuv) {
#ifdef HAVE_EGADS
  double eval[18];
  REF_INT i;
  ego *nodes, *edges, *faces;
  ego object;
  int status;

  if (NULL != ref_geom_facelift(ref_geom)) {
    RSS(ref_facelift_eval_at(ref_geom_facelift(ref_geom), type, id, params, xyz,
                             dxyz_dtuv),
        "facelift eval wrapper");
    return REF_SUCCESS;
  }

  object = (ego)NULL;
  switch (type) {
    case (REF_GEOM_NODE):
      RNS(ref_geom->nodes, "nodes not loaded");
      if (id < 1 || id > ref_geom->nnode) return REF_INVALID;
      nodes = (ego *)(ref_geom->nodes);
      object = nodes[id - 1];
      {
        ego ref, *pchldrn;
        int oclass, mtype, nchild, *psens;
        REIS(EGADS_SUCCESS,
             EG_getTopology(object, &ref, &oclass, &mtype, xyz, &nchild,
                            &pchldrn, &psens),
             "EG topo node");
      }
      return REF_SUCCESS;
    case (REF_GEOM_EDGE):
      RNS(ref_geom->edges, "edges not loaded");
      if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
      edges = (ego *)(ref_geom->edges);
      object = edges[id - 1];
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      if (id < 1 || id > ref_geom->nface) return REF_INVALID;
      faces = (ego *)(ref_geom->faces);
      object = faces[id - 1];
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom");
  }

  status = EG_evaluate(object, params, eval);
  if (EGADS_SUCCESS != status) {
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double range[4] = {-999, -999, -999, -999};
    printf(" %d EG_getTopology\n",
           EG_getTopology(object, &ref, &oclass, &mtype, range, &nchild,
                          &pchldrn, &psens));
    printf("range %f %f %f %f\n", range[0], range[1], range[2], range[3]);
    printf("type %d id %d\n", type, id);
    if (type > 0) printf("param[0] = %f\n", params[0]);
    if (type > 1) printf("param[1] = %f\n", params[1]);
    REIS(EGADS_SUCCESS, status, "eval");
  }
  xyz[0] = eval[0];
  xyz[1] = eval[1];
  xyz[2] = eval[2];
  if (NULL != dxyz_dtuv) {
    for (i = 0; i < 6; i++) dxyz_dtuv[i] = eval[3 + i];
    if (REF_GEOM_FACE == type)
      for (i = 0; i < 9; i++) dxyz_dtuv[6 + i] = eval[9 + i];
  }
  return REF_SUCCESS;
#else
  REF_INT i;
  printf("evaluating to (0,0,0), No EGADS linked for %s\n", __func__);
  printf("type %d id %d\n", type, id);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(params);
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  if (NULL != dxyz_dtuv) {
    for (i = 0; i < 6; i++) {
      dxyz_dtuv[i] = 0.0;
    }
    if (REF_GEOM_FACE == type) {
      for (i = 0; i < 9; i++) dxyz_dtuv[6 + i] = 0.0;
    }
  }
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_inverse_eval(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                                  REF_DBL *xyz, REF_DBL *param) {
#ifdef HAVE_EGADS
  ego object = (ego)NULL;
  int i, guess_status, noguess_status;
  REF_DBL guess_param[2], noguess_param[2];
  REF_DBL guess_closest[3], noguess_closest[3];
  REF_DBL guess_dist, noguess_dist;

  ego ref, *pchldrn;
  int oclass, mtype, nchild, *psens;
  double range[4];

  REF_BOOL guess_in_range, noguess_in_range;

  REF_BOOL verbose = REF_FALSE;

  if (NULL != ref_geom_facelift(ref_geom)) {
    RSS(ref_facelift_inverse_eval(ref_geom_facelift(ref_geom), type, id, xyz,
                                  param),
        "facelift inv eval wrapper");
    return REF_SUCCESS;
  }

  switch (type) {
    case (REF_GEOM_NODE):
      printf("GEOM_NODE ref_geom_inverse_eval not defined\n");
      return REF_IMPLEMENT;
    case (REF_GEOM_EDGE):
      RNS(ref_geom->edges, "edges not loaded");
      if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
      object = ((ego *)(ref_geom->edges))[id - 1];
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      if (id < 1 || id > ref_geom->nface) return REF_INVALID;
      object = ((ego *)(ref_geom->faces))[id - 1];
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom");
  }

  REIS(EGADS_SUCCESS,
       EG_getTopology(object, &ref, &oclass, &mtype, range, &nchild, &pchldrn,
                      &psens),
       "EG topo node");

  for (i = 0; i < type; i++) guess_param[i] = param[i];
  for (i = 0; i < type; i++) noguess_param[i] = param[i];
  guess_status = EG_invEvaluateGuess(object, xyz, guess_param, guess_closest);
  noguess_status = EG_invEvaluate(object, xyz, noguess_param, noguess_closest);
  if (verbose) {
    printf("guess %d noguess %d type %d id %d xyz %f %f %f\n", guess_status,
           noguess_status, type, id, xyz[0], xyz[1], xyz[2]);
    for (i = 0; i < type; i++) {
      printf("%d: start %f guess %f noguess %f\n", i, param[i], guess_param[i],
             noguess_param[i]);
    }
  }
  if (EGADS_SUCCESS != guess_status && EGADS_SUCCESS != noguess_status)
    return REF_FAILURE;

  if (EGADS_SUCCESS == guess_status && EGADS_SUCCESS != noguess_status) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
    return REF_SUCCESS;
  }

  if (EGADS_SUCCESS != guess_status && EGADS_SUCCESS == noguess_status) {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
    return REF_SUCCESS;
  }

  guess_dist = sqrt(pow(guess_closest[0] - xyz[0], 2) +
                    pow(guess_closest[1] - xyz[1], 2) +
                    pow(guess_closest[2] - xyz[2], 2));
  noguess_dist = sqrt(pow(noguess_closest[0] - xyz[0], 2) +
                      pow(noguess_closest[1] - xyz[1], 2) +
                      pow(noguess_closest[2] - xyz[2], 2));

  guess_in_range = REF_TRUE;
  for (i = 0; i < type; i++)
    if (guess_param[i] < range[0 + 2 * i] || range[1 + 2 * i] < guess_param[i])
      guess_in_range = REF_FALSE;

  noguess_in_range = REF_TRUE;
  for (i = 0; i < type; i++)
    if (noguess_param[i] < range[0 + 2 * i] ||
        range[1 + 2 * i] < noguess_param[i])
      noguess_in_range = REF_FALSE;

  if (verbose) {
    printf("guess %e noguess %e\n", guess_dist, noguess_dist);
    printf("guess %d noguess %d in range\n", guess_in_range, noguess_in_range);
  }

  if (guess_in_range && !noguess_in_range) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
    return REF_SUCCESS;
  }

  if (!guess_in_range && noguess_in_range) {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
    return REF_SUCCESS;
  }

  if (guess_dist < noguess_dist) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
  } else {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
  }
  return REF_SUCCESS;

#else
  printf("no-op, No EGADS linked for %s type %d id %d x %f p %f n %d\n",
         __func__, type, id, xyz[0], param[0], ref_geom_n(ref_geom));
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_egads_invert(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                            REF_DBL *xyz0, REF_DBL *uv) {
  REF_DBL xyz[3], error[3];
  REF_DBL dxyz_dtuv[15];
  REF_DBL drsduv[4], duvdrs[4];
  REF_DBL r[3], s[3], n[3];
  REF_DBL drs[2], duv[2];
  REF_INT i;
  REF_DBL tol = 1.0e-14;
  REF_BOOL verbose = REF_FALSE;

  for (i = 0; i < 20; i++) {
    RSS(ref_egads_eval_at(ref_geom, type, id, uv, xyz, dxyz_dtuv), "eval");
    error[0] = xyz[0] - xyz0[0];
    error[1] = xyz[1] - xyz0[1];
    error[2] = xyz[2] - xyz0[2];
    RSS(ref_geom_face_rsn(ref_geom, id, uv, r, s, n), "rsn");
    /* [x_u,y_u,z_u] [x_v,y_v,z_v] */
    /*  drsduv = [r s] * dxyz_dtuv */
    drsduv[0] = dxyz_dtuv[0] * r[0] + dxyz_dtuv[1] * r[1] + dxyz_dtuv[2] * r[2];
    drsduv[1] = dxyz_dtuv[0] * s[0] + dxyz_dtuv[1] * s[1] + dxyz_dtuv[2] * s[2];
    drsduv[2] = dxyz_dtuv[3] * r[0] + dxyz_dtuv[4] * r[1] + dxyz_dtuv[5] * r[2];
    drsduv[3] = dxyz_dtuv[3] * s[0] + dxyz_dtuv[4] * s[1] + dxyz_dtuv[5] * s[2];
    RSS(ref_matrix_inv_gen(2, drsduv, duvdrs), "inv");
    drs[0] = ref_math_dot(r, error);
    drs[1] = ref_math_dot(s, error);
    if (verbose)
      printf(" r %e s %e err %e on %d\n", drs[0], drs[1],
             sqrt(ref_math_dot(error, error)), i);
    duv[0] = duvdrs[0] * drs[0] + duvdrs[2] * drs[1];
    duv[1] = duvdrs[1] * drs[0] + duvdrs[3] * drs[1];
    uv[0] = uv[0] - duv[0];
    uv[1] = uv[1] - duv[1];
    if (sqrt(drs[0] * drs[0] + drs[1] * drs[1]) < tol) {
      return REF_SUCCESS;
    }
  }
  return REF_FAILURE;
}

REF_STATUS ref_egads_gap(REF_GEOM ref_geom, REF_INT node, REF_DBL *gap) {
  REF_INT item, geom, type;
  REF_DBL dist, face_xyz[3], gap_xyz[3];
  REF_BOOL has_node, has_edge;
  *gap = 0.0;

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &has_node), "n");
  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &has_edge), "n");
  if (!has_edge) return REF_SUCCESS;

  if (has_node) {
    type = REF_GEOM_NODE;
  } else {
    type = REF_GEOM_FACE;
  }
  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      RSS(ref_egads_eval(ref_geom, geom, gap_xyz, NULL), "eval");
    }
  }

  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
      RSS(ref_egads_eval(ref_geom, geom, face_xyz, NULL), "eval");
      dist = sqrt(pow(face_xyz[0] - gap_xyz[0], 2) +
                  pow(face_xyz[1] - gap_xyz[1], 2) +
                  pow(face_xyz[2] - gap_xyz[2], 2));
      (*gap) = MAX((*gap), dist);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_egads_feature_size(REF_GRID ref_grid, REF_INT node, REF_DBL *h0,
                                  REF_DBL *dir0, REF_DBL *h1, REF_DBL *dir1,
                                  REF_DBL *h2, REF_DBL *dir2) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT edge_item, face_item, edge_geom, face_geom, edgeid, faceid, iloop;
  ego ref, *cadnodes, *edges, *loops;
  int oclass, mtype, ncadnode, nedge, nloop, *senses;
  double data[18];
  double trange[2];
  REF_DBL xyz[3];
  REF_DBL dxyz_dt[6];
  REF_DBL xyz1[3];
  REF_INT ineligible_cad_node0, ineligible_cad_node1;
  REF_INT cad_node0, cad_node1;
  REF_DBL param[2];
  REF_INT other_edgeid, iedge;
  REF_DBL len, diagonal, hmax;
  REF_DBL tangent[3], dx[3], dot, orth[3];
  REF_DBL gap;
  int status;

  /* initialize isotropic with bounding box */
  RSS(ref_egads_diagonal(ref_geom, REF_EMPTY, &diagonal), "bbox diag init");
  hmax = diagonal;
  hmax /= MAX(1.0, ref_geom_segments_per_bounding_box_diagonal(ref_geom));
  *h0 = hmax;
  dir0[0] = 1.0;
  dir0[1] = 0.0;
  dir0[2] = 0.0;
  *h1 = hmax;
  dir1[0] = 0.0;
  dir1[1] = 1.0;
  dir1[2] = 0.0;
  *h2 = hmax;
  dir2[0] = 0.0;
  dir2[1] = 0.0;
  dir2[2] = 1.0;

  each_ref_geom_having_node(ref_geom, node, edge_item, edge_geom) {
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, edge_geom)) {
      edgeid = ref_geom_id(ref_geom, edge_geom);
      REIS(EGADS_SUCCESS,
           EG_getTopology(((ego *)(ref_geom->edges))[edgeid - 1], &ref, &oclass,
                          &mtype, trange, &ncadnode, &cadnodes, &senses),
           "EG topo edge");
      RAS(mtype != DEGENERATE, "edge interior DEGENERATE");
      RSS(ref_egads_eval(ref_geom, edge_geom, xyz, dxyz_dt), "eval edge");
      tangent[0] = dxyz_dt[0];
      tangent[1] = dxyz_dt[1];
      tangent[2] = dxyz_dt[2];
      if (REF_SUCCESS != ref_math_normalize(tangent)) {
        tangent[0] = 0.0;
        tangent[1] = 0.0;
        tangent[2] = 0.0;
      }
      RAS(0 < ncadnode && ncadnode < 3, "edge children");
      ineligible_cad_node0 = EG_indexBodyTopo(ref_geom->body, cadnodes[0]);
      if (2 == ncadnode) {
        ineligible_cad_node1 = EG_indexBodyTopo(ref_geom->body, cadnodes[1]);
      } else {
        ineligible_cad_node1 = ineligible_cad_node0; /* ONENODE edge */
      }
      each_ref_geom_having_node(ref_geom, node, face_item, face_geom) {
        if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
          faceid = ref_geom_id(ref_geom, face_geom);
          REIS(EGADS_SUCCESS,
               EG_getTopology(((ego *)(ref_geom->faces))[faceid - 1], &ref,
                              &oclass, &mtype, data, &nloop, &loops, &senses),
               "topo");
          for (iloop = 0; iloop < nloop; iloop++) {
            /* loop through all Edges associated with this Loop */
            REIS(EGADS_SUCCESS,
                 EG_getTopology(loops[iloop], &ref, &oclass, &mtype, data,
                                &nedge, &edges, &senses),
                 "topo");
            for (iedge = 0; iedge < nedge; iedge++) {
              other_edgeid =
                  EG_indexBodyTopo((ego)(ref_geom->body), edges[iedge]);
              /* qualified? does not share geom nodes */
              REIS(EGADS_SUCCESS,
                   EG_getTopology(((ego *)(ref_geom->edges))[other_edgeid - 1],
                                  &ref, &oclass, &mtype, trange, &ncadnode,
                                  &cadnodes, &senses),
                   "EG topo node");
              if (mtype == DEGENERATE) continue; /* skip DEGENERATE */
              RAS(0 < ncadnode && ncadnode < 3, "edge children");
              cad_node0 = EG_indexBodyTopo(ref_geom->body, cadnodes[0]);
              if (2 == ncadnode) {
                cad_node1 = EG_indexBodyTopo(ref_geom->body, cadnodes[1]);
              } else {
                cad_node1 = cad_node0; /* ONENODE edge */
              }
              if (cad_node0 == ineligible_cad_node0 ||
                  cad_node0 == ineligible_cad_node1 ||
                  cad_node1 == ineligible_cad_node0 ||
                  cad_node1 == ineligible_cad_node1)
                continue;
              /* inverse projection */
              status =
                  EG_invEvaluate(((ego *)(ref_geom->edges))[other_edgeid - 1],
                                 xyz, param, xyz1);
              REIS(EGADS_SUCCESS, status, "EG_invEvaluate opp edge");
              dx[0] = xyz1[0] - xyz[0];
              dx[1] = xyz1[1] - xyz[1];
              dx[2] = xyz1[2] - xyz[2];
              len = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
              RSS(ref_egads_gap(ref_geom, node, &gap), "edge gap");
              len = MAX(len, gap);
              RSS(ref_math_normalize(dx), "direction across face");
              if (len < *h0) {
                dot = ref_math_dot(tangent, dx);
                orth[0] = tangent[0] - dot * dx[0];
                orth[1] = tangent[1] - dot * dx[1];
                orth[2] = tangent[2] - dot * dx[2];
                if (REF_SUCCESS == ref_math_normalize(orth)) {
                  *h0 = len;
                  dir0[0] = dx[0];
                  dir0[1] = dx[1];
                  dir0[2] = dx[2];
                  RSS(ref_egads_diagonal(ref_geom, edge_geom, h1),
                      "local diag");
                  dir1[0] = orth[0];
                  dir1[1] = orth[1];
                  dir1[2] = orth[2];
                  *h2 = diagonal;
                  ref_math_cross_product(dir0, dir1, dir2);
                } else {
                  *h0 = len;
                  dir0[0] = dx[0];
                  dir0[1] = dx[1];
                  dir0[2] = dx[2];
                  *h1 = diagonal;
                  *h2 = diagonal;
                  RSS(ref_math_orthonormal_system(dir0, dir1, dir2),
                      "arbitrary orthonormal dir1 dir2");
                }
              }
            }
          }
        }
      }
    }
  }
#else
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  RSS(ref_egads_diagonal(ref_geom, REF_EMPTY, h0), "bbox diag init");
  dir0[0] = 1.0;
  dir0[1] = 0.0;
  dir0[2] = 0.0;
  *h1 = *h0;
  *h2 = *h0;
  RSS(ref_math_orthonormal_system(dir0, dir1, dir2),
      "arbitrary orthonormal dir1 dir2");
  SUPRESS_UNUSED_COMPILER_WARNING(node);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_egads_geom_cost(REF_GEOM ref_geom, REF_INT type, REF_INT id) {
#ifdef HAVE_EGADS
  ego object;
  ego geom, *bodies, ref;
  int oclass, mtype, nbody, *senses;
  int *pinfo = NULL;
  double *prv = NULL;
  switch (type) {
    case REF_GEOM_EDGE:
      object = ((ego *)(ref_geom->edges))[id - 1];
      break;
    case REF_GEOM_FACE:
      object = ((ego *)(ref_geom->faces))[id - 1];
      break;
    default:
      return REF_SUCCESS;
  }
  REIS(EGADS_SUCCESS,
       EG_getTopology(object, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                      &senses),
       "EG getTopology");
  REIS(EGADS_SUCCESS, EG_getGeometry(geom, &oclass, &mtype, &ref, &pinfo, &prv),
       "EG getGeometry");
  if (mtype == BEZIER || mtype == BSPLINE)
    printf("%d %d pinfo %d %d %d %d %d %d %d\n", type, id, pinfo[0], pinfo[1],
           pinfo[2], pinfo[3], pinfo[4], pinfo[5], pinfo[6]);
  EG_free(pinfo);
  EG_free(prv);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
#endif
  return REF_SUCCESS;
}

#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE) && \
    defined(HAVE_EGADS_EFFECTIVE)
static REF_STATUS ref_egads_quilt_attributes(ego body, ego ebody) {
  int nface, i;
  ego *faces;
  REF_DICT ref_dict;
  REF_INT key, flag, minflag, maxflag, n, value;
  RSS(ref_dict_create(&ref_dict), "create");
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, FACE, &nface, &faces),
       "EG face topo");
  for (i = 0; i < nface; i++) {
    int len, atype;
    const double *preals;
    const int *pints;
    const char *string;

    if (EGADS_SUCCESS == EG_attributeRet(faces[i], "group", &atype, &len,
                                         &pints, &preals, &string)) {
      if (ATTRINT == atype) RSS(ref_dict_store(ref_dict, i, pints[0]), "store");
      if (ATTRREAL == atype) {
        double dbl = round(preals[0]);
        RSS(ref_dict_store(ref_dict, i, (REF_INT)dbl), "store");
      }
    }
  }
  minflag = REF_INT_MAX;
  maxflag = REF_INT_MIN;
  each_ref_dict_key_value(ref_dict, i, key, flag) {
    minflag = MIN(flag, minflag);
    maxflag = MAX(flag, maxflag);
  }
  if (minflag <= maxflag) printf("flag range %d %d\n", minflag, maxflag);
  for (flag = minflag; flag <= maxflag; flag++) {
    n = 0;
    each_ref_dict_key_value(ref_dict, i, key, value) {
      if (flag == value) {
        n++;
      }
    }
    if (n > 0) {
      ego eface, *group_faces;
      ref_malloc(group_faces, n, ego);
      printf("group %d has %d members", flag, n);
      n = 0;
      each_ref_dict_key_value(ref_dict, i, key, value) {
        if (flag == value) {
          printf(" %d", key + 1);
          group_faces[n] = faces[key];
          n++;
        }
      }
      printf("\n");
      if (n > 1)
        REIS(EGADS_SUCCESS, EG_makeEFace(ebody, n, group_faces, &eface),
             "initEB");
      ref_free(group_faces);
    }
  }
  EG_free(faces);
  RSS(ref_dict_free(ref_dict), "free");
  return REF_SUCCESS;
}
#endif

#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE) && \
    defined(HAVE_EGADS_EFFECTIVE)
static REF_STATUS ref_egads_quilt_angle(REF_GEOM ref_geom, ego body, ego ebody,
                                        double angle, REF_INT *e2f) {
  int nedge, edge, nface, i;
  ego *faces;
  REF_DICT ref_dict;
  REF_INT key, flag, minflag, maxflag, n, value, relaxations;
  REF_BOOL rerun;
  RSS(ref_dict_create(&ref_dict), "create");
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL),
       "EG edge topo");
  REIS(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, FACE, &nface, &faces),
       "EG face topo");
  for (i = 0; i < nface; i++) {
    RSS(ref_dict_store(ref_dict, i, i), "store");
  }
  rerun = REF_TRUE;
  relaxations = 0;
  while (rerun) {
    rerun = REF_FALSE;
    for (edge = 0; edge < nedge; edge++) {
      REF_INT group0, group1;
      REF_DBL min_angle, max_angle;
      if (REF_EMPTY != e2f[0 + 2 * edge] && REF_EMPTY != e2f[1 + 2 * edge]) {
        RSS(ref_dict_value(ref_dict, e2f[0 + 2 * edge] - 1, &group0), "g0");
        RSS(ref_dict_value(ref_dict, e2f[1 + 2 * edge] - 1, &group1), "g1");
        if (group0 != group1) {
          RSS(ref_egads_edge_crease(ref_geom, edge + 1, &min_angle, &max_angle),
              "crease");
          if (max_angle < angle) {
            RSS(ref_dict_store(ref_dict, e2f[0 + 2 * edge] - 1,
                               MIN(group0, group1)),
                "store0");
            RSS(ref_dict_store(ref_dict, e2f[1 + 2 * edge] - 1,
                               MIN(group0, group1)),
                "store1");
            rerun = REF_TRUE;
          }
        }
      }
    }
    relaxations++;
    if (relaxations % 100 == 0) printf("relaxation %d\n", relaxations);
  }

  minflag = REF_INT_MAX;
  maxflag = REF_INT_MIN;
  each_ref_dict_key_value(ref_dict, i, key, flag) {
    minflag = MIN(flag, minflag);
    maxflag = MAX(flag, maxflag);
  }
  if (minflag <= maxflag) printf("flag range %d %d\n", minflag, maxflag);
  for (flag = minflag; flag <= maxflag; flag++) {
    n = 0;
    each_ref_dict_key_value(ref_dict, i, key, value) {
      if (flag == value) {
        n++;
      }
    }
    if (n > 0) {
      ego eface, *group_faces;
      ref_malloc(group_faces, n, ego);
      printf("flag %d has %d members", flag, n);
      n = 0;
      each_ref_dict_key_value(ref_dict, i, key, value) {
        if (flag == value) {
          printf(" %d", key + 1);
          group_faces[n] = faces[key];
          n++;
        }
      }
      printf("\n");
      if (n > 1)
        REIS(EGADS_SUCCESS, EG_makeEFace(ebody, n, group_faces, &eface),
             "initEB");
      ref_free(group_faces);
    }
  }
  EG_free(faces);
  RSS(ref_dict_free(ref_dict), "free");
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_egads_quilt(REF_GEOM ref_geom, REF_INT auto_tparams,
                           REF_DBL *global_params) {
#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE) && \
    defined(HAVE_EGADS_EFFECTIVE)
  ego effective[2];
  ego tess, model;
  double angle;
  REF_BOOL quilt_on_angle = REF_FALSE;

  RAS(ref_geom_model_loaded(ref_geom), "load model before quilting");
  RAS(!ref_geom_effective(ref_geom), "already effective, quilting twice?");

  REIS(EGADS_SUCCESS,
       EG_copyObject((ego)(ref_geom->body), NULL, &(effective[0])),
       "copy body");

  /* need to use copy to build tess so they match */
  ref_geom->body = effective[0];
  if (NULL != ref_geom->faces) EG_free((ego *)(ref_geom->faces));
  if (NULL != ref_geom->edges) EG_free((ego *)(ref_geom->edges));
  if (NULL != ref_geom->nodes) EG_free((ego *)(ref_geom->nodes));
  ref_free(ref_geom->face_seg_per_rad);
  ref_free(ref_geom->face_min_length);
  ref_free(ref_geom->initial_cell_height);
  ref_free(ref_geom->uv_area_sign);
  RSS(ref_egads_cache_body_objects(ref_geom), "cache egads objects");
  RSS(ref_egads_tess_create(ref_geom, &tess, auto_tparams, global_params),
      "create tess object");

  angle = 10.0;
  REIS(EGADS_SUCCESS, EG_initEBody(tess, angle, &effective[1]), "init xEB");

  RSS(ref_egads_quilt_attributes(effective[0], effective[1]), "quilt attr");

  if (quilt_on_angle) {
    REF_INT *e2f;
    RSS(ref_egads_edge_faces(ref_geom, &e2f), "edge2face");
    RSS(ref_egads_quilt_angle(ref_geom, effective[0], effective[1], angle, e2f),
        "quilt attr");
    ref_free(e2f);
  }

  REIS(EGADS_SUCCESS, EG_finishEBody(effective[1]), "complete EB");
  REIS(0, EG_deleteObject(tess), "delete tess copied into ebody");

  REIS(0, EG_deleteObject((ego)(ref_geom->model)), "delete body model");
  REIS(EGADS_SUCCESS,
       EG_makeTopology((ego)(ref_geom->context), NULL, MODEL, 2, NULL, 1,
                       effective, NULL, &model),
       "make Topo Model");
  ref_geom->model = (void *)model;

  ref_geom->body = (void *)effective[1];
  ref_geom_effective(ref_geom) = REF_TRUE;

  if (NULL != ref_geom->faces) EG_free((ego *)(ref_geom->faces));
  if (NULL != ref_geom->edges) EG_free((ego *)(ref_geom->edges));
  if (NULL != ref_geom->nodes) EG_free((ego *)(ref_geom->nodes));
  ref_free(ref_geom->face_seg_per_rad);
  ref_free(ref_geom->face_min_length);
  ref_free(ref_geom->initial_cell_height);
  ref_free(ref_geom->uv_area_sign);
  RSS(ref_egads_cache_body_objects(ref_geom), "cache egads objects");

  return REF_SUCCESS;
#else
  printf("no-op, EGADS not linked with HAVE_EGADS_EFFECTIVE %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(auto_tparams);
  SUPRESS_UNUSED_COMPILER_WARNING(global_params);
  return REF_SUCCESS;
#endif
}

REF_STATUS ref_egads_add_attribute(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                                   const char *name, const char *value) {
#ifdef HAVE_EGADS
#ifdef HAVE_EGADS_LITE
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
  SUPRESS_UNUSED_COMPILER_WARNING(name);
  SUPRESS_UNUSED_COMPILER_WARNING(value);
  RSS(REF_IMPLEMENT, "full EGADS required to add attribute");
#else
  {
    ego object = NULL;
    int ignored_len = -1;
    switch (type) {
      case (REF_GEOM_NODE):
        RNS(ref_geom->nodes, "nodes not loaded");
        if (id < 1 || id > ref_geom->nnode) return REF_INVALID;
        object = ((ego *)(ref_geom->nodes))[id - 1];
        break;
      case (REF_GEOM_EDGE):
        RNS(ref_geom->edges, "edges not loaded");
        if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
        object = ((ego *)(ref_geom->edges))[id - 1];
        break;
      case (REF_GEOM_FACE):
        RNS(ref_geom->faces, "faces not loaded");
        if (id < 1 || id > ref_geom->nface) return REF_INVALID;
        object = ((ego *)(ref_geom->faces))[id - 1];
        break;
      case (REF_GEOM_SOLID):
        RNS(ref_geom->body, "body not loaded");
        object = (ego)(ref_geom->body);
        break;
      default:
        RSS(REF_FAILURE, "unknown type");
    }

    REIS(EGADS_SUCCESS,
         EG_attributeAdd(object, name, ATTRSTRING, ignored_len, NULL, NULL,
                         value),
         "add attribute");
  }
#endif
  return REF_SUCCESS;
#else
  printf("no-op, EGADS not linked with HAVE_EGADS_EFFECTIVE %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
  SUPRESS_UNUSED_COMPILER_WARNING(name);
  SUPRESS_UNUSED_COMPILER_WARNING(value);
  return REF_SUCCESS;
#endif
}

REF_STATUS ref_egads_get_attribute(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                                   const char *name, const char **value) {
#ifdef HAVE_EGADS
  ego object = NULL;
  int attribute_type, len;
  const int *ints;
  const double *reals;
  int egads_status;

  *value = NULL;

  switch (type) {
    case (REF_GEOM_NODE):
      RNS(ref_geom->nodes, "nodes not loaded");
      if (id < 1 || id > ref_geom->nnode) return REF_INVALID;
      object = ((ego *)(ref_geom->nodes))[id - 1];
      break;
    case (REF_GEOM_EDGE):
      RNS(ref_geom->edges, "edges not loaded");
      if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
      object = ((ego *)(ref_geom->edges))[id - 1];
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      if (id < 1 || id > ref_geom->nface) return REF_INVALID;
      object = ((ego *)(ref_geom->faces))[id - 1];
      break;
    case (REF_GEOM_SOLID):
      RNS(ref_geom->body, "body not loaded");
      object = (ego)(ref_geom->body);
      break;
    default:
      RSS(REF_FAILURE, "unknown type");
  }

  egads_status = EG_attributeRet(object, name, &attribute_type, &len, &ints,
                                 &reals, value);
  if (EGADS_NOTFOUND == egads_status) return REF_NOT_FOUND;
  REIS(EGADS_SUCCESS, egads_status, "get/return attribute");
  REIS(ATTRSTRING, attribute_type, "expected string");

  return REF_SUCCESS;
#else
  *value = NULL;
  printf("no-op, EGADS not linked with HAVE_EGADS_EFFECTIVE %s\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(id);
  SUPRESS_UNUSED_COMPILER_WARNING(name);
  return REF_SUCCESS;
#endif
}

REF_STATUS ref_egads_extract_mapbc(REF_GEOM ref_geom, const char *mapbc) {
  FILE *file;
  file = fopen(mapbc, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", mapbc);
  RNS(file, "unable to open file");

  if (ref_geom->manifold) {
    REF_INT face_id;
    fprintf(file, "%d\n", ref_geom->nface);
    for (face_id = 1; face_id <= ref_geom->nface; face_id++) {
      const char *attribute = NULL;
      char *bc_name;
      REF_SIZE len, i;
      RSS(ref_egads_get_attribute(ref_geom, REF_GEOM_FACE, face_id, "bc_name",
                                  &attribute),
          "get");
      RNS(attribute, "attribute NULL");
      len = strlen(attribute);
      RAS(10000 > len, "attribute more than 10000 bytes");
      ref_malloc(bc_name, (REF_LONG)(len + 1), char);
      strcpy(bc_name, attribute);
      for (i = 0; i < len; i++) {
        if ('_' == bc_name[i]) {
          bc_name[i] = ' ';
          break;
        }
      }
      fprintf(file, "%d %s\n", face_id, bc_name);
      ref_free(bc_name);
    }
  } else {
    REF_INT edge_id;
    fprintf(file, "%d\n", 2 + ref_geom->nedge);
    for (edge_id = 1; edge_id <= ref_geom->nedge; edge_id++) {
      const char *attribute = NULL;
      char *bc_name;
      REF_SIZE len, i;
      RSS(ref_egads_get_attribute(ref_geom, REF_GEOM_EDGE, edge_id, "bc_name",
                                  &attribute),
          "get");
      RNS(attribute, "attribute NULL");
      len = strlen(attribute);
      RAS(10000 > len, "attribute more than 10000 bytes");
      ref_malloc(bc_name, (REF_LONG)(len + 1), char);
      strcpy(bc_name, attribute);
      for (i = 0; i < len; i++) {
        if ('_' == bc_name[i]) {
          bc_name[i] = ' ';
          break;
        }
      }
      fprintf(file, "%d %s\n", edge_id, bc_name);
      ref_free(bc_name);
    }
    edge_id = ref_geom->nedge + 1;
    fprintf(file, "%d %s\n", edge_id, "6662 symmetry-y-min");
    edge_id = ref_geom->nedge + 2;
    fprintf(file, "%d %s\n", edge_id, "6662 symmetry-y-max");
  }
  fclose(file);

  return REF_SUCCESS;
}
