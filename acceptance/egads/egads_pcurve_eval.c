/*

gcc-7  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized
-I/Users/mpark/esp/EngSketchPad/include -o egads_pcurve_eval egads_pcurve_eval.c
-Wl,-rpath,/Users/mpark/esp/EngSketchPad/lib -L/Users/mpark/esp/EngSketchPad/lib
-legads   -lm

*/

#include <stdio.h>
#include <stdlib.h>

#include "egads.h"

#define is_equal(a, b, msg)                                                  \
  {                                                                          \
    int ref_private_status_ai, ref_private_status_bi;                        \
    ref_private_status_ai = (a);                                             \
    ref_private_status_bi = (b);                                             \
    if (ref_private_status_ai != ref_private_status_bi) {                    \
      printf("%s: %d: %s: %s\nexpected %d was %d\n", __FILE__, __LINE__,     \
             __func__, (msg), ref_private_status_ai, ref_private_status_bi); \
      return 1;                                                              \
    }                                                                        \
  }

#define is_true(a, msg)                                                \
  {                                                                    \
    if (!(a)) {                                                        \
      printf("%s: %d: %s: %s\n", __FILE__, __LINE__, __func__, (msg)); \
      return 1;                                                        \
    }                                                                  \
  }

int main(void) {
  ego context;
  ego model = NULL;
  ego geom, *bodies, *children;
  int oclass, mtype, nbody, *senses, nchild;
  ego solid;
  int nface;
  ego *faces;
  int faceid;
  int nedge;
  ego *edges;
  int edgeid;
  ego esurf, *eloops;
  int nloop;
  double data[18];
  int sense;
  double t, param[2];
  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 3) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, "hemisphere.egads", &model),
           "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
  is_equal(1, nbody, "expected 1 body");
  solid = bodies[0];
  is_equal(EGADS_SUCCESS,
           EG_getTopology(solid, &geom, &oclass, &mtype, NULL, &nchild,
                          &children, &senses),
           "EG topo body type");
  is_equal(SOLIDBODY, mtype, "expected SOLIDBODY");

  is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, FACE, &nface, &faces),
           "EG face topo");

  faceid = 1;
  is_equal(EGADS_SUCCESS,
           EG_getTopology(faces[faceid - 1], &esurf, &oclass, &mtype, data,
                          &nloop, &eloops, &senses),
           "topo");
  printf("faceid %d uv box ([%f,%f],[%f,%f])\n", faceid, data[0], data[1],
         data[2], data[3]);

  is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EDGE, &nedge, &edges),
           "EG face topo");

  edgeid = 3;
  sense = 0; /* when edge is used twice in loop, this is its index */
  t = 0.5;
  is_equal(EGADS_SUCCESS,
           EG_getEdgeUV(faces[faceid - 1], edges[edgeid - 1], sense, t, param),
           "eval edge face uv");

  printf("edgeid %d t %f pcurve eval on faceid %d is (%f,%f)\n", edgeid, t,
         faceid, param[0], param[1]);

  is_equal(EGADS_SUCCESS,
           EG_getTopology(edges[edgeid - 1], &esurf, &oclass, &mtype, data,
                          &nloop, &eloops, &senses),
           "topo");
  printf("edgeid %d t box [%f,%f]\n", edgeid, data[0], data[1]);

  return 0;
}
