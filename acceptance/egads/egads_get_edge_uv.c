/*

gcc-10  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I/Users/mpark/local/pkgs/EngSketchPad/include -o egads_get_edge_uv-beta \
egads_get_edge_uv.c -Wl,-rpath,/Users/mpark/local/pkgs/EngSketchPad/lib \
-L/Users/mpark/local/pkgs/EngSketchPad/lib -legads   -lm

gcc-10  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I/Users/mpark/local/pkgs/EGADS/trunk/include -o egads_get_edge_uv-svn \
egads_get_edge_uv.c -Wl,-rpath,/Users/mpark/local/pkgs/EGADS/trunk/lib \
-L/Users/mpark/local/pkgs/EGADS/trunk/lib -legads   -lm

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
  double t;
  double uv[2];
  ego face_ego, edge_ego;
  int sense;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 3) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, "c40f.egads", &model),
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
  is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EDGE, &nedge, &edges),
           "EG edge topo");

  printf("nface %d nedge %d\n", nface, nedge);

  faceid = 168;
  edgeid = 723;

  edge_ego = edges[edgeid - 1];
  face_ego = faces[faceid - 1];
  t = 287.984998;
  sense = 0;

  is_equal(EGADS_SUCCESS, EG_getEdgeUV(face_ego, edge_ego, sense, t, uv),
           "eval edge face uv");

  return 0;
}
