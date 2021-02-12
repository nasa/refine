/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/EGADS/trunk/include -o egads_get_edge_uv \
egads_get_edge_uv.c -Wl,-rpath,${HOME}/local/pkgs/EGADS/trunk/lib \
-L${HOME}/local/pkgs/EGADS/trunk/lib -legads   -lm

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
  ego body = NULL;
  ego geom, *bodies;
  int oclass, mtype, nbody, *senses, nchild, nego;
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

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, "c40f-eff.egads", &body),
           "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(body, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
  is_equal(1, nbody, "expected 1 body");
  body = bodies[0];

  {
    int ibody;
    int bodyclass, bodytype;
    ego owner, prev, next;
    for (ibody = 0; ibody < nego; ibody++) {
      is_equal(EGADS_SUCCESS,
               EG_getInfo(bodies[ibody], &bodyclass, &bodytype, &owner, &prev,
                          &next),
               "info");
      if (EBODY == bodyclass) {
        body = bodies[ibody];
        printf("ego %d is an EBODY\n", ibody);
      }
    }
  }

  is_equal(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EFACE, &nface, &faces),
           "EG face topo");
  is_equal(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EEDGE, &nedge, &edges),
           "EG edge topo");

  printf("nface %d nedge %d\n", nface, nedge);

  faceid = 13;
  edgeid = 45;

  edge_ego = edges[edgeid - 1];
  face_ego = faces[faceid - 1];
  sense = 0;

  {
    ego ref;
    double trange[2];
    ego *pchldrn;
    int *psens;
    is_equal(EGADS_SUCCESS,
             EG_getTopology(edge_ego, &ref, &oclass, &mtype, trange, &nchild,
                            &pchldrn, &psens),
             "EG topo edge");
    printf("edge %d mtype %d (5 is DEGENERATE)\n", edgeid, mtype);
    t = trange[0];
    is_equal(EGADS_SUCCESS, EG_getEdgeUV(face_ego, edge_ego, sense, t, uv),
             "eval edge face uv");
    printf("edge %d t %.18e face %d uv %.18e %.18e\n", edgeid, t, faceid, uv[0],
           uv[1]);
    t = trange[1];
    is_equal(EGADS_SUCCESS, EG_getEdgeUV(face_ego, edge_ego, sense, t, uv),
             "eval edge face uv");
    printf("edge %d t %.18e face %d uv %.18e %.18e\n", edgeid, t, faceid, uv[0],
           uv[1]);
  }
  return 0;
}
