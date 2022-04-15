/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/ESPbeta-2022-02-02/include -o egads_loop_dup \
egads_loop_dup.c -Wl,-rpath,${HOME}/local/pkgs/ESPbeta-2022-02-02/lib \
-L${HOME}/local/pkgs/ESPbeta-2022-02-02/lib -legads   -lm

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

int main(int argc, char *argv[]) {
  ego context;
  ego body = NULL;
  ego geom, *bodies;
  int oclass, nbody, *senses, mtype;
  int nface;
  ego *faces, face;
  ego surface, *loops;
  int faceclass, facetype, nloop, *loopsenses;
  double facebounds[4];
  ego loop_ref;
  int iloop;
  int loopclass, looptype, nchild, *children_senses;
  ego *children;

  is_equal(2, argc, "usage: egads_loop_dup project.egads");

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 3) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, argv[1], &body), "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(body, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses),
           "EG topo body");
  is_equal(1, nbody, "expected 1 body");
  body = bodies[0];

  is_equal(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, FACE, &nface, &faces),
           "EG face topo");
  is_equal(1, nface, "expected 1 face");
  face = faces[0];

  is_equal(EGADS_SUCCESS,
           EG_getTopology(face, &surface, &faceclass, &facetype, facebounds,
                          &nloop, &loops, &loopsenses),
           "topo");
  is_equal(1, nloop, "expected 1 loop");

  iloop = 0;
  is_equal(EGADS_SUCCESS,
           EG_getTopology(loops[iloop], &loop_ref, &loopclass, &looptype, NULL,
                          &nchild, &children, &children_senses),
           "topo");
  printf("loop_ref %p surface %p\n", (void *)loop_ref, (void *)surface);
  return 0;
}
