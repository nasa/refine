/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/EGADS/trunk/include -o egads_copy_leak \
egads_copy_leak.c -Wl,-rpath,${HOME}/local/pkgs/EGADS/trunk/lib \
-L${HOME}/local/pkgs/EGADS/trunk/lib -legads   -lm \
&& valgrind --leak-check=full --show-leak-kinds=all ./egads_copy_leak

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
  ego box1, box2;
  ego body;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  {
    int stype = BOX;
    double data[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 0.5};
    is_equal(EGADS_SUCCESS, EG_makeSolidBody(context, stype, data, &box1),
             "make solid body");
  }
  {
    int stype = BOX;
    double data[6] = {0.0, 0.0, 0.5, 1.0, 1.0, 0.5};
    is_equal(EGADS_SUCCESS, EG_makeSolidBody(context, stype, data, &box2),
             "make solid body");
  }
  {
    ego model, geom, *bodies;
    int oclass, nego, nbody, *senses;
    is_equal(EGADS_SUCCESS, EG_generalBoolean(box1, box2, FUSION, 0.0, &model),
             "make solid body");
    is_equal(EGADS_SUCCESS,
             EG_getTopology(model, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                            &senses),
             "EG topo bodies");
    is_equal(1, nbody, "expected 1 body");
    is_equal(EGADS_SUCCESS, EG_copyObject(bodies[0], NULL, &body), "copy body");
    is_equal(0, EG_deleteObject(model), "delete temp model");
    is_equal(0, EG_deleteObject(box1), "delete box1");
    is_equal(0, EG_deleteObject(box2), "delete box2");
  }

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");

  return 0;
}
