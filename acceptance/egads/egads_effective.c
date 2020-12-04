/*

gcc-10  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized
-I/Users/mpark/local/pkgs/EGADS/trunk/include -o egads_effective
egads_effective.c -Wl,-rpath,/Users/mpark/local/pkgs/EGADS/trunk/lib
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
  ego geom, *bodies;
  int oclass, mtype, nbody, *senses;
  ego solid;
  double params[3];
  ego tess;
  int tess_status, nvert;
  double angle;
  ego ebody;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 3) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, "boxbox.egads", &model),
           "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
  is_equal(1, nbody, "expected 1 body");
  solid = bodies[0];

  params[0] = 0.1;
  params[1] = 0.01;
  params[2] = 20.0;

  is_equal(EGADS_SUCCESS, EG_makeTessBody(solid, params, &tess), "EG tess");
  is_equal(EGADS_SUCCESS, EG_statusTessBody(tess, &geom, &tess_status, &nvert),
           "EG tess");
  is_equal(1, tess_status, "tess not closed");

  angle = 10.0;
  is_equal(EGADS_SUCCESS, EG_virtualize(tess, angle, &ebody), "virt");
  is_equal(EGADS_SUCCESS, EG_finalize(ebody), "fin")

      is_equal(EGADS_SUCCESS, EG_saveModel(model, "boxboxeff.egads"),
               "EG save eff");

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");
  return 0;
}
