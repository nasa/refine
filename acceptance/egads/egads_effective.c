/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/EGADS/trunk/include -o egads_effective \
egads_effective.c -Wl,-rpath,${HOME}/local/pkgs/EGADS/trunk/lib \
-L${HOME}/local/pkgs/EGADS/trunk/lib -legads   -lm \
&& ./egads_effective

*/

#include <math.h>
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
  ego model = NULL;
  ego geom, *bodies;
  int oclass, mtype, nbody, *senses;
  ego newBodies[2], tess;
  double params[3], diag, box[6];
  int tess_status, nvert;
  double angle;
  int neface;
  ego *efaces;

  if (argc < 2) {
    printf("usage: \n %s project.egads\n", argv[0]);
    return 1;
  }

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, argv[1], &model), "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
  is_equal(1, nbody, "expected 1 body");

  /* copy the Body so we can use/save it later */
  is_equal(EGADS_SUCCESS, EG_copyObject(bodies[0], NULL, &newBodies[0]),
           "EG copy object");
  EG_deleteObject(model);

  /* make the tessellation object */
  is_equal(EGADS_SUCCESS, EG_getBoundingBox(newBodies[0], box), "EG bbox");
  diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
              (box[1] - box[4]) * (box[1] - box[4]) +
              (box[2] - box[5]) * (box[2] - box[5]));

  params[0] = 0.1 * diag;
  params[1] = 0.01 * diag;
  params[2] = 20.0;

  is_equal(EGADS_SUCCESS, EG_makeTessBody(newBodies[0], params, &tess),
           "EG tess");
  is_equal(EGADS_SUCCESS, EG_statusTessBody(tess, &geom, &tess_status, &nvert),
           "EG tess");
  is_equal(1, tess_status, "tess not closed");

  /* make the effective topology object */
  angle = 10.0;
  is_equal(EGADS_SUCCESS, EG_initEBody(tess, angle, &newBodies[1]), "initEB");
  is_equal(EGADS_SUCCESS,
           EG_makeAttrEFaces(newBodies[1], "group", &neface, &efaces), "finEB");
  EG_free(efaces);
  is_equal(EGADS_SUCCESS, EG_finishEBody(newBodies[1]), "finEB");

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");
  return 0;
}
