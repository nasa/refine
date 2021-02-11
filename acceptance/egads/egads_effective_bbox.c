/*

gcc-10  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/EGADS/trunk/include -o egads_effective_bbox \
egads_effective_bbox.c -Wl,-rpath,${HOME}/local/pkgs/EGADS/trunk/lib \
-L${HOME}/local/pkgs/EGADS/trunk/lib -legads   -lm \
&& ./egads_effective_bbox

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
  int oclass, nego, nbody, *senses;
  ego ebody = NULL;

  if (argc < 2) {
    printf("usage: \n %s project-eff.egads\n", argv[0]);
    return 1;
  }

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS, EG_loadModel(context, 0, argv[1], &model), "EG load");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(model, &geom, &oclass, &nego, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
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
        ebody = bodies[ibody];
        printf("EBODY extracted from model\n");
      }
    }
  }

  {
    ego *faces;
    int nface, iface;
    ego face;
    double diag, box[6];
    is_equal(EGADS_SUCCESS, EG_getBodyTopos(ebody, NULL, EFACE, &nface, &faces),
             "EG face topo");
    iface = 16;
    face = faces[iface - 1];
    is_equal(EGADS_SUCCESS, EG_getBoundingBox(face, box), "EG bounding box");
    diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                (box[1] - box[4]) * (box[1] - box[4]) +
                (box[2] - box[5]) * (box[2] - box[5]));
    printf("face %d diagoal %e\n", iface, diag);
    printf("face %d bbox %e %e %e %e %e %e\n", iface, box[0], box[1], box[2],
           box[3], box[4], box[5]);
    EG_free(faces);

    is_equal(EGADS_SUCCESS, EG_getBoundingBox(ebody, box), "EG bounding box");
    diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
                (box[1] - box[4]) * (box[1] - box[4]) +
                (box[2] - box[5]) * (box[2] - box[5]));
    printf("ebody diagoal %e\n", diag);
    printf("ebody bbox %e %e %e %e %e %e\n", box[0], box[1], box[2], box[3],
           box[4], box[5]);
  }

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");
  return 0;
}
