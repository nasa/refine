/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/EGADS/trunk/include -o egads_save \
egads_save.c -Wl,-rpath,${HOME}/local/pkgs/EGADS/trunk/lib \
-L${HOME}/local/pkgs/EGADS/trunk/lib -legads   -lm \
&& ./egads_save

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
  ego box;
  ego model;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  {
    int stype = BOX;
    double data[] = {-0.5, -0.5, -0.5, 1.0, 1.0, 0.5};
    is_equal(EGADS_SUCCESS, EG_makeSolidBody(context, stype, data, &box),
             "make solid body");
  }

  is_equal(
      EGADS_SUCCESS,
      EG_makeTopology(context, NULL, MODEL, 0, NULL, 1, &box, NULL, &model),
      "make Topo Model");

  remove("egads_save_box.egads");
  is_equal(EGADS_SUCCESS, EG_saveModel(model, "egads_save_box.egads"),
           "EG save");

  EG_free(model);
  model = NULL;

  is_equal(EGADS_SUCCESS,
           EG_loadModel(context, 0, "egads_save_box.egads", &model), "EG load");
  {
    ego geom, *children;
    int oclass, mtype, nchild, *senses;
    is_equal(EGADS_SUCCESS,
             EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nchild,
                            &children, &senses),
             "EG topo bodies");
    printf("oclass %d mtype %d nchild %d\n", oclass, mtype, nchild);
  }

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");
  return 0;
}
