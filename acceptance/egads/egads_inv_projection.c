/*

gcc-7  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized
-I/Users/mpark/esp/EngSketchPad/include -o egads_inv_projection
egads_inv_projection.c -Wl,-rpath,/Users/mpark/esp/EngSketchPad/lib
-L/Users/mpark/esp/EngSketchPad/lib -legads   -lm

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
  double input_xyz[3];
  double param[2];
  double output_xyz[3];

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 3) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS,
           EG_loadModel(context, 0, "onera-m6-sharp-te.egads", &model),
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

  faceid = 2;

  input_xyz[0] = 25.0;
  input_xyz[1] = 0.0;
  input_xyz[2] = 25.0;
  is_equal(EGADS_SUCCESS,
           EG_invEvaluate(faces[faceid - 1], input_xyz, param, output_xyz),
           "EG inv  eval");

  printf(" input %f %f %f\n", input_xyz[0], input_xyz[1], input_xyz[2]);
  printf("output %f %f %f\n", output_xyz[0], output_xyz[1], output_xyz[2]);
  printf(" param %f %f\n", param[0], param[1]);

  input_xyz[0] = 25.0;
  input_xyz[1] = 0.0;
  input_xyz[2] = -25.0;
  is_equal(EGADS_SUCCESS,
           EG_invEvaluate(faces[faceid - 1], input_xyz, param, output_xyz),
           "EG inv  eval");

  printf(" input %f %f %f\n", input_xyz[0], input_xyz[1], input_xyz[2]);
  printf("output %f %f %f\n", output_xyz[0], output_xyz[1], output_xyz[2]);
  printf(" param %f %f\n", param[0], param[1]);

  return 0;
}
