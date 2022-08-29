/*

gcc  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I${HOME}/local/pkgs/ESP121/EngSketchPad/include -o egads_edge \
egads_edge.c -Wl,-rpath,${HOME}/local/pkgs/ESP121/EngSketchPad/lib \
-L${HOME}/local/pkgs/ESP121/EngSketchPad/lib -legads -lm \
&& ./egads_edge

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
  ego model;
  size_t cad_data_size;
  char *cad_data;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  is_equal(EGADS_SUCCESS,
           EG_loadModel(context, 0, "extract-face1.egads", &model), "EG load");

  /* entry point NOT in egads.h */
  int EG_exportModel(ego mobject, size_t * nbytes, char *stream[]);

  is_equal(EGADS_SUCCESS, EG_exportModel(model, &cad_data_size, &cad_data),
           "EG stream");
  {
    FILE *file;
    file = fopen("extract-face1.bin", "w");
    if (NULL == (void *)file) printf("unable to open file\n");
    is_equal(1, fwrite(&cad_data_size, sizeof(size_t), 1, file),
             "write egadslite data size");
    is_equal(cad_data_size, fwrite(cad_data, sizeof(char), cad_data_size, file),
             "write egadslite data");
    fclose(file);
  }
  {
    ego geom, *children, body, *edges;
    int oclass, mtype, nchild, *senses;
    int nedge, i;
    double params[2], eval[18];
    is_equal(EGADS_SUCCESS,
             EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nchild,
                            &children, &senses),
             "EG topo bodies");
    printf("oclass %d mtype %d nchild %d\n", oclass, mtype, nchild);
    is_equal(1, nchild, "expected 1 body");
    body = children[0];
    is_equal(EGADS_SUCCESS, EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges),
             "EG edge topo");
    printf("nedge %d\n", nedge);
    params[0] = 0.25;
    is_equal(EGADS_SUCCESS, EG_evaluate(edges[0], params, eval), "eval");
    for (i = 0; i < 9; i++) {
      printf("%d %f\n", i, eval[i]);
    }
    EG_free(edges);
  }

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");

  return 0;
}
