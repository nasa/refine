/*

gcc-10  -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized \
-I/Users/mpark/local/pkgs/EGADS/trunk/include -o egads_effective_edgenode \
egads_effective_edgenode.c -Wl,-rpath,/Users/mpark/local/pkgs/EGADS/trunk/lib \
-L/Users/mpark/local/pkgs/EGADS/trunk/lib -legads   -lm \
&& ./egads_effective_edgenode

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

int main(void) {
  ego context;
  ego model = NULL;
  ego geom, *bodies;
  int oclass, mtype, nbody, *senses;
  ego newBodies[3], newModel;
  double params[3];
  int tess_status, nvert;
  double angle;
  ego solid;

  is_equal(EGADS_SUCCESS, EG_open(&context), "EG open");
  /* Success returns the old output level. (0-silent to 3-debug) */
  is_true(EG_setOutLevel(context, 2) >= 0, "make verbose");

  {
    ego box1, box2;
    {
      int stype = BOX;
      double data[] = {-0.5, -0.5, -0.5, 1.0, 1.0, 0.5};
      is_equal(EGADS_SUCCESS, EG_makeSolidBody(context, stype, data, &box1),
               "make solid body");
    }
    {
      int stype = BOX;
      double offset = 0.01;
      double data[] = {-0.5 + offset, -0.5, 0.0, 1.0, 1.0, 0.5};
      is_equal(EGADS_SUCCESS, EG_makeSolidBody(context, stype, data, &box2),
               "make solid body");
    }
    is_equal(EGADS_SUCCESS, EG_generalBoolean(box1, box2, FUSION, 0.1, &model),
             "make solid body");
  }
  is_equal(EGADS_SUCCESS,
           EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                          &senses),
           "EG topo bodies");
  printf("oclass %d mtype %d nbody %d\n", oclass, mtype, nbody);
  is_equal(1, nbody, "expected 1 body");

  remove("egads_effective_edgenode_orig.egads");
  is_equal(EGADS_SUCCESS,
           EG_saveModel(bodies[0], "egads_effective_edgenode_orig.egads"),
           "EG save eff");

  /* copy the Body so we can use/save it later */
  is_equal(EGADS_SUCCESS, EG_copyObject(bodies[0], NULL, &newBodies[0]),
           "EG copy object");
  EG_deleteObject(model);

  /* make the tessellation object */
  params[0] = 0.1;
  params[1] = 0.01;
  params[2] = 20.0;

  is_equal(EGADS_SUCCESS, EG_makeTessBody(newBodies[0], params, &newBodies[1]),
           "EG tess");
  is_equal(EGADS_SUCCESS,
           EG_statusTessBody(newBodies[1], &geom, &tess_status, &nvert),
           "EG tess");
  is_equal(1, tess_status, "tess not closed");

  /* make the effective topology object */
  angle = 10.0;
  is_equal(EGADS_SUCCESS, EG_initEBody(newBodies[1], angle, &newBodies[2]),
           "initEB");
  is_equal(EGADS_SUCCESS, EG_finishEBody(newBodies[2]), "finEB");

  /* make the model with the body, tessellation and effective topology body
     notes: 1) mtype = 3 is the total number of objects in the model
               Body Objects must be first
            2) nchild = 1 are the number of actual Body Objects */
  is_equal(EGADS_SUCCESS,
           EG_makeTopology(context, NULL, MODEL, 3, NULL, 1, newBodies, NULL,
                           &newModel),
           "make Topo Model");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(newModel, &geom, &oclass, &mtype, NULL, &nbody,
                          &bodies, &senses),
           "EG topo bodies");
  printf("oclass %d mtype %d nbody %d\n", oclass, mtype, nbody);

  remove("egads_effective_edgenode_eff.egads");
  is_equal(EGADS_SUCCESS,
           EG_saveModel(newModel, "egads_effective_edgenode_eff.egads"),
           "EG save eff");

  is_equal(EGADS_SUCCESS,
           EG_getTopology(newModel, &geom, &oclass, &mtype, NULL, &nbody,
                          &bodies, &senses),
           "EG topo bodies");
  printf("oclass %d mtype %d nbody %d\n", oclass, mtype, nbody);
  is_equal(MODEL, oclass, "not model");

  solid = NULL;
  {
    int ibody;
    int bodyclass, bodytype;
    ego owner, prev, next;
    for (ibody = 0; ibody < mtype; ibody++) {
      is_equal(EGADS_SUCCESS,
               EG_getInfo(bodies[ibody], &bodyclass, &bodytype, &owner, &prev,
                          &next),
               "info");
      if (BODY == bodyclass) printf("BODY ");
      if (TESSELLATION == bodyclass) printf("TESSELLATION ");
      if (EBODY == bodyclass) printf("EBODY ");
      printf("body %d oclass %d mtype %d\n", ibody, bodyclass, bodytype);
      if (EBODY == bodyclass) solid = bodies[ibody];
    }
  }

  {
    int nnode, nedge, nface;
    ego *nodes, *edges, *faces;
    ego node, edge;
    double t, node_eval[18], edge_eval[18], dist;
    double trange[2];
    
    is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, NODE, &nnode, &nodes),
             "EG edge topo");
    is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EEDGE, &nedge, &edges),
             "EG edge topo");
    is_equal(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EFACE, &nface, &faces),
             "EG face topo");
    printf("nnode %d effective nedge %d nface %d\n", nnode, nedge, nface);

    /* edge 7 t=0 node 2 */
    node = nodes[2 - 1];
    edge = edges[7 - 1];
    t = 0;

      {
        ego tempref, *tempchldrn;
        int tempclass, temptype, tempchild, *tempsens;
        is_equal(EGADS_SUCCESS,
             EG_getTopology(node, &tempref, &tempclass, &temptype, node_eval, &tempchild,
                            &tempchldrn, &tempsens),
             "EG topo node");
        is_equal(EGADS_SUCCESS,
             EG_getTopology(edge, &tempref, &tempclass, &temptype, trange, &tempchild,
                            &tempchldrn, &tempsens),
             "EG topo edge");
      }
    is_equal(EGADS_SUCCESS, EG_evaluate(edge, &t, edge_eval), "EG eval edge");

      dist = sqrt(pow(node_eval[0] - edge_eval[0], 2) +
                  pow(node_eval[1] - edge_eval[1], 2) +
                  pow(node_eval[2] - edge_eval[2], 2));
      printf("t %.3f trange %f %f dist %.3e node x=%.3f edge x=%.3f\n", t,
             trange[0], trange[1], dist, node_eval[0], edge_eval[0]);

    /* edge 13 t=0 node 3 */
    node = nodes[3 - 1];
    edge = edges[13 - 1];
    t = 0;

      {
        ego tempref, *tempchldrn;
        int tempclass, temptype, tempchild, *tempsens;
        is_equal(EGADS_SUCCESS,
             EG_getTopology(node, &tempref, &tempclass, &temptype, node_eval, &tempchild,
                            &tempchldrn, &tempsens),
             "EG topo node");
        is_equal(EGADS_SUCCESS,
             EG_getTopology(edge, &tempref, &tempclass, &temptype, trange, &tempchild,
                            &tempchldrn, &tempsens),
             "EG topo edge");
      }
    is_equal(EGADS_SUCCESS, EG_evaluate(edge, &t, edge_eval), "EG eval edge");

      dist = sqrt(pow(node_eval[0] - edge_eval[0], 2) +
                  pow(node_eval[1] - edge_eval[1], 2) +
                  pow(node_eval[2] - edge_eval[2], 2));
      printf("t %.3f trange %f %f dist %.3e node x=%.3f edge x=%.3f\n", t,
             trange[0], trange[1], dist, node_eval[0], edge_eval[0]);

    
    EG_free(faces);
    EG_free(edges);
    EG_free(nodes);
  }

  is_equal(EGADS_SUCCESS, EG_close(context), "EG close");
  return 0;
}
