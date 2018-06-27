
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_EGADS
#include "egads.h"
#endif

#include "ref_geom.h"

#include "ref_cell.h"
#include "ref_export.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_node.h"

#include "ref_dict.h"
#include "ref_sort.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"

REF_STATUS ref_geom_initialize(REF_GEOM ref_geom) {
  REF_INT geom;
  ref_geom_n(ref_geom) = 0;
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++) {
    ref_geom_type(ref_geom, geom) = REF_EMPTY;
    ref_geom_id(ref_geom, geom) = geom + 1;
  }
  ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
  ref_geom_blank(ref_geom) = 0;
  if (NULL != (void *)(ref_geom->ref_adj))
    RSS(ref_adj_free(ref_geom->ref_adj), "free to prevent leak");
  RSS(ref_adj_create(&(ref_geom->ref_adj)), "create ref_adj for ref_geom");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_create(REF_GEOM *ref_geom_ptr) {
  REF_GEOM ref_geom;

  (*ref_geom_ptr) = NULL;

  ref_malloc(*ref_geom_ptr, 1, REF_GEOM_STRUCT);

  ref_geom = (*ref_geom_ptr);

  ref_geom_max(ref_geom) = 10;

  ref_malloc(ref_geom->descr, 3 * ref_geom_max(ref_geom), REF_INT);
  ref_malloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);
  ref_geom->ref_adj = (REF_ADJ)NULL;
  RSS(ref_geom_initialize(ref_geom), "init geom list");

  ref_geom->uv_area_sign = NULL;
  ref_geom->segments_per_radian_of_curvature = 10.0;

  ref_geom->nnode = REF_EMPTY;
  ref_geom->nedge = REF_EMPTY;
  ref_geom->nface = REF_EMPTY;
  ref_geom->context = NULL;
#ifdef HAVE_EGADS
  {
    ego context;
    REIS(EGADS_SUCCESS, EG_open(&context), "EG open");
    /* Success returns the old output level. (0-silent to 3-debug) */
    RAS(EG_setOutLevel(context, 0) >= 0, "make silent");
    ref_geom->context = (void *)context;
  }
#endif
  ref_geom->solid = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_free(REF_GEOM ref_geom) {
  if (NULL == (void *)ref_geom) return REF_NULL;
  ref_free(ref_geom->cad_data);
#ifdef HAVE_EGADS
  if (NULL != ref_geom->faces) EG_free((ego *)(ref_geom->faces));
  if (NULL != ref_geom->edges) EG_free((ego *)(ref_geom->edges));
  if (NULL != ref_geom->nodes) EG_free((ego *)(ref_geom->nodes));
  /* solid is not freeable */
  if (NULL != ref_geom->context)
    REIS(EGADS_SUCCESS, EG_close((ego)(ref_geom->context)), "EG close");
#endif
  RSS(ref_adj_free(ref_geom->ref_adj), "adj free");
  ref_free(ref_geom->uv_area_sign);
  ref_free(ref_geom->param);
  ref_free(ref_geom->descr);
  ref_free(ref_geom);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_deep_copy(REF_GEOM *ref_geom_ptr, REF_GEOM original) {
  REF_GEOM ref_geom;
  REF_INT geom, i;
  (*ref_geom_ptr) = NULL;

  ref_malloc(*ref_geom_ptr, 1, REF_GEOM_STRUCT);

  ref_geom = (*ref_geom_ptr);

  ref_geom_n(ref_geom) = ref_geom_n(original);
  ref_geom_max(ref_geom) = ref_geom_max(original);

  ref_malloc(ref_geom->descr, 3 * ref_geom_max(ref_geom), REF_INT);
  ref_malloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);
  ref_geom->uv_area_sign = NULL;
  ref_geom->segments_per_radian_of_curvature =
      original->segments_per_radian_of_curvature;

  for (geom = 0; geom < ref_geom_max(ref_geom); geom++)
    for (i = 0; i < 3; i++)
      ref_geom_descr(ref_geom, i, geom) = ref_geom_descr(original, i, geom);
  ref_geom_blank(ref_geom) = ref_geom_blank(original);
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++)
    for (i = 0; i < 2; i++)
      ref_geom_param(ref_geom, i, geom) = ref_geom_param(original, i, geom);

  RSS(ref_adj_deep_copy(&(ref_geom->ref_adj), original->ref_adj),
      "deep copy ref_adj for ref_geom");

  ref_geom->nnode = REF_EMPTY;
  ref_geom->nedge = REF_EMPTY;
  ref_geom->nface = REF_EMPTY;
  ref_geom->context = NULL;
  ref_geom->solid = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_pack(REF_GEOM ref_geom, REF_INT *o2n) {
  REF_INT geom, new, i;
  new = 0;
  each_ref_geom(ref_geom, geom) {
    for (i = 0; i < 2; i++)
      ref_geom_descr(ref_geom, i, new) = ref_geom_descr(ref_geom, i, geom);
    ref_geom_descr(ref_geom, 2, new) = o2n[ref_geom_descr(ref_geom, 2, geom)];
    for (i = 0; i < 2; i++)
      ref_geom_param(ref_geom, i, new) = ref_geom_param(ref_geom, i, geom);
    new ++;
  }
  REIS(new, ref_geom_n(ref_geom), "count mismatch");
  if (ref_geom_n(ref_geom) < ref_geom_max(ref_geom)) {
    for (geom = ref_geom_n(ref_geom); geom < ref_geom_max(ref_geom); geom++) {
      ref_geom_type(ref_geom, geom) = REF_EMPTY;
      ref_geom_id(ref_geom, geom) = geom + 1;
    }
    ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
    ref_geom_blank(ref_geom) = ref_geom_n(ref_geom);
  } else {
    ref_geom_blank(ref_geom) = REF_EMPTY;
  }
  RSS(ref_adj_free(ref_geom->ref_adj), "free to prevent leak");
  RSS(ref_adj_create(&(ref_geom->ref_adj)), "create ref_adj for ref_geom");

  each_ref_geom(ref_geom, geom) {
    RSS(ref_adj_add(ref_geom->ref_adj, ref_geom_node(ref_geom, geom), geom),
        "register geom");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_save(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  FILE *file;
  REF_INT geom, global;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  RSS(ref_node_synchronize_globals(ref_node), "sync");
  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");
  fprintf(file, "%d\n", ref_cell_n(ref_cell));
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    nodes[0] = ref_node_global(ref_node, nodes[0]);
    nodes[1] = ref_node_global(ref_node, nodes[1]);
    fprintf(file, "%d %d %d\n", nodes[0], nodes[1], nodes[2]);
  }
  fprintf(file, "%d\n", ref_geom_n(ref_geom));
  each_ref_geom(ref_geom, geom) {
    global = ref_node_global(ref_node, ref_geom_node(ref_geom, geom));
    fprintf(file, "%d %d %d %.18e %.18e\n", ref_geom_type(ref_geom, geom),
            ref_geom_id(ref_geom, geom), global,
            ref_geom_param(ref_geom, 0, geom),
            ref_geom_param(ref_geom, 1, geom));
  }
  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_load(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT edge, nedge;
  REF_INT new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT geom, ngeom;
  REF_INT global, node, gtype, id;
  REF_DBL param[2];
  RSS(ref_node_synchronize_globals(ref_node), "sync");
  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");
  REIS(1, fscanf(file, "%d", &nedge), "nedge");
  for (edge = 0; edge < nedge; edge++) {
    REIS(3, fscanf(file, "%d %d %d", &(nodes[0]), &(nodes[1]), &(nodes[2])),
         "edge nodes");
    RSS(ref_node_local(ref_node, nodes[0], &node), "g2l");
    nodes[0] = node;
    RSS(ref_node_local(ref_node, nodes[1], &node), "g2l");
    nodes[1] = node;
    RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell), "new edge");
  }
  REIS(1, fscanf(file, "%d", &ngeom), "ngeom");
  for (geom = 0; geom < ngeom; geom++) {
    REIS(5,
         fscanf(file, "%d %d %d %lf %lf\n", &gtype, &id, &global, &(param[0]),
                &(param[1])),
         "gr");
    RSS(ref_node_local(ref_node, global, &node), "g2l");
    RSS(ref_geom_add(ref_geom, node, gtype, id, param), "add geom");
  }
  if (!feof(file)) {
    printf("expected end of %s file after %d edge %d geom\n", filename, nedge,
           ngeom);
    fclose(file);
    return REF_FAILURE;
  }
  fclose(file);
  RSS(ref_geom_uv_area_report(ref_grid), "report");
  return REF_SUCCESS;
}

REF_STATUS ref_geom_node_faces(REF_GRID ref_grid, REF_ADJ *ref_adj_arg) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj;
  REF_INT *e2f, id, toponode;
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double trange[2];
  RSS(ref_adj_create(ref_adj_arg), "create ref_adj");
  ref_adj = *ref_adj_arg;
  RSS(ref_geom_edge_faces(ref_grid, &e2f), "edge2face");
  for (id = 1; id <= ref_geom->nedge; id++) {
    object = ((ego *)(ref_geom->edges))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");

    if (0 < nchild) {
      toponode = EG_indexBodyTopo(ref_geom->solid, pchldrn[0]);
      if (0 < e2f[0 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[0 + 2 * (id - 1)]),
            "add");
      }
      if (0 < e2f[1 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[1 + 2 * (id - 1)]),
            "add");
      }
    }

    if (1 < nchild) {
      toponode = EG_indexBodyTopo(ref_geom->solid, pchldrn[1]);
      if (0 < e2f[0 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[0 + 2 * (id - 1)]),
            "add");
      }
      if (0 < e2f[1 + 2 * (id - 1)]) {
        RSS(ref_adj_add_uniquely(ref_adj, toponode, e2f[1 + 2 * (id - 1)]),
            "add");
      }
    }
  }
  ref_free(e2f);
  return REF_SUCCESS;
#else
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (REF_EMPTY == ref_geom->nnode) printf("No EGADS loaded\n");
  *ref_adj_arg = (REF_ADJ)NULL;
  printf("No EGADS linked for %s\n", __func__);

  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_edge_faces(REF_GRID ref_grid, REF_INT **edge_face_arg) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT *e2f, *nface;
  REF_INT face, edge;
  ref_malloc_init(*edge_face_arg, 2 * (ref_geom->nedge), REF_INT, REF_EMPTY);
  e2f = *edge_face_arg;
  ref_malloc_init(nface, (ref_geom->nedge), REF_INT, 0);
  printf("nedge %d nface %d\n", ref_geom->nedge, ref_geom->nface);
  for (face = 0; face < (ref_geom->nface); face++) {
    ego esurf, *eloops, eref;
    int oclass, mtype, nloop, *senses, *pinfo;
    double data[18], *preal;
    ego ecurve, *eedges;
    int iloop, iedge, nedge;

    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->faces))[face], &esurf, &oclass,
                        &mtype, data, &nloop, &eloops, &senses),
         "topo");
    printf("faceid %d uv ([%f,%f],[%f,%f])\n", face + 1, data[0], data[1],
           data[2], data[3]);
    REIS(EGADS_SUCCESS,
         EG_getGeometry(esurf, &oclass, &mtype, &eref, &pinfo, &preal), "geom");
    EG_free(pinfo);
    EG_free(preal);
    for (iloop = 0; iloop < nloop; iloop++) {
      /* loop through all Edges associated with this Loop */
      REIS(EGADS_SUCCESS,
           EG_getTopology(eloops[iloop], &ecurve, &oclass, &mtype, data, &nedge,
                          &eedges, &senses),
           "topo");
      for (iedge = 0; iedge < nedge; iedge++) {
        edge = EG_indexBodyTopo((ego)(ref_geom->solid), eedges[iedge]) - 1;
        e2f[nface[edge] + 2 * edge] = face + 1;
        nface[edge]++;
      }
    }
  }

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    REF_INT i;
    printf("edge %4d has %2d faces:", edge + 1, nface[edge]);
    for (i = 0; i < nface[edge]; i++) printf(" %4d", e2f[i + 2 * edge]);
    printf("\n");
  }
  ref_free(nface);
  return REF_SUCCESS;
#else
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (REF_EMPTY == ref_geom->nnode) printf("No EGADS loaded\n");
  *edge_face_arg = (REF_INT *)NULL;
  printf("No EGADS linked for %s\n", __func__);

  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_face_surface(REF_GEOM ref_geom, REF_INT faceid,
                                 REF_GEOM_SURFACE *surface) {
#ifdef HAVE_EGADS
  ego esurf, *eloops, eref;
  int oclass, mtype, nloop, *senses, *pinfo;
  double data[18], *preal;

  RNS(ref_geom->faces, "faces not loaded");
  if (faceid < 1 || faceid > ref_geom->nface) return REF_INVALID;

  REIS(EGADS_SUCCESS,
       EG_getTopology(((ego *)(ref_geom->faces))[faceid - 1], &esurf, &oclass,
                      &mtype, data, &nloop, &eloops, &senses),
       "topo");
  REIS(EGADS_SUCCESS,
       EG_getGeometry(esurf, &oclass, &mtype, &eref, &pinfo, &preal), "geom");
  EG_free(pinfo);
  EG_free(preal);
  switch (mtype) {
    case PLANE:
      *surface = REF_GEOM_PLANE;
      break;
    case SPHERICAL:
      *surface = REF_GEOM_SPHERICAL;
      break;
    case CYLINDRICAL:
      *surface = REF_GEOM_CYLINDRICAL;
      break;
    case REVOLUTION:
      *surface = REF_GEOM_REVOLUTION;
      break;
    case TOROIDAL:
      *surface = REF_GEOM_TOROIDAL;
      break;
    case TRIMMED:
      *surface = REF_GEOM_TRIMMED;
      break;
    case BEZIER:
      *surface = REF_GEOM_BEZIER;
      break;
    case BSPLINE:
      *surface = REF_GEOM_BSPLINE;
      break;
    case OFFSET:
      *surface = REF_GEOM_OFFSET;
      break;
    case CONICAL:
      *surface = REF_GEOM_CONICAL;
      break;
    case EXTRUSION:
      *surface = REF_GEOM_EXTRUSION;
      break;
    default:
      printf("mtype %d\n", mtype);
      RSS(REF_IMPLEMENT, "unknown surface type");
  }
  return REF_SUCCESS;
#else
  if (REF_EMPTY == ref_geom->nnode) printf("No EGADS loaded face %d\n", faceid);
  *surface = REF_GEOM_SURFACE_LAST;
  printf("No EGADS linked for %s\n", __func__);
  return REF_IMPLEMENT;
#endif
}

#ifdef HAVE_EGADS
static REF_STATUS ref_geom_recon_nodes(REF_GRID ref_grid, REF_INT **cad_nodes) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_ADJ n2f;
  REF_INT id, item, best_node, node, faceid;
  REF_DBL best_dist, dist;
  REF_INT max_faceids = 50;
  REF_INT *grid_faceids, grid_nfaceids;
  REF_INT *cad_faceids, cad_nfaceids;
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double xyz[3];
  REF_BOOL show_xyz = REF_FALSE;
  REF_INT i, j;
  REF_BOOL found, all_found;
  RSS(ref_geom_node_faces(ref_grid, &n2f), "build n2f");
  ref_malloc(grid_faceids, max_faceids, REF_INT);
  ref_malloc(cad_faceids, max_faceids, REF_INT);
  printf("searching for %d topo nodes\n", ref_geom->nnode);
  ref_malloc(*cad_nodes, ref_geom->nnode, REF_INT);
  for (id = 1; id <= ref_geom->nnode; id++) {
    object = ((ego *)(ref_geom->nodes))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, xyz, &nchild, &pchldrn,
                        &psens),
         "EG topo node");
    cad_nfaceids = 0;
    each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
      cad_faceids[cad_nfaceids] = faceid;
      cad_nfaceids++;
    }
    best_node = REF_EMPTY;
    best_dist = 1.0e20;
    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node, max_faceids,
                                  &grid_nfaceids, grid_faceids),
          "count faceids");
      if (grid_nfaceids != cad_nfaceids) continue;
      all_found = REF_TRUE;
      for (i = 0; i < cad_nfaceids; i++) {
        found = REF_FALSE;
        for (j = 0; j < grid_nfaceids; j++) {
          found = found || (cad_faceids[i] == grid_faceids[j]);
        }
        all_found = all_found && found;
      }
      if (!all_found) continue;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      if (dist < best_dist) {
        best_node = node;
        best_dist = dist;
      }
    }
    if (REF_EMPTY == best_node) {
      printf("cad node %d expects faces", id);
      each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
        printf(" %d", faceid);
      }
      printf("\n");
      printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      THROW("unable to find a grid vertex that matches cad topo");
    }
    printf(" topo node id %3d node %6d dist %.4e fid", id, best_node,
           best_dist);
    RSS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node, max_faceids,
                                &grid_nfaceids, grid_faceids),
        "count faceids");
    for (i = 0; i < grid_nfaceids; i++) printf(" %d", grid_faceids[i]);
    printf(" expects");
    each_ref_adj_node_item_with_ref(n2f, id, item, faceid) {
      printf(" %d", faceid);
    }
    printf("\n");
    if (show_xyz) {
      printf(" d %23.15e %23.15e %23.15e\n",
             ref_node_xyz(ref_node, 0, best_node),
             ref_node_xyz(ref_node, 1, best_node),
             ref_node_xyz(ref_node, 2, best_node));
      printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
    }
    (*cad_nodes)[id - 1] = best_node;
    RSS(ref_geom_add(ref_geom, best_node, REF_GEOM_NODE, id, NULL), "node");
  }
  ref_free(cad_faceids);
  ref_free(grid_faceids);
  ref_adj_free(n2f);
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_geom_recon(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double xyz[9], trange[2];
  REF_INT node, id, best_node;
  REF_DBL best_dist, dist, best_param;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *cad_nodes;
  double t;
  double param[2];
  REF_INT i, cell, edge_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pass, updates;
  REF_BOOL show_xyz = REF_FALSE;
  REF_BOOL debug = REF_FALSE;
#define REF_GEOM_MAX_FACEIDS (50)
  REF_INT nfaceid, faceids[REF_GEOM_MAX_FACEIDS];
  REF_INT nfaceid0, faceids0[REF_GEOM_MAX_FACEIDS];
  REF_INT nfaceid1, faceids1[REF_GEOM_MAX_FACEIDS];
  REF_INT *e2f;
  REF_INT degree, max_node = 50, *node_list;

  /* to allow recon after meshb load times */
  RSS(ref_geom_initialize(ref_geom), "clear out previous assoc");
  RSS(ref_cell_free(ref_grid_edg(ref_grid)), "clear out edge");
  RSS(ref_cell_create(&ref_grid_edg(ref_grid), 2, REF_TRUE), "edg");

  RSS(ref_geom_recon_nodes(ref_grid, &cad_nodes), "recover nodes");

  RSS(ref_geom_edge_faces(ref_grid, &e2f), "compute edge faces");

  ref_malloc(node_list, max_node, REF_INT);

  for (id = 1; id <= ref_geom->nedge; id++) {
    object = ((ego *)(ref_geom->edges))[id - 1];
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");
    if (mtype == DEGENERATE) {
      printf(" topo edge id %3d degen\n", id);
    } else {
      int toponode0, toponode1;
      REF_INT node0, node1;
      double closest[3];
      REF_INT next_node, current_node;
      REF_INT geom;
      REIS(TWONODE, mtype, "ONENODE edge not implemented");
      REIS(2, nchild, "expect two topo node for edge");
      toponode0 = EG_indexBodyTopo(ref_geom->solid, pchldrn[0]);
      toponode1 = EG_indexBodyTopo(ref_geom->solid, pchldrn[1]);
      node0 = cad_nodes[toponode0 - 1];
      node1 = cad_nodes[toponode1 - 1];
      printf(" topo edge id %3d fid %d %d trange [%f,%f]\n", id,
             e2f[0 + 2 * (id - 1)], e2f[1 + 2 * (id - 1)], trange[0],
             trange[1]);
      REIS(EGADS_SUCCESS, EG_evaluate(object, &(trange[0]), xyz), "EG eval");
      node = node0;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      printf("  node0 id %2d index %3d t %f dist %e fid", toponode0, node0,
             trange[0], dist);
      RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node0,
                                  REF_GEOM_MAX_FACEIDS, &nfaceid0, faceids0),
          REF_INCREASE_LIMIT, "count faceids");
      for (i = 0; i < nfaceid0; i++) printf(" %d", faceids0[i]);
      printf("\n");
      if (show_xyz) {
        printf(" d %23.15e %23.15e %23.15e\n", ref_node_xyz(ref_node, 0, node),
               ref_node_xyz(ref_node, 1, node),
               ref_node_xyz(ref_node, 2, node));
        printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      }
      REIS(EGADS_SUCCESS, EG_evaluate(object, &(trange[1]), xyz), "EG eval");
      node = node1;
      dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                  pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                  pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
      printf("  node1 id %2d index %3d t %f dist %e fid", toponode1, node1,
             trange[1], dist);
      RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), node1,
                                  REF_GEOM_MAX_FACEIDS, &nfaceid1, faceids1),
          REF_INCREASE_LIMIT, "count faceids");
      for (i = 0; i < nfaceid1; i++) printf(" %d", faceids1[i]);
      printf("\n");
      if (show_xyz) {
        printf(" d %23.15e %23.15e %23.15e\n", ref_node_xyz(ref_node, 0, node),
               ref_node_xyz(ref_node, 1, node),
               ref_node_xyz(ref_node, 2, node));
        printf(" c %23.15e %23.15e %23.15e\n", xyz[0], xyz[1], xyz[2]);
      }
      current_node = node0;
      t = trange[0];
      {
        best_node = node0;
        param[0] = trange[0];
        RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
                                  ref_node_xyz_ptr(ref_node, best_node), param),
            "inv wrapper");
        REIS(EGADS_SUCCESS, EG_evaluate(object, &(param[0]), closest),
             "EG eval");
        dist = sqrt(pow(closest[0] - ref_node_xyz(ref_node, 0, best_node), 2) +
                    pow(closest[1] - ref_node_xyz(ref_node, 1, best_node), 2) +
                    pow(closest[2] - ref_node_xyz(ref_node, 2, best_node), 2));
        printf("   best_node %5d t %f best_dist %e fid", best_node, param[0],
               dist);
        RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node,
                                    REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
            REF_INCREASE_LIMIT, "count faceids");
        for (i = 0; i < nfaceid; i++) printf(" %d", faceids[i]);
        printf("\n");
        if (param[0] < trange[0] || trange[1] < param[0])
          THROW("edge t out of range");
        if (show_xyz) {
          printf(" d %23.15e %23.15e %23.15e\n",
                 ref_node_xyz(ref_node, 0, best_node),
                 ref_node_xyz(ref_node, 1, best_node),
                 ref_node_xyz(ref_node, 2, best_node));
          printf(" c %23.15e %23.15e %23.15e\n", closest[0], closest[1],
                 closest[2]);
        }
      }
      while (current_node != node1) {
        RSS(ref_cell_node_list_around(ref_grid_tri(ref_grid), current_node,
                                      max_node, &degree, node_list),
            "next node");
        best_node = REF_EMPTY;
        best_dist = 1.0e20;
        best_param = 1.0e20;
        for (i = 0; i < degree; i++) {
          REF_INT fid;
          REF_BOOL have_faceid0, have_faceid1;
          next_node = node_list[i];
          if (REF_SUCCESS ==
              ref_geom_find(ref_geom, next_node, REF_GEOM_EDGE, id, &geom))
            continue; /* this canidate is already part of the edge */
          RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), next_node,
                                      REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
              REF_INCREASE_LIMIT, "count faceids");
          if (nfaceid < 2) continue; /* should be between faces */
          have_faceid0 = REF_FALSE;
          for (fid = 0; fid < nfaceid; fid++)
            if (e2f[0 + 2 * (id - 1)] == faceids[fid]) have_faceid0 = REF_TRUE;
          have_faceid1 = REF_FALSE;
          for (fid = 0; fid < nfaceid; fid++)
            if (e2f[1 + 2 * (id - 1)] == faceids[fid]) have_faceid1 = REF_TRUE;
          if (!have_faceid0 || !have_faceid1)
            continue; /* must have expected faceids */
          param[0] = t;
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_EDGE, id,
                                    ref_node_xyz_ptr(ref_node, next_node),
                                    param),
              "inv wrapper");
          REIS(EGADS_SUCCESS, EG_evaluate(object, &(param[0]), closest),
               "EG eval");
          dist =
              sqrt(pow(closest[0] - ref_node_xyz(ref_node, 0, next_node), 2) +
                   pow(closest[1] - ref_node_xyz(ref_node, 1, next_node), 2) +
                   pow(closest[2] - ref_node_xyz(ref_node, 2, next_node), 2));
          if (dist < best_dist) {
            best_node = next_node;
            best_dist = dist;
            best_param = param[0];
          }
        }
        if (REF_EMPTY == best_node) {
          RSS(REF_FAILURE, "count not find next node");
        }
        param[0] = best_param;
        printf("   best_node %5d t %f best_dist %e fid", best_node, param[0],
               best_dist);
        RXS(ref_cell_id_list_around(ref_grid_tri(ref_grid), best_node,
                                    REF_GEOM_MAX_FACEIDS, &nfaceid, faceids),
            REF_INCREASE_LIMIT, "count faceids");
        for (i = 0; i < nfaceid; i++) printf(" %d", faceids[i]);
        printf("\n");
        if (param[0] < trange[0] || trange[1] < param[0])
          THROW("edge t out of range");
        if (show_xyz) {
          printf(" d %23.15e %23.15e %23.15e\n",
                 ref_node_xyz(ref_node, 0, best_node),
                 ref_node_xyz(ref_node, 1, best_node),
                 ref_node_xyz(ref_node, 2, best_node));
          printf(" c %23.15e %23.15e %23.15e\n", closest[0], closest[1],
                 closest[2]);
        }
        param[1] = t;
        nodes[0] = best_node;
        nodes[1] = current_node;
        nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &cell), "add e");
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, id, &(param[0])),
            "add geom edge 0");
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_EDGE, id, &(param[1])),
            "add geom edge 1");
        current_node = best_node;
        t = param[0];
      }
    }
    if (debug) ref_geom_tec(ref_grid, "ref_geom_each_edge.tec");
  }
  ref_free(cad_nodes);
  ref_free(node_list);
  ref_free(e2f);

  printf("seting face UV under edges\n");
  each_ref_cell_valid_cell_with_nodes(ref_grid_edg(ref_grid), cell,
                                      edge_nodes) {
    REF_INT edgeid, faceid;
    int sense;
    REF_INT ncell, tri_list[2], tri_nodes[REF_CELL_MAX_SIZE_PER];
    edgeid = edge_nodes[2];
    RSS(ref_cell_list_with2(ref_grid_tri(ref_grid), edge_nodes[0],
                            edge_nodes[1], 2, &ncell, tri_list),
        "tri list for edge");
    REIS(2, ncell, "expect two tri for edge");
    for (i = 0; i < ncell; i++) {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri_list[i], tri_nodes),
          "nodes");
      faceid = tri_nodes[3];
      sense = 0; /* when edge is used twice in loop, this is its index */
      RSS(ref_geom_tuv(ref_geom, edge_nodes[0], REF_GEOM_EDGE, edgeid, &t),
          "edge t0");
      REIS(
          EGADS_SUCCESS,
          EG_getEdgeUV(((ego *)(ref_geom->faces))[faceid - 1],
                       ((ego *)(ref_geom->edges))[edgeid - 1], sense, t, param),
          "eval edge face uv");
      RSS(ref_geom_add(ref_geom, edge_nodes[0], REF_GEOM_FACE, faceid, param),
          "add geom face for edge");
      RSS(ref_geom_tuv(ref_geom, edge_nodes[1], REF_GEOM_EDGE, edgeid, &t),
          "edge t0");
      REIS(
          EGADS_SUCCESS,
          EG_getEdgeUV(((ego *)(ref_geom->faces))[faceid - 1],
                       ((ego *)(ref_geom->edges))[edgeid - 1], sense, t, param),
          "eval edge face uv");
      RSS(ref_geom_add(ref_geom, edge_nodes[1], REF_GEOM_FACE, faceid, param),
          "add geom face for edge");
    }
  }

  RSS(ref_geom_degen_param(ref_grid), "UV cleanup at singularities");

  pass = 0;
  updates = REF_EMPTY;
  while (0 != updates) {
    REF_INT geom0, geom1, geom2;
    REF_INT faceid;
    REF_GEOM_SURFACE surface_type;
    double closest[3];
    REF_BOOL inv_eval_wrapper = REF_TRUE;
    updates = 0;
    pass++;
    each_ref_cell_valid_cell_with_nodes(ref_grid_tri(ref_grid), cell, nodes) {
      faceid = nodes[3];
      RSS(ref_geom_face_surface(ref_geom, faceid, &surface_type), "styp");
      RXS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_FACE, faceid, &geom0),
          REF_NOT_FOUND, "find0");
      RXS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_FACE, faceid, &geom1),
          REF_NOT_FOUND, "find1");
      RXS(ref_geom_find(ref_geom, nodes[2], REF_GEOM_FACE, faceid, &geom2),
          REF_NOT_FOUND, "find2");
      if (REF_EMPTY != geom0 && REF_EMPTY == geom1) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[0], REF_GEOM_FACE, faceid, param),
            "geom0");
        if (inv_eval_wrapper && surface_type == REF_GEOM_PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[1]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[1]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_FACE, faceid, param),
            "add face");
      }
      if (REF_EMPTY != geom1 && REF_EMPTY == geom2) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[1], REF_GEOM_FACE, faceid, param),
            "geom1");
        if (inv_eval_wrapper && surface_type == REF_GEOM_PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[2]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[2]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[2], REF_GEOM_FACE, faceid, param),
            "add face");
      }
      if (REF_EMPTY != geom2 && REF_EMPTY == geom0) {
        updates++;
        RSS(ref_geom_tuv(ref_geom, nodes[2], REF_GEOM_FACE, faceid, param),
            "geom2");
        if (inv_eval_wrapper && surface_type == REF_GEOM_PLANE) {
          for (i = 0; i < 3; i++)
            closest[i] = ref_node_xyz(ref_node, i, nodes[0]);
          RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, closest,
                                    param),
              "inv eval");
        } else {
          REIS(EGADS_SUCCESS,
               EG_invEvaluate(((ego *)(ref_geom->faces))[faceid - 1],
                              ref_node_xyz_ptr(ref_node, nodes[0]), param,
                              closest),
               "EG eval");
        }
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_FACE, faceid, param),
            "add face");
      }
    }
    printf(" pass %3d updates %d\n", pass, updates);
  }
  return REF_SUCCESS;
#else
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (REF_EMPTY == ref_geom->nnode) printf("No EGADS loaded\n");
  printf("No EGADS linked for %s\n", __func__);
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_uv_area(REF_GEOM ref_geom, REF_INT *nodes,
                            REF_DBL *uv_area) {
  REF_INT id = nodes[3];
  REF_DBL uv0[2], uv1[2], uv2[2];
  RSS(ref_geom_tuv(ref_geom, nodes[0], REF_GEOM_FACE, id, uv0), "uv0");
  RSS(ref_geom_tuv(ref_geom, nodes[1], REF_GEOM_FACE, id, uv1), "uv1");
  RSS(ref_geom_tuv(ref_geom, nodes[2], REF_GEOM_FACE, id, uv2), "uv2");
  *uv_area = 0.5 * (-uv1[0] * uv0[1] + uv2[0] * uv0[1] + uv0[0] * uv1[1] -
                    uv2[0] * uv1[1] - uv0[0] * uv2[1] + uv1[0] * uv2[1]);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_sign(REF_GRID ref_grid, REF_INT id, REF_DBL *sign) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (NULL == ((ref_geom)->uv_area_sign)) {
    REF_CELL ref_cell = ref_grid_tri(ref_grid);
    REF_INT face;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL uv_area;
    if (REF_EMPTY == ref_geom->nface)
      RSS(ref_geom_infer_nedge_nface(ref_grid), "infer counts");
    ref_malloc_init(ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 0.0);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      face = nodes[3];
      if (face < 1 || ref_geom->nface < face) continue;
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
      if (uv_area < 0.0) {
        ((ref_geom)->uv_area_sign)[face - 1] -= 1.0;
      } else {
        ((ref_geom)->uv_area_sign)[face - 1] += 1.0;
      }
    }
    for (face = 0; face < ref_geom->nface; face++) {
      if (((ref_geom)->uv_area_sign)[face] < 0.0) {
        ((ref_geom)->uv_area_sign)[face] = -1.0;
      } else {
        ((ref_geom)->uv_area_sign)[face] = 1.0;
      }
    }
  }

  if (id < 1 || id > ref_geom->nface) return REF_INVALID;
  *sign = ((ref_geom)->uv_area_sign)[id - 1];

  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_report(REF_GRID ref_grid) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT geom, id, min_id, max_id;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL no_cell;
  REF_DBL uv_area, total_uv_area, min_uv_area, max_uv_area, sign_uv_area;
  REF_INT n_neg, n_pos;

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++) {
    no_cell = REF_TRUE;
    total_uv_area = 0.0;
    min_uv_area = 0.0;
    max_uv_area = 0.0;
    n_neg = 0;
    n_pos = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                        nodes) if (id == nodes[3]) {
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
      total_uv_area += uv_area;
      if (no_cell) {
        min_uv_area = uv_area;
        max_uv_area = uv_area;
        no_cell = REF_FALSE;
      } else {
        min_uv_area = MIN(min_uv_area, uv_area);
        max_uv_area = MAX(max_uv_area, uv_area);
      }
      if (uv_area < 0.0) {
        n_neg++;
      } else {
        n_pos++;
      }
    }
    if (!no_cell) {
      RSS(ref_geom_uv_area_sign(ref_grid, id, &sign_uv_area), "sign");
      printf("face%5d: %4.1f %9.2e total (%10.3e,%10.3e) %d + %d -\n", id,
             sign_uv_area, total_uv_area, min_uv_area, max_uv_area, n_pos,
             n_neg);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_inspect(REF_GEOM ref_geom) {
  REF_INT geom;
  printf("ref_geom = %p\n", (void *)ref_geom);
  printf(" n = %d, max = %d\n", ref_geom_n(ref_geom), ref_geom_max(ref_geom));
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++) {
    switch (ref_geom_type(ref_geom, geom)) {
      case REF_GEOM_NODE:
        printf("%d node: %d global, %d id\n", geom, ref_geom_id(ref_geom, geom),
               ref_geom_node(ref_geom, geom));
        break;
      case REF_GEOM_EDGE:
        printf("%d edge: %d id, %d global, t=%e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom));
        break;
      case REF_GEOM_FACE:
        printf("%d face: %d id, %d global, uv= %e %e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom),
               ref_geom_param(ref_geom, 1, geom));
        break;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tattle(REF_GEOM ref_geom, REF_INT node) {
  REF_INT item, geom;

  printf(" tattle on node = %d\n", node);
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    switch (ref_geom_type(ref_geom, geom)) {
      case REF_GEOM_NODE:
        printf("%d node: %d global, %d id\n", geom, ref_geom_id(ref_geom, geom),
               ref_geom_node(ref_geom, geom));
        break;
      case REF_GEOM_EDGE:
        printf("%d edge: %d id, %d global, t=%e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom));
        break;
      case REF_GEOM_FACE:
        printf("%d face: %d id, %d global, uv= %e %e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom),
               ref_geom_param(ref_geom, 1, geom));
        break;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_supported(REF_GEOM ref_geom, REF_INT node,
                              REF_BOOL *has_support) {
  *has_support = !ref_adj_empty(ref_geom_adj(ref_geom), node);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_add(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param) {
  REF_INT geom;
  REF_INT orig, chunk;
  REF_INT max_limit = REF_INT_MAX / 3;
  REF_STATUS status;

  if (type < 0 || 2 < type) return REF_INVALID;

  status = ref_geom_find(ref_geom, node, type, id, &geom);
  RXS(status, REF_NOT_FOUND, "find failed");

  if (REF_SUCCESS == status) {
    if (type > 0) ref_geom_param(ref_geom, 0, geom) = param[0];
    if (type > 1) ref_geom_param(ref_geom, 1, geom) = param[1];
    return REF_SUCCESS;
  }

  if (REF_EMPTY == ref_geom_blank(ref_geom)) {
    RAS(ref_geom_max(ref_geom) != max_limit,
        "the number of geoms is too large for integers, cannot grow");
    orig = ref_geom_max(ref_geom);
    /* geometric growth for efficiency */
    chunk = MAX(1000, (REF_INT)(1.5 * (REF_DBL)orig));

    /* try to keep under 32-bit limit */
    RAS(max_limit - orig > 0, "chunk limit at max");
    chunk = MIN(chunk, max_limit - orig);

    ref_geom_max(ref_geom) = orig + chunk;

    ref_realloc(ref_geom->descr, 3 * ref_geom_max(ref_geom), REF_INT);
    ref_realloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);

    for (geom = orig; geom < ref_geom_max(ref_geom); geom++) {
      ref_geom_type(ref_geom, geom) = REF_EMPTY;
      ref_geom_id(ref_geom, geom) = geom + 1;
    }
    ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
    ref_geom_blank(ref_geom) = orig;
  }

  geom = ref_geom_blank(ref_geom);
  ref_geom_blank(ref_geom) = ref_geom_id(ref_geom, geom);

  ref_geom_type(ref_geom, geom) = type;
  ref_geom_id(ref_geom, geom) = id;
  ref_geom_node(ref_geom, geom) = node;

  ref_geom_param(ref_geom, 0, geom) = 0.0;
  ref_geom_param(ref_geom, 1, geom) = 0.0;
  if (type > 0) ref_geom_param(ref_geom, 0, geom) = param[0];
  if (type > 1) ref_geom_param(ref_geom, 1, geom) = param[1];

  RSS(ref_adj_add(ref_geom->ref_adj, node, geom), "register geom");

  ref_geom_n(ref_geom)++;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_remove(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                           REF_INT id) {
  REF_INT geom;
  REF_STATUS status;

  status = ref_geom_find(ref_geom, node, type, id, &geom);
  RXS(status, REF_NOT_FOUND, "find failed");

  if (REF_SUCCESS == status) {
    RSS(ref_adj_remove(ref_geom_adj(ref_geom), ref_geom_node(ref_geom, geom),
                       geom),
        "unregister geom");

    ref_geom_type(ref_geom, geom) = REF_EMPTY;
    ref_geom_id(ref_geom, geom) = ref_geom_blank(ref_geom);
    ref_geom_blank(ref_geom) = geom;
    ref_geom_n(ref_geom)--;
  }

  return status;
}

REF_STATUS ref_geom_remove_all(REF_GEOM ref_geom, REF_INT node) {
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;

  item = ref_adj_first(ref_adj, node);
  while (ref_adj_valid(item)) {
    geom = ref_adj_item_ref(ref_adj, item);
    RSS(ref_adj_remove(ref_adj, node, geom), "unregister geom");

    ref_geom_type(ref_geom, geom) = REF_EMPTY;
    ref_geom_id(ref_geom, geom) = ref_geom_blank(ref_geom);
    ref_geom_blank(ref_geom) = geom;
    ref_geom_n(ref_geom)--;

    item = ref_adj_first(ref_adj, node);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_is_a(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_BOOL *it_is) {
  REF_INT item, geom;
  *it_is = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      *it_is = REF_TRUE;
      return REF_SUCCESS;
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_unique_id(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                              REF_INT *id) {
  REF_INT item, geom;
  REF_BOOL found_one;
  found_one = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      if (found_one) return REF_INVALID; /* second one makes invalid */
      found_one = REF_TRUE;
      *id = ref_geom_id(ref_geom, geom);
    }
  }
  if (found_one) return REF_SUCCESS;
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_find(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_INT id, REF_INT *found) {
  REF_INT item, geom;
  *found = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom) &&
        id == ref_geom_id(ref_geom, geom)) {
      *found = geom;
      return REF_SUCCESS;
    }
  }
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_tuv(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param) {
  REF_INT geom;

  RSS(ref_geom_find(ref_geom, node, type, id, &geom), "not found");

  if (type > 0) param[0] = ref_geom_param(ref_geom, 0, geom);
  if (type > 1) param[1] = ref_geom_param(ref_geom, 1, geom);

  return REF_SUCCESS;
}

static REF_STATUS ref_geom_eval_edge_face_uv(REF_GEOM ref_geom,
                                             REF_INT edge_geom) {
#ifdef HAVE_EGADS
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT node, item, face_geom;
  double t;
  double uv[2];
  int sense = 0;
  ego *edges, *faces;
  ego edge, face;
  int status;
  REF_INT faceid;
  REF_DBL xyz[3];

  if (edge_geom < 0 || ref_geom_max(ref_geom) <= edge_geom) return REF_INVALID;
  if (REF_GEOM_EDGE != ref_geom_type(ref_geom, edge_geom)) return REF_INVALID;

  edges = (ego *)(ref_geom->edges);
  edge = edges[ref_geom_id(ref_geom, edge_geom) - 1];

  t = ref_geom_param(ref_geom, 0, edge_geom);

  node = ref_geom_node(ref_geom, edge_geom);

  faces = (ego *)(ref_geom->faces);
  each_ref_adj_node_item_with_ref(ref_adj, node, item, face_geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
      faceid = ref_geom_id(ref_geom, face_geom);
      face = faces[faceid - 1];
      status = EG_getEdgeUV(face, edge, sense, t, uv);
      if (EGADS_TOPOERR == status) {
        RSS(ref_geom_eval(ref_geom, edge_geom, xyz, NULL), "eval edge");
        uv[0] = 0.0;
        uv[1] = 0.0;
        RSS(ref_geom_inverse_eval(ref_geom, REF_GEOM_FACE, faceid, xyz, uv),
            "inv eval");
        ref_geom_param(ref_geom, 0, face_geom) = uv[0];
        ref_geom_param(ref_geom, 1, face_geom) = uv[1];
      } else {
        REIS(EGADS_SUCCESS, status, "eval edge face uv");
        ref_geom_param(ref_geom, 0, face_geom) = uv[0];
        ref_geom_param(ref_geom, 1, face_geom) = uv[1];
      }
    }
  }

  return REF_SUCCESS;
#else
  if (edge_geom < 0 || ref_geom_max(ref_geom) <= edge_geom) return REF_INVALID;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_xyz_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_DBL *xyz) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT item0, item1;
  REF_INT geom0, geom1;
  REF_INT type, id;
  REF_DBL param[2], param0[2], param1[2];
  REF_BOOL has_id;
  REF_INT i;

  for (i = 0; i < 3; i++)
    xyz[i] = 0.5 * (ref_node_xyz(ref_node, i, node0) +
                    ref_node_xyz(ref_node, i, node1));

  /* edge between */
  type = REF_GEOM_EDGE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0) {
    each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                    geom1) {
      if (ref_geom_type(ref_geom, geom0) == type &&
          ref_geom_type(ref_geom, geom1) == type &&
          ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
        id = ref_geom_id(ref_geom, geom0);
        RSS(ref_cell_side_has_id(ref_grid_edg(ref_grid), node0, node1, id,
                                 &has_id),
            "has edge id");
        if (has_id) {
          RSS(ref_geom_tuv(ref_geom, node0, type, id, param0), "node0");
          RSS(ref_geom_tuv(ref_geom, node1, type, id, param1), "node1");
          param[0] = 0.5 * (param0[0] + param1[0]);
          if (ref_geom_model_loaded(ref_geom))
            RSB(ref_geom_inverse_eval(ref_geom, type, id, xyz, param),
                "inv eval edge",
                ref_geom_tec(ref_grid, "ref_geom_split_edge.tec"));
          return REF_SUCCESS;
        }
      }
    }
  }

  /* insert face between */
  type = REF_GEOM_FACE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0) {
    each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                    geom1) {
      if (ref_geom_type(ref_geom, geom0) == type &&
          ref_geom_type(ref_geom, geom1) == type &&
          ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
        id = ref_geom_id(ref_geom, geom0);
        RSS(ref_cell_side_has_id(ref_grid_tri(ref_grid), node0, node1, id,
                                 &has_id),
            "has edge id");
        if (has_id) {
          RSS(ref_geom_tuv(ref_geom, node0, type, id, param0), "node0");
          RSS(ref_geom_tuv(ref_geom, node1, type, id, param1), "node1");
          param[0] = 0.5 * (param0[0] + param1[0]);
          param[1] = 0.5 * (param0[1] + param1[1]);
          if (ref_geom_model_loaded(ref_geom)) {
            RSB(ref_geom_inverse_eval(ref_geom, type, id, xyz, param),
                "inv eval face",
                ref_geom_tec(ref_grid, "ref_geom_split_face.tec"));
          }
          return REF_SUCCESS;
        }
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_add_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_INT new_node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT item0, item1;
  REF_INT geom0, geom1;
  REF_INT type, id;
  REF_DBL param[2], param0[2], param1[2];
  REF_DBL uv_min[2], uv_max[2];
  REF_BOOL has_id;
  REF_BOOL has_edge_support;
  REF_INT edge_geom;

  /* insert edge between */
  has_edge_support = REF_FALSE;
  edge_geom = REF_EMPTY;
  type = REF_GEOM_EDGE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0) {
    each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                    geom1) {
      if (ref_geom_type(ref_geom, geom0) == type &&
          ref_geom_type(ref_geom, geom1) == type &&
          ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
        id = ref_geom_id(ref_geom, geom0);
        RSS(ref_cell_side_has_id(ref_grid_edg(ref_grid), node0, node1, id,
                                 &has_id),
            "has edge id");
        if (has_id) {
          if (REF_EMPTY != edge_geom) { /* should only be one */
            printf("previous edge geom %d id %d\n", edge_geom,
                   ref_geom_id(ref_geom, edge_geom));
            printf("this edge geom %d %d id %d\n", geom0, geom1,
                   ref_geom_id(ref_geom, geom0));
            RSS(ref_geom_tec(ref_grid, "ref_geom_edge_edge.tec"), "tec");
          }
          REIS(REF_EMPTY, edge_geom, "found a second edge spanning nodes");
          RAS(!has_edge_support, "already has support? second edge?");
          RSS(ref_geom_tuv(ref_geom, node0, type, id, param0), "node0");
          RSS(ref_geom_tuv(ref_geom, node1, type, id, param1), "node1");
          param[0] = 0.5 * (param0[0] + param1[0]);
          if (ref_geom_model_loaded(ref_geom))
            RSB(ref_geom_inverse_eval(ref_geom, type, id,
                                      ref_node_xyz_ptr(ref_node, new_node),
                                      param),
                "inv eval edge",
                ref_geom_tec(ref_grid, "ref_geom_split_edge.tec"));
          /* enforce bounding box and use midpoint as full-back */
          if (param[0] < MIN(param0[0], param1[0]) ||
              MAX(param0[0], param1[0]) < param[0])
            param[0] = 0.5 * (param0[0] + param1[0]);
          RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
          has_edge_support = REF_TRUE;
          RSS(ref_geom_find(ref_geom, new_node, type, id, &edge_geom),
              "find the new edge for later face uv evaluation");
        }
      }
    }
  }

  /* insert face between */
  type = REF_GEOM_FACE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0) {
    each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                    geom1) {
      if (ref_geom_type(ref_geom, geom0) == type &&
          ref_geom_type(ref_geom, geom1) == type &&
          ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
        id = ref_geom_id(ref_geom, geom0);
        RSS(ref_cell_side_has_id(ref_grid_tri(ref_grid), node0, node1, id,
                                 &has_id),
            "has edge id");
        if (has_id) {
          RSS(ref_geom_tuv(ref_geom, node0, type, id, param0), "node0");
          RSS(ref_geom_tuv(ref_geom, node1, type, id, param1), "node1");
          param[0] = 0.5 * (param0[0] + param1[0]);
          param[1] = 0.5 * (param0[1] + param1[1]);
          RSS(ref_geom_is_a(ref_geom, new_node, REF_GEOM_EDGE,
                            &has_edge_support),
              "look for edge");
          if (ref_geom_model_loaded(ref_geom) && !has_edge_support) {
            RSB(ref_geom_inverse_eval(ref_geom, type, id,
                                      ref_node_xyz_ptr(ref_node, new_node),
                                      param),
                "inv eval face",
                ref_geom_tec(ref_grid, "ref_geom_split_face.tec"));
            /* enforce bounding box of node0 and try midpoint */
            RSS(ref_geom_tri_uv_bounding_box2(ref_grid, node0, node1, id,
                                              uv_min, uv_max),
                "bb");
            if (param[0] < uv_min[0] || uv_max[0] < param[0] ||
                param[1] < uv_min[1] || uv_max[1] < param[1]) {
              param[0] = 0.5 * (param0[0] + param1[0]);
              param[1] = 0.5 * (param0[1] + param1[1]);
            }
          }
          RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
        }
      }
    }
  }

  /* if there is an edge between, set the face uv based on edge t */
  if (ref_geom_model_loaded(ref_geom) && has_edge_support) {
    RSS(ref_geom_eval_edge_face_uv(ref_geom, edge_geom),
        "set face uv under edge");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_support_between(REF_GRID ref_grid, REF_INT node0,
                                    REF_INT node1, REF_BOOL *needs_support) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT item0, item1;
  REF_INT geom0, geom1;
  REF_INT type, id;
  REF_BOOL has_id;

  *needs_support = REF_FALSE;
  /* assume face check is sufficient */
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0)
      each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                      geom1) {
    type = REF_GEOM_FACE;
    if (ref_geom_type(ref_geom, geom0) == type &&
        ref_geom_type(ref_geom, geom1) == type &&
        ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
      id = ref_geom_id(ref_geom, geom0);
      RSS(ref_cell_side_has_id(ref_grid_tri(ref_grid), node0, node1, id,
                               &has_id),
          "has edge id");
      if (has_id) {
        *needs_support = REF_TRUE;
        return REF_SUCCESS;
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_uv_bounding_box(REF_GRID ref_grid, REF_INT node,
                                        REF_DBL *uv_min, REF_DBL *uv_max) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node, id, iuv;
  REF_DBL uv[2];

  /* get face id and initialize min and max */
  RSS(ref_geom_unique_id(ref_geom, node, REF_GEOM_FACE, &id), "id");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_min), "uv_min");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_max), "uv_max");

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_geom_tuv(ref_geom, ref_cell_c2n(ref_cell, cell_node, cell),
                       REF_GEOM_FACE, id, uv),
          "uv");
      for (iuv = 0; iuv < 2; iuv++) uv_min[iuv] = MIN(uv_min[iuv], uv[iuv]);
      for (iuv = 0; iuv < 2; iuv++) uv_max[iuv] = MAX(uv_max[iuv], uv[iuv]);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_uv_bounding_box2(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_INT id,
                                         REF_DBL *uv_min, REF_DBL *uv_max) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, cell_node, iuv;
  REF_DBL uv[2];
  REF_INT i, ncell, cells[2];

  /* get face id and initialize min and max */
  RSS(ref_geom_tuv(ref_geom, node0, REF_GEOM_FACE, id, uv_min), "uv_min");
  RSS(ref_geom_tuv(ref_geom, node0, REF_GEOM_FACE, id, uv_max), "uv_max");

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells),
      "get list");

  for (i = 0; i < ncell; i++) {
    cell = cells[i];
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_geom_tuv(ref_geom, ref_cell_c2n(ref_cell, cell_node, cell),
                       REF_GEOM_FACE, id, uv),
          "uv");
      for (iuv = 0; iuv < 2; iuv++) uv_min[iuv] = MIN(uv_min[iuv], uv[iuv]);
      for (iuv = 0; iuv < 2; iuv++) uv_max[iuv] = MAX(uv_max[iuv], uv[iuv]);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_constrain(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;
  REF_BOOL have_geom_node;
  REF_BOOL have_geom_edge;
  REF_BOOL have_geom_face;
  REF_INT node_geom;
  REF_INT edge_geom;
  REF_INT face_geom;
  REF_DBL xyz[3];

  /* put sloppy geom handling here */

  /* no geom, do nothing */
  if (ref_adj_empty(ref_adj, node)) return REF_SUCCESS;

  have_geom_node = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_NODE == ref_geom_type(ref_geom, geom)) {
      have_geom_node = REF_TRUE;
      node_geom = geom;
      break;
    }
  }

  if (have_geom_node) { /* update T of edges? update UV of (degen) faces? */
    RSS(ref_geom_eval(ref_geom, node_geom, xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, node_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    return REF_SUCCESS;
  }

  have_geom_edge = REF_FALSE;
  edge_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
      have_geom_edge = REF_TRUE;
      edge_geom = geom;
      break;
    }
  }

  /* edge geom, evaluate edge and update face uv */
  if (have_geom_edge) {
    RSS(ref_geom_eval(ref_geom, edge_geom, xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    RSS(ref_geom_eval_edge_face_uv(ref_geom, edge_geom), "resol edge uv");
    return REF_SUCCESS;
  }

  /* look for face geom */
  have_geom_face = REF_FALSE;
  face_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
      have_geom_face = REF_TRUE;
      face_geom = geom;
      break;
    }
  }

  /* face geom, evaluate on face uv */
  if (have_geom_face) {
    RSS(ref_geom_eval(ref_geom, face_geom, xyz, NULL), "eval face");
    node = ref_geom_node(ref_geom, face_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}
/*
   [x_t,y_t,z_t] edge
   [x_tt,y_tt,z_tt]
   [x_u,y_u,z_u] [x_v,y_v,z_v] face
   [x_uu,y_uu,z_uu] [x_uv,y_uv,z_uv] [x_vv,y_vv,z_vv]
*/
REF_STATUS ref_geom_eval(REF_GEOM ref_geom, REF_INT geom, REF_DBL *xyz,
                         REF_DBL *dxyz_dtuv) {
#ifdef HAVE_EGADS
  double eval[18];
  double params[2];
  REF_INT i;
  ego *nodes, *edges, *faces;
  ego object;
  int status;
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  params[0] = 0.0;
  params[1] = 0.0;
  object = (ego)NULL;
  switch (ref_geom_type(ref_geom, geom)) {
    case (REF_GEOM_NODE):
      RNS(ref_geom->nodes, "nodes not loaded");
      if (ref_geom_id(ref_geom, geom) < 1 ||
          ref_geom_id(ref_geom, geom) > ref_geom->nnode)
        return REF_INVALID;
      nodes = (ego *)(ref_geom->nodes);
      object = nodes[ref_geom_id(ref_geom, geom) - 1];
      {
        ego ref, *pchldrn;
        int oclass, mtype, nchild, *psens;
        REIS(EGADS_SUCCESS,
             EG_getTopology(object, &ref, &oclass, &mtype, xyz, &nchild,
                            &pchldrn, &psens),
             "EG topo node");
      }
      return REF_SUCCESS;
      break;
    case (REF_GEOM_EDGE):
      RNS(ref_geom->edges, "edges not loaded");
      if (ref_geom_id(ref_geom, geom) < 1 ||
          ref_geom_id(ref_geom, geom) > ref_geom->nedge)
        return REF_INVALID;
      edges = (ego *)(ref_geom->edges);
      object = edges[ref_geom_id(ref_geom, geom) - 1];
      params[0] = ref_geom_param(ref_geom, 0, geom);
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      if (ref_geom_id(ref_geom, geom) < 1 ||
          ref_geom_id(ref_geom, geom) > ref_geom->nface)
        return REF_INVALID;
      faces = (ego *)(ref_geom->faces);
      object = faces[ref_geom_id(ref_geom, geom) - 1];
      params[0] = ref_geom_param(ref_geom, 0, geom);
      params[1] = ref_geom_param(ref_geom, 1, geom);
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom");
  }

  status = EG_evaluate(object, params, eval);
  if (EGADS_SUCCESS != status) {
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double trange[2];
    printf("geom %d type %d id %d\n", geom, ref_geom_type(ref_geom, geom),
           ref_geom_id(ref_geom, geom));
    if (ref_geom_type(ref_geom, geom) > 0) printf("param[0] = %f\n", params[0]);
    if (ref_geom_type(ref_geom, geom) > 1) printf("param[1] = %f\n", params[1]);
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, trange, &nchild,
                        &pchldrn, &psens),
         "EG topo node");
    printf("trange %f %f\n", trange[0], trange[1]);
    REIS(EGADS_SUCCESS, status, "eval");
  }
  xyz[0] = eval[0];
  xyz[1] = eval[1];
  xyz[2] = eval[2];
  if (NULL != dxyz_dtuv) {
    for (i = 0; i < 6; i++) dxyz_dtuv[i] = eval[3 + i];
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom))
      for (i = 0; i < 9; i++) dxyz_dtuv[6 + i] = eval[9 + i];
  }
  return REF_SUCCESS;
#else
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  printf("evaluating to (0,0,0), No EGADS linked for %s\n", __func__);
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  if (NULL != dxyz_dtuv) dxyz_dtuv[0] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_inverse_eval(REF_GEOM ref_geom, REF_INT type, REF_INT id,
                                 REF_DBL *xyz, REF_DBL *param) {
#ifdef HAVE_EGADS
  ego object = (ego)NULL;
  int i, guess_status, noguess_status;
  REF_DBL guess_param[2], noguess_param[2];
  REF_DBL guess_closest[3], noguess_closest[3];
  REF_DBL guess_dist, noguess_dist;

  ego ref, *pchldrn;
  int oclass, mtype, nchild, *psens;
  double range[4];

  REF_BOOL guess_in_range, noguess_in_range;

  REF_BOOL verbose = REF_FALSE;

  switch (type) {
    case (REF_GEOM_NODE):
      printf("GEOM_NODE ref_geom_inverse_eval not defined\n");
      return REF_IMPLEMENT;
      break;
    case (REF_GEOM_EDGE):
      RNS(ref_geom->edges, "edges not loaded");
      if (id < 1 || id > ref_geom->nedge) return REF_INVALID;
      object = ((ego *)(ref_geom->edges))[id - 1];
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      if (id < 1 || id > ref_geom->nface) return REF_INVALID;
      object = ((ego *)(ref_geom->faces))[id - 1];
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom");
  }

  REIS(EGADS_SUCCESS,
       EG_getTopology(object, &ref, &oclass, &mtype, range, &nchild, &pchldrn,
                      &psens),
       "EG topo node");

  for (i = 0; i < type; i++) guess_param[i] = param[i];
  guess_status = EG_invEvaluateGuess(object, xyz, guess_param, guess_closest);
  noguess_status = EG_invEvaluate(object, xyz, noguess_param, noguess_closest);
  if (verbose) printf("guess %d noguess %d\n", guess_status, noguess_status);
  if (EGADS_SUCCESS != guess_status && EGADS_SUCCESS != noguess_status)
    return REF_FAILURE;

  if (EGADS_SUCCESS == guess_status && EGADS_SUCCESS != noguess_status) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
    return REF_SUCCESS;
  }

  if (EGADS_SUCCESS != guess_status && EGADS_SUCCESS == noguess_status) {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
    return REF_SUCCESS;
  }

  guess_dist = sqrt(pow(guess_closest[0] - xyz[0], 2) +
                    pow(guess_closest[1] - xyz[1], 2) +
                    pow(guess_closest[2] - xyz[2], 2));
  noguess_dist = sqrt(pow(noguess_closest[0] - xyz[0], 2) +
                      pow(noguess_closest[1] - xyz[1], 2) +
                      pow(noguess_closest[2] - xyz[2], 2));

  guess_in_range = REF_TRUE;
  for (i = 0; i < type; i++)
    if (guess_param[i] < range[0 + 2 * i] || range[1 + 2 * i] < guess_param[i])
      guess_in_range = REF_FALSE;

  noguess_in_range = REF_TRUE;
  for (i = 0; i < type; i++)
    if (noguess_param[i] < range[0 + 2 * i] ||
        range[1 + 2 * i] < noguess_param[i])
      noguess_in_range = REF_FALSE;

  if (verbose) {
    printf("guess %e noguess %e\n", guess_dist, noguess_dist);
    for (i = 0; i < type; i++)
      printf("%d:guess %f noguess %f\n", i, guess_param[i], noguess_param[i]);
    printf("guess %d noguess %d\n", guess_in_range, noguess_in_range);
  }

  if (guess_in_range && !noguess_in_range) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
    return REF_SUCCESS;
  }

  if (!guess_in_range && noguess_in_range) {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
    return REF_SUCCESS;
  }

  if (guess_dist < noguess_dist) {
    for (i = 0; i < type; i++) param[i] = guess_param[i];
  } else {
    for (i = 0; i < type; i++) param[i] = noguess_param[i];
  }
  return REF_SUCCESS;

#else
  printf("no-op, No EGADS linked for %s type %d id %d x %f p %f n %d\n",
         __func__, type, id, xyz[0], param[0], ref_geom_n(ref_geom));
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_curvature(REF_GEOM ref_geom, REF_INT geom, REF_DBL *kr,
                              REF_DBL *r, REF_DBL *ks, REF_DBL *s) {
#ifdef HAVE_EGADS
  double curvature[8];
  double params[2];
  ego *faces;
  ego object;
  int egads_status;
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  params[0] = 0.0;
  params[1] = 0.0;
  object = (ego)NULL;
  switch (ref_geom_type(ref_geom, geom)) {
    case (REF_GEOM_NODE):
      RSS(REF_IMPLEMENT, "geom node");
      break;
    case (REF_GEOM_EDGE):
      RSS(REF_IMPLEMENT, "geom edge");
      break;
    case (REF_GEOM_FACE):
      RNS(ref_geom->faces, "faces not loaded");
      faces = (ego *)(ref_geom->faces);
      object = faces[ref_geom_id(ref_geom, geom) - 1];
      params[0] = ref_geom_param(ref_geom, 0, geom);
      params[1] = ref_geom_param(ref_geom, 1, geom);
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom");
  }
  RNS(object, "EGADS object is NULL. Has the geometry been loaded?");

  egads_status = EG_curvature(object, params, curvature);
  if (EGADS_DEGEN == egads_status) {
    REF_DBL xyz[3], dxyz_duv[15], du, dv;
    ego ref, *pchldrn;
    int oclass, mtype, nchild, *psens;
    double uv_range[4];
    REF_DBL shift = 1.0e-2;
    RSS(ref_geom_eval(ref_geom, geom, xyz, dxyz_duv), "eval");
    du = sqrt(ref_math_dot(&(dxyz_duv[0]), &(dxyz_duv[0])));
    dv = sqrt(ref_math_dot(&(dxyz_duv[3]), &(dxyz_duv[3])));
    REIS(EGADS_SUCCESS,
         EG_getTopology(object, &ref, &oclass, &mtype, uv_range, &nchild,
                        &pchldrn, &psens),
         "EG topo face");
    if (du > dv) {
      params[0] =
          (1.0 - shift) * params[0] + shift * 0.5 * (uv_range[0] + uv_range[1]);
    } else {
      params[1] =
          (1.0 - shift) * params[1] + shift * 0.5 * (uv_range[2] + uv_range[3]);
    }
    egads_status = EG_curvature(object, params, curvature);
  }
  REIS(EGADS_SUCCESS, egads_status, "curve");
  *kr = curvature[0];
  r[0] = curvature[1];
  r[1] = curvature[2];
  r[2] = curvature[3];
  *ks = curvature[4];
  s[0] = curvature[5];
  s[1] = curvature[6];
  s[2] = curvature[7];
  return REF_SUCCESS;
#else
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  printf("curvature 0, 0: No EGADS linked for %s\n", __func__);
  *kr = 0.0;
  r[0] = 1.0;
  r[1] = 0.0;
  r[2] = 0.0;
  *ks = 0.0;
  s[0] = 0.0;
  s[1] = 1.0;
  s[2] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_uv_rsn(REF_DBL *uv, REF_DBL *r, REF_DBL *s, REF_DBL *n,
                           REF_DBL *drsduv) {
  REF_INT i;
  REF_DBL dot;
  REF_DBL len;

  for (i = 0; i < 3; i++) r[i] = uv[i];
  drsduv[0] = 1.0;
  drsduv[1] = 0.0;
  for (i = 0; i < 3; i++) s[i] = uv[i + 3];
  drsduv[2] = 0.0;
  drsduv[3] = 1.0;
  len = sqrt(ref_math_dot(r, r));
  drsduv[0] /= len;
  drsduv[1] /= len;
  RSS(ref_math_normalize(r), "norm r (u)");
  len = sqrt(ref_math_dot(s, s));
  drsduv[2] /= len;
  drsduv[3] /= len;
  RSS(ref_math_normalize(s), "norm s (v)");

  dot = ref_math_dot(r, s);
  for (i = 0; i < 3; i++) s[i] -= dot * r[i];
  drsduv[2] -= dot * drsduv[0];
  drsduv[3] -= dot * drsduv[1];

  len = sqrt(ref_math_dot(s, s));
  drsduv[2] /= len;
  drsduv[3] /= len;
  RSS(ref_math_normalize(s), "norm s (v)");

  ref_math_cross_product(r, s, n);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_rsn(REF_GEOM ref_geom, REF_INT geom, REF_DBL *r, REF_DBL *s,
                        REF_DBL *n) {
  REF_DBL kr, ks;
  RSS(ref_geom_curvature(ref_geom, geom, &kr, r, &ks, s), "eval face");
  ref_math_cross_product(r, s, n);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_param(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom;
  REF_INT node;
  REF_DBL xyz[3];
  REF_DBL dist, max, max_node, max_edge, global_max;
  REF_BOOL node_constraint, edge_constraint;

  max = 0.0;
  each_ref_geom_node(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_geom_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    max = MAX(max, dist);
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max node");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo node max eval dist %e\n", max);

  max = 0.0;
  max_node = 0.0;
  each_ref_geom_edge(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_geom_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &node_constraint), "n");
    if (node_constraint) {
      max_node = MAX(max_node, dist);
    } else {
      max = MAX(max, dist);
    }
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max edge");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo edge max eval dist %e\n", max);
  RSS(ref_mpi_max(ref_mpi, &max_node, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_node = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo edge node tol %e\n", max_node);

  max = 0.0;
  max_node = 0.0;
  max_edge = 0.0;
  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_geom_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &node_constraint), "n");
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &edge_constraint), "n");
    if (node_constraint) {
      max_node = MAX(max_node, dist);
    } else {
      if (edge_constraint) {
        max_edge = MAX(max_edge, dist);
      } else {
        max = MAX(max, dist);
      }
    }
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max face");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face max eval dist %e\n", max);
  RSS(ref_mpi_max(ref_mpi, &max_edge, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_edge = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face edge tol %e\n", max_edge);
  RSS(ref_mpi_max(ref_mpi, &max_node, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_node = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face node tol %e\n", max_node);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_topo(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT item, geom;
  REF_BOOL geom_node, geom_edge, geom_face;
  REF_BOOL no_face, no_edge;
  REF_BOOL found_one;
  REF_BOOL found_too_many;
  REF_INT cell, ncell, cell_list[2];

  for (node = 0; node < ref_node_max(ref_node); node++) {
    if (ref_node_valid(ref_node, node)) {
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node");
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge");
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &geom_face), "face");
      no_face = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_qua(ref_grid), node);
      no_edge = ref_cell_node_empty(ref_grid_edg(ref_grid), node);
      if (geom_node) {
        if (no_edge && ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) {
          THROW("geom node missing edge");
        }
        if (no_face && ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) {
          THROW("geom node missing tri or qua");
        }
      }
      if (geom_edge) {
        if (no_edge && ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) {
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("geom edge missing edge");
        }
        if (no_face && ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) {
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("geom edge missing tri or qua");
        }
      }
      if (geom_face) {
        if (no_face && ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) {
          printf("no face for geom\n");
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("geom face missing tri or qua");
        }
      }
      if (!no_face) {
        if (!geom_face) {
          printf("no geom for face\n");
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("geom face missing tri or qua");
        }
      }
      if (geom_edge && !geom_node) {
        found_one = REF_FALSE;
        found_too_many = REF_FALSE;
        each_ref_geom_having_node(ref_geom, node, item, geom) {
          if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
            if (found_one) found_too_many = REF_TRUE;
            found_one = REF_TRUE;
          }
        }
        if (!found_one || found_too_many) {
          if (!found_one) printf("none found\n");
          if (found_too_many) printf("found too many\n");
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("multiple geom edge away from geom node");
        }
      }
      if (geom_face && !geom_edge) {
        found_one = REF_FALSE;
        found_too_many = REF_FALSE;
        each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item,
                                        geom) {
          if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
            if (found_one) found_too_many = REF_TRUE;
            found_one = REF_TRUE;
          }
        }
        if (!found_one || found_too_many) {
          if (!found_one) printf("none found\n");
          if (found_too_many) printf("found too many\n");
          RSS(ref_node_location(ref_node, node), "loc");
          RSS(ref_geom_tattle(ref_geom, node), "tatt");
          RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
          THROW("multiple geom face away from geom edge");
        }
      }
    } else {
      if (!ref_adj_empty(ref_geom_adj(ref_geom), node))
        THROW("invalid node has geom");
    }
  }

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_list_with2(ref_cell, ref_cell_c2n(ref_cell, 0, cell),
                            ref_cell_c2n(ref_cell, 1, cell), 2, &ncell,
                            cell_list),
        "edge list for edge");
    if (2 == ncell) {
      printf("error: two edg found with same nodes\n");
      printf("edg %d n %d %d id %d\n", cell_list[0],
             ref_cell_c2n(ref_cell, 0, cell_list[0]),
             ref_cell_c2n(ref_cell, 1, cell_list[0]),
             ref_cell_c2n(ref_cell, 2, cell_list[0]));
      printf("edg %d n %d %d id %d\n", cell_list[1],
             ref_cell_c2n(ref_cell, 0, cell_list[1]),
             ref_cell_c2n(ref_cell, 1, cell_list[1]),
             ref_cell_c2n(ref_cell, 2, cell_list[1]));
      RSS(ref_node_location(ref_node, ref_cell_c2n(ref_cell, 0, cell)), "loc");
      RSS(ref_node_location(ref_node, ref_cell_c2n(ref_cell, 1, cell)), "loc");
      RSS(ref_geom_tattle(ref_geom, ref_cell_c2n(ref_cell, 0, cell)), "tatt");
      RSS(ref_geom_tattle(ref_geom, ref_cell_c2n(ref_cell, 1, cell)), "tatt");
      RSS(ref_geom_tec(ref_grid, "ref_geom_typo_error.tec"), "geom tec");
    }
    REIS(1, ncell, "expect only one edge cell for two nodes");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tetgen_volume(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  char *ugrid_name = "ref_geom_test.ugrid";
  char *smesh_name = "ref_geom_test.smesh";
  char *node_name = "ref_geom_test.1.node";
  char *ele_name = "ref_geom_test.1.ele";
  char command[1024];
  FILE *file;
  REF_INT nnode, ndim, attr, mark;
  REF_INT ntet, node_per;
  REF_INT node, nnode_surface, item, new_node;
  REF_DBL xyz[3], dist;
  REF_INT cell, new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  int system_status;

  printf("%d surface nodes %d triangles\n", ref_node_n(ref_node),
         ref_cell_n(ref_grid_tri(ref_grid)));

  printf("tec360 ref_geom_test_tetgen_geom.tec\n");
  RSS(ref_geom_tec(ref_grid, "ref_geom_test_tetgen_geom.tec"), "dbg geom");
  printf("tec360 ref_geom_test_tetgen_surf.tec\n");
  RSS(ref_export_tec_surf(ref_grid, "ref_geom_test_tetgen_surf.tec"),
      "dbg surf");
  RSS(ref_export_by_extension(ref_grid, ugrid_name), "ugrid");
  RSS(ref_export_smesh(ref_grid, smesh_name), "smesh");
  sprintf(command, "tetgen -pYq1.0/0z %s > %s.out", smesh_name, smesh_name);
  printf("%s\n", command);
  system_status = system(command);
  REIS(0, system_status, "tetgen failed");

  file = fopen(node_name, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", node_name);
  RNS(file, "unable to open file");

  REIS(1, fscanf(file, "%d", &nnode), "node header nnode");
  REIS(1, fscanf(file, "%d", &ndim), "node header ndim");
  REIS(3, ndim, "not 3D");
  REIS(1, fscanf(file, "%d", &attr), "node header attr");
  REIS(0, attr, "nodes have attribute 3D");
  REIS(1, fscanf(file, "%d", &mark), "node header mark");
  REIS(0, mark, "nodes have mark");

  /* verify surface nodes */
  nnode_surface = ref_node_n(ref_node);
  for (node = 0; node < nnode_surface; node++) {
    REIS(1, fscanf(file, "%d", &item), "node item");
    RES(node, item, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    dist = sqrt((xyz[0] - ref_node_xyz(ref_node, 0, node)) *
                    (xyz[0] - ref_node_xyz(ref_node, 0, node)) +
                (xyz[1] - ref_node_xyz(ref_node, 1, node)) *
                    (xyz[1] - ref_node_xyz(ref_node, 1, node)) +
                (xyz[2] - ref_node_xyz(ref_node, 2, node)) *
                    (xyz[2] - ref_node_xyz(ref_node, 2, node)));
    if (dist > 1.0e-12) {
      printf("node %d off by %e\n", node, dist);
      THROW("tetgen moved node");
    }
  }

  /* interior nodes */
  for (node = nnode_surface; node < nnode; node++) {
    REIS(1, fscanf(file, "%d", &item), "node item");
    RES(node, item, "node index");
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    ref_node_xyz(ref_node, 0, new_node) = xyz[0];
    ref_node_xyz(ref_node, 1, new_node) = xyz[1];
    ref_node_xyz(ref_node, 2, new_node) = xyz[2];
  }

  fclose(file);

  /* check faces when paranoid, but tetgen -z should not mess with them */

  file = fopen(ele_name, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", ele_name);
  RNS(file, "unable to open file");

  REIS(1, fscanf(file, "%d", &ntet), "ele header ntet");
  REIS(1, fscanf(file, "%d", &node_per), "ele header node_per");
  REIS(4, node_per, "expected tets");
  REIS(1, fscanf(file, "%d", &mark), "ele header mark");
  REIS(0, mark, "ele have mark");

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    REIS(1, fscanf(file, "%d", &item), "tet item");
    RES(cell, item, "node index");
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_import_ugrid_tets(REF_GRID ref_grid,
                                        const char *filename) {
  REF_CELL ref_cell;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_DBL xyz[3];
  REF_INT new_node, orig_nnode, node, tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua;
  REF_INT face_id;
  REF_INT cell, new_cell;

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fscanf(file, "%d", &nnode), "nnode");
  RES(1, fscanf(file, "%d", &ntri), "ntri");
  RES(1, fscanf(file, "%d", &nqua), "nqua");
  RES(1, fscanf(file, "%d", &ntet), "ntet");
  RES(1, fscanf(file, "%d", &npyr), "npyr");
  RES(1, fscanf(file, "%d", &npri), "npri");
  RES(1, fscanf(file, "%d", &nhex), "nhex");

  orig_nnode = ref_node_n(ref_node);

  for (node = 0; node < nnode; node++) {
    REIS(1, fscanf(file, "%lf", &(xyz[0])), "x");
    REIS(1, fscanf(file, "%lf", &(xyz[1])), "y");
    REIS(1, fscanf(file, "%lf", &(xyz[2])), "z");
    if (node >= orig_nnode) {
      RSS(ref_node_add(ref_node, node, &new_node), "new_node");
      REIS(node, new_node, "node index");
      ref_node_xyz(ref_node, 0, new_node) = xyz[0];
      ref_node_xyz(ref_node, 1, new_node) = xyz[1];
      ref_node_xyz(ref_node, 2, new_node) = xyz[2];
    }
  }

  for (tri = 0; tri < ntri; tri++) {
    for (node = 0; node < 3; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
  }
  for (qua = 0; qua < nqua; qua++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "qua");
  }

  for (tri = 0; tri < ntri; tri++) {
    RES(1, fscanf(file, "%d", &face_id), "tri id");
  }

  for (qua = 0; qua < nqua; qua++) {
    RES(1, fscanf(file, "%d", &face_id), "qua id");
  }

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_aflr_volume(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  char *surface_ugrid_name = "ref_geom_test_surface.ugrid";
  char *volume_ugrid_name = "ref_geom_test_volume.ugrid";
  char command[1024];
  int system_status;

  printf("%d surface nodes %d triangles\n", ref_node_n(ref_node),
         ref_cell_n(ref_grid_tri(ref_grid)));

  printf("tec360 ref_geom_test_aflr_geom.tec\n");
  RSS(ref_geom_tec(ref_grid, "ref_geom_test_aflr_geom.tec"), "dbg geom");
  printf("tec360 ref_geom_test_aflr_surf.tec\n");
  RSS(ref_export_tec_surf(ref_grid, "ref_geom_test_aflr_surf.tec"), "dbg surf");
  RSS(ref_export_by_extension(ref_grid, surface_ugrid_name), "ugrid");
  sprintf(command, "aflr3 -igrid %s -ogrid %s -mrecrbf=0 -angqbf=180 > %s.out",
          surface_ugrid_name, volume_ugrid_name, volume_ugrid_name);
  printf("%s\n", command);
  system_status = system(command);
  REIS(0, system_status, "aflr failed");

  RSS(ref_import_ugrid_tets(ref_grid, volume_ugrid_name), "tets only");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_load(REF_GEOM ref_geom, const char *filename) {
#ifdef HAVE_EGADS
  ego context;
  ego model = NULL;
  ego geom, *bodies, *children;
  int oclass, mtype, nbody, *senses, nchild;
  ego solid, *faces, *edges, *nodes;
  int nface, nedge, nnode;

  context = (ego)(ref_geom->context);

#ifdef HAVE_EGADS_LITE
  {
    /* entry point NOT in egads.h */
    int EG_importModel(egObject * context, const size_t nbytes,
                       const char stream[], egObject **model);

    SUPRESS_UNUSED_COMPILER_WARNING(filename);

    RAS(0 < ref_geom_cad_data_size(ref_geom), "zero size cad_data");
    RNS(ref_geom_cad_data(ref_geom), "cad_data NULL");
    REIS(EGADS_SUCCESS,
         EG_importModel(context, (size_t)ref_geom_cad_data_size(ref_geom),
                        (char *)ref_geom_cad_data(ref_geom), &model),
         "EG load");
  }
#else
  REIS(EGADS_SUCCESS, EG_loadModel(context, 0, filename, &model), "EG load");

  {
    /* entry point NOT in egads.h */
    int EG_exportModel(ego mobject, size_t * nbytes, char *stream[]);

    size_t cad_data_size;
    REF_BYTE *cad_data;

    REIS(EGADS_SUCCESS, EG_exportModel(model, &cad_data_size, &cad_data),
         "EG stream");
    ref_geom_cad_data_size(ref_geom) = cad_data_size;
    /* safe non-NULL free, if already allocated, to prevent memory leaks */
    ref_free(ref_geom->cad_data);
    ref_malloc(ref_geom_cad_data(ref_geom), ref_geom_cad_data_size(ref_geom),
               REF_BYTE);
    memcpy(ref_geom_cad_data(ref_geom), cad_data,
           ref_geom_cad_data_size(ref_geom));
    EG_free(cad_data);
  }
#endif

  REIS(EGADS_SUCCESS,
       EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody, &bodies,
                      &senses),
       "EG topo bodies");
  REIS(1, nbody, "expected 1 body");
  solid = bodies[0];
  REIS(EGADS_SUCCESS,
       EG_getTopology(solid, &geom, &oclass, &mtype, NULL, &nchild, &children,
                      &senses),
       "EG topo body type");
  REIS(SOLIDBODY, mtype, "expected SOLIDBODY");
  ref_geom->solid = (void *)solid;

  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, NODE, &nnode, &nodes),
       "EG node topo");
  ref_geom->nnode = nnode;
  ref_geom->nodes = (void *)nodes;
  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, EDGE, &nedge, &edges),
       "EG edge topo");
  ref_geom->nedge = nedge;
  ref_geom->edges = (void *)edges;
  REIS(EGADS_SUCCESS, EG_getBodyTopos(solid, NULL, FACE, &nface, &faces),
       "EG face topo");
  ref_geom->nface = nface;
  ref_geom->faces = (void *)faces;

#else
  printf("returning empty grid from %s, No EGADS linked for %s\n", __func__,
         filename);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_geom_infer_nedge_nface(REF_GRID ref_grid) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT min_id, max_id;
  RSS(ref_geom_faceid_range(ref_grid, &min_id, &max_id), "face range");
  REIS(1, min_id, "first face id not 1");
  ref_geom->nface = max_id;
  RSS(ref_geom_edgeid_range(ref_grid, &min_id, &max_id), "edge range");
  REIS(1, min_id, "first edge id not 1");
  ref_geom->nedge = max_id;
  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_diagonal(REF_GEOM ref_geom, REF_DBL *diag) {
#ifdef HAVE_EGADS
  ego solid;
  double box[6];
  solid = (ego)(ref_geom->solid);

  RNS(solid, "EGADS solid object is NULL. Has the geometry been loaded?");

  REIS(EGADS_SUCCESS, EG_getBoundingBox(solid, box), "EG bounding box");
  *diag = sqrt((box[0] - box[3]) * (box[0] - box[3]) +
               (box[1] - box[4]) * (box[1] - box[4]) +
               (box[2] - box[5]) * (box[2] - box[5]));

#else
  printf("returning 1.0 from %s, No EGADS\n", __func__);
  *diag = 1.0;
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_tess(REF_GRID ref_grid, REF_DBL *params) {
#ifdef HAVE_EGADS
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri, new_cell;
  REF_DBL param[2];
  double size;
  ego geom;
  ego solid, tess;
  int tess_status, nvert;
  int face, edge, plen, tlen;
  const double *points, *uv, *t;
  const int *ptype, *pindex, *tris, *tric;
  int node, new_node, pty, pin;
  double verts[3];

  solid = (ego)(ref_geom->solid);

  RSS(ref_geom_egads_diagonal(ref_geom, &size), "bbox diag");
  printf("suggested params %f %f %f, size %e\n", 0.25 * size, 0.001 * size,
         15.0, size);
  printf("   actual params %f %f %f\n", params[0], params[1], params[2]);
  /* maximum length of an EDGE segment or triangle side (in physical space) */
  /* curvature-based value that looks locally at the deviation between
     the centroid of the discrete object and the underlying geometry */
  /* maximum interior dihedral angle (in degrees) */

  /* unlimited override to make coarse init grids
     params[0] = 1.0;
     params[1] = 1.0;
     params[2] = 180.0;
  */

  REIS(EGADS_SUCCESS, EG_makeTessBody(solid, params, &tess), "EG tess");
  REIS(EGADS_SUCCESS, EG_statusTessBody(tess, &geom, &tess_status, &nvert),
       "EG tess");
  REIS(1, tess_status, "tess not closed");

  for (node = 0; node < nvert; node++) {
    REIS(EGADS_SUCCESS, EG_getGlobal(tess, node + 1, &pty, &pin, verts),
         "global node info");
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
    ref_node_xyz(ref_node, 0, node) = verts[0];
    ref_node_xyz(ref_node, 1, node) = verts[1];
    ref_node_xyz(ref_node, 2, node) = verts[2];
    if (0 == pty)
      RSS(ref_geom_add(ref_geom, node, REF_GEOM_NODE, pin, NULL), "node");
  }

  for (face = 0; face < (ref_geom->nface); face++) {
    REIS(EGADS_SUCCESS,
         EG_getTessFace(tess, face + 1, &plen, &points, &uv, &ptype, &pindex,
                        &tlen, &tris, &tric),
         "tess query face");

    for (node = 0; node < plen; node++) {
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, node + 1, &(nodes[0])), "l2g0");
      nodes[0] -= 1;
      param[0] = uv[0 + 2 * node];
      param[1] = uv[1 + 2 * node];
      RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_FACE, face + 1, param),
          "face uv");
    }
    for (tri = 0; tri < tlen; tri++) {
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[0 + 3 * tri], &(nodes[0])),
           "l2g0");
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[1 + 3 * tri], &(nodes[1])),
           "l2g1");
      REIS(EGADS_SUCCESS,
           EG_localToGlobal(tess, face + 1, tris[2 + 3 * tri], &(nodes[2])),
           "l2g2");
      nodes[0] -= 1;
      nodes[1] -= 1;
      nodes[2] -= 1;
      nodes[3] = face + 1;
      RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell), "new tri");
    }
  }

  RSS(ref_geom_degen_param(ref_grid), "UV cleanup at singularities");

  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    int egads_status;
    REF_BOOL degenerate;
    degenerate = REF_FALSE;
    REIS(EGADS_SUCCESS, EG_getTessEdge(tess, edge + 1, &plen, &points, &t),
         "tess query edge");
    for (node = 0; node < plen; node++) {
      egads_status = EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0]));
      if (EGADS_DEGEN == egads_status) {
        degenerate = REF_TRUE;
      } else {
        REIS(EGADS_SUCCESS, egads_status, "l2g0");
        nodes[0] -= 1;
        param[0] = t[node];
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, edge + 1, param),
            "edge t");
      }
    }
    if (!degenerate)
      for (node = 0; node < (plen - 1); node++) {
        /* assue edge index is 1-bias */
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 1, &(nodes[0])),
             "l2g0");
        REIS(EGADS_SUCCESS,
             EG_localToGlobal(tess, -(edge + 1), node + 2, &(nodes[1])),
             "l2g1");
        nodes[0] -= 1;
        nodes[1] -= 1;
        nodes[2] = edge + 1;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell), "new edge");
      }
  }

#else
  printf("returning empty grid from %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  SUPRESS_UNUSED_COMPILER_WARNING(params);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_geom_degen_param(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT face;

  for (face = 0; face < (ref_geom->nface); face++) {
    ego esurf, *eloops, eref;
    int oclass, mtype, nloop, *senses, *pinfo;
    double trange[4], data[18], *preal;
    REIS(EGADS_SUCCESS,
         EG_getTopology(((ego *)(ref_geom->faces))[face], &esurf, &oclass,
                        &mtype, data, &nloop, &eloops, &senses),
         "topo");
    REIS(EGADS_SUCCESS,
         EG_getGeometry(esurf, &oclass, &mtype, &eref, &pinfo, &preal), "geom");
    EG_free(pinfo);
    EG_free(preal);
    if (mtype != PLANE) {
      ego ecurve, *eedges, *echilds;
      int iloop, iedge, nedge, nchild, inode;
      REF_INT node;
      REF_DBL param[2];
      /* loop through all Loops associated with this Face */
      for (iloop = 0; iloop < nloop; iloop++) {
        /* loop through all Edges associated with this Loop */
        REIS(EGADS_SUCCESS,
             EG_getTopology(eloops[iloop], &ecurve, &oclass, &mtype, data,
                            &nedge, &eedges, &senses),
             "topo");
        for (iedge = 0; iedge < nedge; iedge++) {
          REIS(EGADS_SUCCESS,
               EG_getTopology(eedges[iedge], &ecurve, &oclass, &mtype, trange,
                              &nchild, &echilds, &senses),
               "tp");
          if (mtype == DEGENERATE) {
            double uvmin[6], uvmax[6];
            REF_INT geom;
            printf("face id %d has degen\n", face + 1);
            /* find index of bounding Node */
            inode = EG_indexBodyTopo((ego)(ref_geom->solid), echilds[0]);
            REIS(EGADS_SUCCESS,
                 EG_evaluate(eedges[iedge + nedge], &trange[0], uvmin),
                 "eval min");
            REIS(EGADS_SUCCESS,
                 EG_evaluate(eedges[iedge + nedge], &trange[1], uvmax),
                 "eval max");
            each_ref_geom_node(ref_geom,
                               geom) if (inode == ref_geom_id(ref_geom, geom)) {
              node = ref_geom_node(ref_geom, geom);
              printf("tess node index %d\n", node);
              ref_node_location(ref_grid_node(ref_grid), node);
              RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, face + 1, param),
                  "face uv");
              printf("u tess %f tmin %f tmax %f\n", param[0], uvmin[0],
                     uvmax[0]);
              printf("v tess %f tmin %f tmax %f\n", param[1], uvmin[1],
                     uvmax[1]);
              param[0] = 0.5 * (uvmin[0] + uvmax[0]);
              param[1] = 0.5 * (uvmin[1] + uvmax[1]);
              RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, face + 1, param),
                  "face uv");
            }
          }
        }
      }
    }
  }
#else
  printf("unable to %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_tess_repair_topo(REF_GRID ref_grid) {
#ifdef HAVE_EGADS
  REF_CELL ref_edg = ref_grid_edg(ref_grid);
  REF_CELL ref_tri = ref_grid_tri(ref_grid);
  REF_INT edg, nedg, edg_list[2];
  REF_INT ntri, tri_list[2];
  REF_INT node0, node1;
  REF_INT *e2f;

  RSS(ref_geom_edge_faces(ref_grid, &e2f), "edge2face");
  
  each_ref_cell_valid_cell(ref_edg, edg) {
    node0 = ref_cell_c2n(ref_edg, 0, edg);
    node1 = ref_cell_c2n(ref_edg, 1, edg);
    RSS(ref_cell_list_with2(ref_edg, node0, node1, 2, &nedg, edg_list),
        "edge list for nodes");
    if (2 == nedg) {
      printf("two edg found with same nodes %d %d from %d\n",node0,node1,edg);
      printf("edg %d n %d %d id %d\n", edg_list[0],
             ref_cell_c2n(ref_edg, 0, edg_list[0]),
             ref_cell_c2n(ref_edg, 1, edg_list[0]),
             ref_cell_c2n(ref_edg, 2, edg_list[0]));
      printf("edg %d n %d %d id %d\n", edg_list[1],
             ref_cell_c2n(ref_edg, 0, edg_list[1]),
             ref_cell_c2n(ref_edg, 1, edg_list[1]),
             ref_cell_c2n(ref_edg, 2, edg_list[1]));
      RSS(ref_cell_list_with2(ref_tri, node0, node1, 2, &ntri,
                              tri_list),
          "tri list for nodes");
    }
  }

  ref_free(e2f);
  
#else
  printf("unable to %s, No EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_geom_edge_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nedg;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_edge(ref_geom, geom) if (id == ref_geom_id(ref_geom, geom))
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
  nnode = ref_dict_n(ref_dict);

  nedg = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) if (id == nodes[2])
      nedg++;

  if (0 == nnode || 0 == nedg) /* skip degenerate */
  {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"edge%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, nedg, "point", "felineseg");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), 0.0, 0.0,
            ref_geom_param(ref_geom, 0, geom));
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                      nodes) if (id == nodes[2]) {
    for (node = 0; node < 2; node++) {
      RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
      fprintf(file, " %d", local + 1);
    }
    fprintf(file, "\n");
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_face_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_face(ref_geom, geom) if (id == ref_geom_id(ref_geom, geom))
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
  nnode = ref_dict_n(ref_dict);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) if (id == nodes[3])
      ntri++;

  if (0 == nnode || 0 == ntri) /* skip degenerate */
  {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"face%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, ntri, "point", "fetriangle");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), ref_geom_param(ref_geom, 0, geom),
            ref_geom_param(ref_geom, 1, geom), 0.0);
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                      nodes) if (id == nodes[3]) {
    for (node = 0; node < 3; node++) {
      RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
      fprintf(file, " %d", local + 1);
    }
    fprintf(file, "\n");
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_norm_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL r[3], s[3], n[3];
  REF_DBL area_sign;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_face(ref_geom, geom) if (id == ref_geom_id(ref_geom, geom))
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
  nnode = ref_dict_n(ref_dict);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) if (id == nodes[3])
      ntri++;

  if (0 == nnode || 0 == ntri) /* skip degenerate */
  {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"norm%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, ntri, "point", "fetriangle");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom), "not found");
    RSS(ref_geom_rsn(ref_geom, geom, r, s, n), "rsn");
    RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign), "a sign");
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), -area_sign * n[0],
            -area_sign * n[1], -area_sign * n[2]);
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                      nodes) if (id == nodes[3]) {
    for (node = 0; node < 3; node++) {
      RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
      fprintf(file, " %d", local + 1);
    }
    fprintf(file, "\n");
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_curve_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL kr, r[3], ks, s[3];

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_face(ref_geom, geom) if (id == ref_geom_id(ref_geom, geom))
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
  nnode = ref_dict_n(ref_dict);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) if (id == nodes[3])
      ntri++;

  if (0 == nnode || 0 == ntri) /* skip degenerate */
  {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(file,
          "zone t=\"curve%d\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          id, nnode, ntri, "point", "fetriangle");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom), "not found");
    RSS(ref_geom_curvature(ref_geom, geom, &kr, r, &ks, s), "curve");
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), ABS(kr), ABS(ks), 0.0);
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                      nodes) if (id == nodes[3]) {
    for (node = 0; node < 3; node++) {
      RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
      fprintf(file, " %d", local + 1);
    }
    fprintf(file, "\n");
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tec(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *file;
  REF_INT geom, id, min_id, max_id;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine cad coupling in tecplot format\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"t\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_edge(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_geom_edge_tec_zone(ref_grid, id, file), "tec edge");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_geom_face_tec_zone(ref_grid, id, file), "tec face");

  if (REF_FALSE && ref_geom_model_loaded(ref_geom)) {
    for (id = min_id; id <= max_id; id++)
      RSS(ref_geom_norm_tec_zone(ref_grid, id, file), "tec norm");
    for (id = min_id; id <= max_id; id++)
      RSS(ref_geom_curve_tec_zone(ref_grid, id, file), "tec curve");
  }

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_ghost(REF_GEOM ref_geom, REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_nnode, *b_nnode;
  REF_INT a_nnode_total, b_nnode_total;
  REF_INT *a_global, *b_global;
  REF_INT *a_part, *b_part;
  REF_INT *a_ngeom, *b_ngeom;
  REF_INT a_ngeom_total, b_ngeom_total;
  REF_INT *a_tgi, *b_tgi;
  REF_DBL *a_param, *b_param;
  REF_INT part, node, degree;
  REF_INT *a_next, *b_next;
  REF_INT local, item, geom;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_next, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_next, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(a_nnode, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_nnode, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(a_ngeom, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_ngeom, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) if (ref_mpi_rank(ref_mpi) !=
                                               ref_node_part(ref_node, node))
      a_nnode[ref_node_part(ref_node, node)]++;

  RSS(ref_mpi_alltoall(ref_mpi, a_nnode, b_nnode, REF_INT_TYPE),
      "alltoall nnodes");

  a_nnode_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_nnode_total += a_nnode[part];
  ref_malloc(a_global, a_nnode_total, REF_INT);
  ref_malloc(a_part, a_nnode_total, REF_INT);

  b_nnode_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_nnode_total += b_nnode[part];
  ref_malloc(b_global, b_nnode_total, REF_INT);
  ref_malloc(b_part, b_nnode_total, REF_INT);

  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_nnode[part - 1];

  each_ref_node_valid_node(ref_node, node) if (ref_mpi_rank(ref_mpi) !=
                                               ref_node_part(ref_node, node)) {
    part = ref_node_part(ref_node, node);
    a_global[a_next[part]] = ref_node_global(ref_node, node);
    a_part[a_next[part]] = ref_mpi_rank(ref_mpi);
    a_next[ref_node_part(ref_node, node)]++;
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_nnode, b_global, b_nnode, 1,
                        REF_INT_TYPE),
      "alltoallv global");
  RSS(ref_mpi_alltoallv(ref_mpi, a_part, a_nnode, b_part, b_nnode, 1,
                        REF_INT_TYPE),
      "alltoallv global");

  for (node = 0; node < b_nnode_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
    part = b_part[node];
    RSS(ref_adj_degree(ref_geom_adj(ref_geom), local, &degree), "deg");
    /* printf("%d: node %d global %d local %d part %d degree %d\n",
       ref_mpi_rank(ref_mpi), node,b_global[node], local, part, degree); */
    b_ngeom[part] += degree;
  }

  RSS(ref_mpi_alltoall(ref_mpi, b_ngeom, a_ngeom, REF_INT_TYPE),
      "alltoall ngeoms");

  a_ngeom_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_ngeom_total += a_ngeom[part];
  ref_malloc(a_tgi, 3 * a_ngeom_total, REF_INT);
  ref_malloc(a_param, 2 * a_ngeom_total, REF_DBL);

  b_ngeom_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_ngeom_total += b_ngeom[part];
  ref_malloc(b_tgi, 3 * b_ngeom_total, REF_INT);
  ref_malloc(b_param, 2 * b_ngeom_total, REF_DBL);

  b_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) b_next[part] =
      b_next[part - 1] + b_ngeom[part - 1];

  for (node = 0; node < b_nnode_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
    part = b_part[node];
    each_ref_geom_having_node(ref_geom, local, item, geom) {
      b_tgi[0 + 3 * b_next[part]] = ref_geom_type(ref_geom, geom);
      b_tgi[1 + 3 * b_next[part]] =
          ref_node_global(ref_node, ref_geom_node(ref_geom, geom));
      b_tgi[2 + 3 * b_next[part]] = ref_geom_id(ref_geom, geom);
      b_param[0 + 2 * b_next[part]] = ref_geom_param(ref_geom, 0, geom);
      b_param[1 + 2 * b_next[part]] = ref_geom_param(ref_geom, 1, geom);
      b_next[part]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_tgi, b_ngeom, a_tgi, a_ngeom, 3,
                        REF_INT_TYPE),
      "alltoallv tgi");
  RSS(ref_mpi_alltoallv(ref_mpi, b_param, b_ngeom, a_param, a_ngeom, 2,
                        REF_DBL_TYPE),
      "alltoallv param");

  for (geom = 0; geom < a_ngeom_total; geom++) {
    RSS(ref_node_local(ref_node, a_tgi[1 + 3 * geom], &local), "g2l");
    RSS(ref_geom_add(ref_geom, local, a_tgi[0 + 3 * geom], a_tgi[2 + 3 * geom],
                     &(a_param[2 * geom])),
        "add ghost");
  }

  free(b_param);
  free(b_tgi);
  free(a_param);
  free(a_tgi);
  free(b_part);
  free(b_global);
  free(a_part);
  free(a_global);
  free(b_ngeom);
  free(a_ngeom);
  free(b_nnode);
  free(a_nnode);
  free(b_next);
  free(a_next);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_faceid_range(REF_GRID ref_grid, REF_INT *min_faceid,
                                 REF_INT *max_faceid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_faceid = REF_INT_MAX;
  *max_faceid = REF_INT_MIN;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_node_per(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_node_per(ref_cell)]);
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_node_per(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_node_per(ref_cell)]);
  }

  if (ref_mpi_para(ref_mpi)) {
    REF_INT global;

    RSS(ref_mpi_min(ref_mpi, min_faceid, &global, REF_INT_TYPE),
        "mpi min face");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi min face");
    *min_faceid = global;

    RSS(ref_mpi_max(ref_mpi, max_faceid, &global, REF_INT_TYPE),
        "mpi max face");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi max face");
    *max_faceid = global;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_edgeid_range(REF_GRID ref_grid, REF_INT *min_edgeid,
                                 REF_INT *max_edgeid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_edgeid = REF_INT_MAX;
  *max_edgeid = REF_INT_MIN;

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_edgeid = MIN(*min_edgeid, nodes[ref_cell_node_per(ref_cell)]);
    *max_edgeid = MAX(*max_edgeid, nodes[ref_cell_node_per(ref_cell)]);
  }

  if (ref_mpi_para(ref_mpi)) {
    REF_INT global;

    RSS(ref_mpi_min(ref_mpi, min_edgeid, &global, REF_INT_TYPE),
        "mpi min edge");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi min edge");
    *min_edgeid = global;

    RSS(ref_mpi_max(ref_mpi, max_edgeid, &global, REF_INT_TYPE),
        "mpi max edge");
    RSS(ref_mpi_bcast(ref_mpi, &global, 1, REF_INT_TYPE), "mpi max edge");
    *max_edgeid = global;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_face_match(REF_GRID ref_grid) {
#if defined(HAVE_EGADS) && !defined(HAVE_EGADS_LITE)
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  double *cad_box;
  double *face_box;
  double *cad_cga;
  double *face_cga;
  double massprop[14];
  REF_DBL area, centroid[3];
  ego face_ego;
  REF_INT face, faceid, min_faceid, max_faceid, nfaceid;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT i, j;
  REF_DBL *face_norm;
  REF_DBL norm;
  REF_INT *candidates;
  REF_DICT ref_dict;
  REF_INT old_faceid, new_faceid;

  ref_malloc(cad_box, 6 * ref_geom->nface, double);
  ref_malloc(cad_cga, 4 * ref_geom->nface, double);
  for (face = 0; face < (ref_geom->nface); face++) {
    face_ego = ((ego *)(ref_geom->faces))[face];
    REIS(EGADS_SUCCESS, EG_getBoundingBox(face_ego, &(cad_box[6 * face])),
         "EG bounding box");
    REIS(EGADS_SUCCESS, EG_getMassProperties(face_ego, massprop),
         "EG mass properties");
    cad_cga[0 + 4 * face] = massprop[2];
    cad_cga[1 + 4 * face] = massprop[3];
    cad_cga[2 + 4 * face] = massprop[4];
    cad_cga[3 + 4 * face] = massprop[1];
  }

  RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "id range");
  nfaceid = max_faceid - min_faceid + 1;
  ref_malloc(face_box, 6 * nfaceid, double);
  ref_malloc(face_cga, 4 * nfaceid, double);
  for (face = 0; face < nfaceid; face++) {
    faceid = face + min_faceid;
    for (j = 0; j < 3; j++) {
      face_box[j + 0 * 3 + 6 * face] = 1.0e200;
      face_box[j + 1 * 3 + 6 * face] = -1.0e200;
    }
    for (j = 0; j < 4; j++) {
      face_cga[j + 4 * face] = 0.0;
    }
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (faceid == nodes[ref_cell_node_per(ref_cell)]) {
        for (i = 0; i < ref_cell_node_per(ref_cell); i++) {
          for (j = 0; j < 3; j++) {
            face_box[j + 0 * 3 + 6 * face] =
                MIN(face_box[j + 0 * 3 + 6 * face],
                    ref_node_xyz(ref_node, j, nodes[i]));
            face_box[j + 1 * 3 + 6 * face] =
                MAX(face_box[j + 1 * 3 + 6 * face],
                    ref_node_xyz(ref_node, j, nodes[i]));
          }
        }
        RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
        RSS(ref_node_tri_centroid(ref_node, nodes, centroid), "cent");
        for (j = 0; j < 3; j++) {
          face_cga[j + 4 * face] += centroid[j] * area;
        }
        face_cga[3 + 4 * face] += area;
      }
    }
    for (j = 0; j < 3; j++) {
      face_cga[j + 4 * face] /= face_cga[3 + 4 * face];
    }
  }

  for (face = 0; face < (ref_geom->nface); face++) {
    printf("%4d cmin %10.6f %10.6f %10.6f\n", face + 1, cad_box[0 + 6 * face],
           cad_box[1 + 6 * face], cad_box[2 + 6 * face]);
    printf("%4d cmax %10.6f %10.6f %10.6f\n", face + 1, cad_box[3 + 6 * face],
           cad_box[4 + 6 * face], cad_box[5 + 6 * face]);
    printf("%4d cmpp %10.6f %10.6f %10.6f %10.6f\n", face + 1,
           cad_cga[0 + 4 * face], cad_cga[1 + 4 * face], cad_cga[2 + 4 * face],
           cad_cga[3 + 4 * face]);
  }

  for (face = 0; face < nfaceid; face++) {
    faceid = face + min_faceid;
    printf("%4d fmin %10.6f %10.6f %10.6f\n", faceid, face_box[0 + 6 * face],
           face_box[1 + 6 * face], face_box[2 + 6 * face]);
    printf("%4d fmax %10.6f %10.6f %10.6f\n", faceid, face_box[3 + 6 * face],
           face_box[4 + 6 * face], face_box[5 + 6 * face]);
    printf("%4d fmpp %10.6f %10.6f %10.6f %10.6f\n", face + 1,
           face_cga[0 + 4 * face], face_cga[1 + 4 * face],
           face_cga[2 + 4 * face], face_cga[3 + 4 * face]);
  }

  ref_malloc(face_norm, nfaceid, REF_DBL);
  ref_malloc(candidates, nfaceid, REF_INT);
  RSS(ref_dict_create(&ref_dict), "create dict");
  for (i = 0; i < (ref_geom->nface); i++) {
    /* bbox */
    for (face = 0; face < nfaceid; face++) {
      norm = 0.0;
      for (j = 0; j < 6; j++) {
        norm += pow(face_box[j + 6 * face] - cad_box[j + 6 * i], 2);
      }
      norm = sqrt(norm);
      face_norm[face] = norm;
    }
    RSS(ref_sort_heap_dbl(nfaceid, face_norm, candidates), "sort");
    printf("%4d bbox %4.2f", i + 1,
           face_norm[candidates[0]] / face_norm[candidates[1]]);
    if (0.5 < face_norm[candidates[0]] / face_norm[candidates[1]]) {
      printf(" *");
    } else {
      if (candidates[0] + min_faceid == i + 1) {
        printf(" m");
      } else {
        printf("  ");
      }
    }
    for (j = 0; j < 3; j++) {
      printf(" %4d %10.3e", min_faceid + candidates[j],
             face_norm[candidates[j]]);
    }
    printf("\n");
    /* centroid */
    for (face = 0; face < nfaceid; face++) {
      norm = 0.0;
      for (j = 0; j < 3; j++) {
        norm += pow(face_cga[j + 4 * face] - cad_cga[j + 4 * i], 2);
      }
      norm = sqrt(norm);
      face_norm[face] = norm;
    }
    RSS(ref_sort_heap_dbl(nfaceid, face_norm, candidates), "sort");
    printf("%4d cent %4.2f", i + 1,
           face_norm[candidates[0]] / face_norm[candidates[1]]);
    if (0.5 < face_norm[candidates[0]] / face_norm[candidates[1]]) {
      printf(" *");
    } else {
      if (candidates[0] + min_faceid == i + 1) {
        printf(" m");
      } else {
        printf("  ");
      }
    }
    for (j = 0; j < 3; j++) {
      printf(" %4d %10.3e", min_faceid + candidates[j],
             face_norm[candidates[j]]);
    }
    printf("\n");

    /* last best canidate */
    RSS(ref_dict_store(ref_dict, candidates[0] + min_faceid, i + 1),
        "best guess at new number");
  }

  RSS(ref_dict_inspect(ref_dict), "list map");
  printf("replacing %d faces\n", ref_dict_n(ref_dict));
  each_ref_cell_valid_cell(ref_cell, cell) {
    old_faceid = ref_cell_c2n(ref_cell, ref_cell_node_per(ref_cell), cell);
    RSS(ref_dict_value(ref_dict, old_faceid, &new_faceid),
        "map old to new faceid");
    ref_cell_c2n(ref_cell, ref_cell_node_per(ref_cell), cell) = new_faceid;
  }

  ref_dict_free(ref_dict);
  ref_free(candidates);
  ref_free(face_norm);
  ref_free(face_cga);
  ref_free(face_box);
  ref_free(cad_cga);
  ref_free(cad_box);

#else
  printf("unable to %s, Need full EGADS linked.\n", __func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}
