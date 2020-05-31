
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

#include "ref_meshlink.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MESHLINK
#define IS64BIT
#include "GeomKernel_Geode_c.h"
#include "MeshAssociativity_c.h"
#include "MeshLinkParser_xerces_c.h"
#define REF_MESHLINK_MAX_STRING_SIZE 256
#endif

#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_math.h"

REF_STATUS ref_meshlink_open(REF_GRID ref_grid, const char *xml_filename) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  if (NULL == xml_filename) return REF_SUCCESS;
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  MLINT iFile;
  MLINT numGeomFiles;
  MeshLinkFileConstObj geom_file;
  char geom_fname[REF_MESHLINK_MAX_STRING_SIZE];

  REIS(0, ML_createMeshAssociativityObj(&mesh_assoc),
       "Error creating Mesh Associativity Object");
  printf("have mesh_assoc\n");
  ref_geom->meshlink = (void *)mesh_assoc;
  /* Read Geometry-Mesh associativity */
  {
    /* NULL schema filename uses schemaLocation in meshlink file */
    const char *schema_filename = NULL;
    /* Xerces MeshLink XML parser */
    MeshLinkParserObj parser;
    REIS(0, ML_createMeshLinkParserXercesObj(&parser), "create parser");
    printf("validate %s\n", xml_filename);
    REIS(0, ML_parserValidateFile(parser, xml_filename, schema_filename),
         "validate");
    printf("parse %s\n", xml_filename);
    REIS(0, ML_parserReadMeshLinkFile(parser, xml_filename, mesh_assoc),
         "parse");
    ML_freeMeshLinkParserXercesObj(&parser);
  }
  printf("populated mesh_assoc\n");

  printf("extracting geom_kernel\n");
  REIS(0, ML_createGeometryKernelGeodeObj(&geom_kernel),
       "Error creating Geometry Kernel Object");
  printf("have geom kernel\n");

  printf("activate geode\n");
  REIS(0, ML_addGeometryKernel(mesh_assoc, geom_kernel),
       "Error adding Geometry Kernel Object");
  REIS(0, ML_setActiveGeometryKernelByName(mesh_assoc, "Geode"),
       "Error adding Geometry Kernel Object");
  printf("active geom kernel\n");

  numGeomFiles = ML_getNumGeometryFiles(mesh_assoc);
  printf("geom files %" MLINT_FORMAT "\n", numGeomFiles);
  for (iFile = 0; iFile < numGeomFiles; ++iFile) {
    REIS(0, ML_getGeometryFileObj(mesh_assoc, iFile, &geom_file),
         "Error getting Geometry File");
    REIS(0, ML_getFilename(geom_file, geom_fname, REF_MESHLINK_MAX_STRING_SIZE),
         "Error getting Geometry File Name");
    printf("geom file %" MLINT_FORMAT " %s\n", iFile, geom_fname);
    REIS(0, ML_readGeomFile(geom_kernel, geom_fname),
         "Error reading Geometry File");
  }

#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_parse(REF_GRID ref_grid, const char *geom_filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *f = NULL;
  char line[1024];
  REF_INT tri, ntri, edge, nedge, gref, new_cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  int status;
  REF_DBL param[2] = {0.0, 0.0};
  if (NULL == geom_filename) return REF_SUCCESS;
  printf("parsing %s\n", geom_filename);
  f = fopen(geom_filename, "r");
  if (NULL == (void *)f) printf("unable to open %s\n", geom_filename);
  RNS(f, "unable to open file");
  while (!feof(f)) {
    status = fscanf(f, "%s", line);
    if (EOF == status) break;
    REIS(1, status, "line read failed");

    if (0 == strcmp("sheet", line)) {
      REIS(2, fscanf(f, "%d %d", &ntri, &gref), "sheet size gref");
      printf("sheet ntri %d gref %d\n", ntri, gref);
      for (tri = 0; tri < ntri; tri++) {
        REIS(3, fscanf(f, "%d %d %d", &(nodes[0]), &(nodes[1]), &(nodes[2])),
             "tri nodes");
        (nodes[0])--;
        (nodes[1])--;
        (nodes[2])--;
        nodes[3] = gref;
        RSS(ref_cell_with(ref_grid_tri(ref_grid), nodes, &new_cell),
            "tri for sheet missing");
        ref_cell_c2n(ref_grid_tri(ref_grid), 3, new_cell) = nodes[3];
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_FACE, nodes[3], param),
            "face uv");
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_FACE, nodes[3], param),
            "face uv");
        RSS(ref_geom_add(ref_geom, nodes[2], REF_GEOM_FACE, nodes[3], param),
            "face uv");
      }
    }

    if (0 == strcmp("string", line)) {
      REIS(2, fscanf(f, "%d %d", &nedge, &gref), "string size gref");
      printf("sheet ntri %d gref %d\n", nedge, gref);
      for (edge = 0; edge < nedge; edge++) {
        REIS(2, fscanf(f, "%d %d", &(nodes[0]), &(nodes[1])), "edge nodes");
        (nodes[0])--;
        (nodes[1])--;
        nodes[2] = gref;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell),
            "add edg for string");
        RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE, nodes[2], param),
            "edge t");
        RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_EDGE, nodes[2], param),
            "edge t");
      }
    }
  }
  fclose(f);

  { /* mark cad nodes */
    REF_INT node;
    REF_INT edges[2];
    REF_INT id;
    each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
      RXS(ref_cell_id_list_around(ref_grid_edg(ref_grid), node, 2, &nedge,
                                  edges),
          REF_INCREASE_LIMIT, "count faceids");
      if (nedge > 1) {
        id = (REF_INT)ref_node_global(ref_grid_node(ref_grid), node);
        RSS(ref_geom_add(ref_geom, node, REF_GEOM_NODE, id, param), "node");
      }
    }
  }

  return REF_SUCCESS;
}

#ifdef HAVE_MESHLINK
static REF_STATUS ref_swap_same_faceid(REF_GRID ref_grid, REF_INT node0,
                                       REF_INT node1, REF_BOOL *same) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT ncell;
  REF_INT cells[2];
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT id0, id1;
  REF_BOOL has_edg, has_tri, has_qua;

  *same = REF_FALSE;

  RSB(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells),
      "more then two", {
        ref_node_location(ref_grid_node(ref_grid), node0);
        ref_node_location(ref_grid_node(ref_grid), node1);
      });

  REIB(2, ncell, "there should be two triangles for manifold", {
    ref_node_location(ref_grid_node(ref_grid), node0);
    ref_node_location(ref_grid_node(ref_grid), node1);
  });

  RSS(ref_cell_nodes(ref_cell, cells[0], nodes), "nodes tri0");
  id0 = nodes[ref_cell_node_per(ref_cell)];
  RSS(ref_cell_nodes(ref_cell, cells[1], nodes), "nodes tri1");
  id1 = nodes[ref_cell_node_per(ref_cell)];

  if (id0 == id1) {
    *same = REF_TRUE;
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}
#endif

REF_STATUS ref_meshlink_cache(REF_GRID ref_grid, const char *block_name) {
  if (NULL == block_name) return REF_SUCCESS;
  printf("extracting mesh_model %s\n", block_name);
#ifdef HAVE_MESHLINK
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  MeshAssociativityObj mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  MeshModelObj mesh_model;

  REIS(0, ML_getMeshModelByName(mesh_assoc, block_name, &mesh_model),
       "Error creating Mesh Model Object");

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

#ifdef HAVE_MESHLINK
static REF_STATUS ref_meshlink_tattle_point(REF_NODE ref_node, REF_INT node,
                                            MeshAssociativityObj mesh_assoc,
                                            MeshPointObj mesh_point) {
  char ref[1024];
  char name[1024];
  MLINT gref;
  MLINT mid;
  const MLINT sizeAttIDs = 24;
  MLINT attIDs[24];
  MLINT numAttIDs;
  ParamVertexConstObj paramVert;

  REIS(0,
       ML_getMeshPointInfo(mesh_assoc, mesh_point, ref, 1024, name, 1024, &gref,
                           &mid, attIDs, sizeAttIDs, &numAttIDs, &paramVert),
       "bad point info");
  {
    GeometryKernelObj geom_kernel = NULL;
    ProjectionDataObj projection_data = NULL;
    GeometryGroupObj geom_group = NULL;
    MLVector3D point;
    MLVector3D projected_point;
    MLVector2D uv;
    char entity_name[1024];
    REF_DBL dist;
    point[0] = ref_node_xyz(ref_node, 0, node);
    point[1] = ref_node_xyz(ref_node, 1, node);
    point[2] = ref_node_xyz(ref_node, 2, node);
    REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
    REIS(0, ML_createProjectionDataObj(geom_kernel, &projection_data), "make");
    REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");
    REIS(0, ML_projectPoint(geom_kernel, geom_group, point, projection_data),
         "prj");
    REIS(0,
         ML_getProjectionInfo(geom_kernel, projection_data, projected_point, uv,
                              entity_name, 1024),
         "info");
    dist = sqrt(pow(projected_point[0] - point[0], 2) +
                pow(projected_point[1] - point[1], 2) +
                pow(projected_point[2] - point[2], 2));
    printf("node %d geom ref %" MLINT_FORMAT "", node, gref);
    printf(" dist %e to %s\n", dist, entity_name);
    ML_freeProjectionDataObj(&projection_data);
  }
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_meshlink_examine(REF_GRID ref_grid, const char *block_name) {
  if (NULL == block_name) return REF_SUCCESS;
  printf("extracting mesh_model %s\n", block_name);
#ifdef HAVE_MESHLINK
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  MeshAssociativityObj mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  MeshModelObj mesh_model;

  REIS(0, ML_getMeshModelByName(mesh_assoc, block_name, &mesh_model),
       "Error creating Mesh Model Object");
  printf("have mesh model\n");
  /* low is edge, high is face */
  {
    MeshPointObj low_mesh_point = NULL;
    MeshPointObj high_mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findLowestTopoPointByInd(mesh_model, i, &low_mesh_point)) {
        n = i;
        break;
      }
      if (0 != ML_findHighestTopoPointByInd(mesh_model, i, &high_mesh_point)) {
        n = i;
        break;
      }
      RSS(ref_meshlink_tattle_point(ref_node, i - 1, mesh_assoc,
                                    low_mesh_point),
          "tattle");
      RSS(ref_meshlink_tattle_point(ref_node, i - 1, mesh_assoc,
                                    high_mesh_point),
          "tattle");
    }
    printf("%d numpoints\n", n);
  }

  {
    MeshPointObj mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findMeshEdgePointByInd(mesh_model, i, &mesh_point)) {
        n = i;
        break;
      }
    }
    printf("%d num edge points\n", n);
  }

  {
    MeshPointObj mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findMeshFacePointByInd(mesh_model, i, &mesh_point)) {
        n = i;
        break;
      }
    }
    printf("%d num face points\n", n);
  }

  {
    REF_CELL ref_cell = ref_grid_tri(ref_grid);
    REF_EDGE ref_edge;
    REF_INT edge, node0, node1;
    REF_INT n;
    REF_BOOL tri_side;
    RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");
    n = 0;
    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      node0 = ref_edge_e2n(ref_edge, 0, edge);
      node1 = ref_edge_e2n(ref_edge, 1, edge);
      RSS(ref_cell_has_side(ref_cell, node0, node1, &tri_side), "is tri side");
      if (tri_side) {
        MLINT edge_indexes[2];
        MeshEdgeObj mesh_edge = NULL;
        edge_indexes[0] = ref_node_global(ref_node, node0) + 1;
        edge_indexes[1] = ref_node_global(ref_node, node1) + 1;
        REIS(0,
             ML_findLowestTopoEdgeByInds(mesh_model, edge_indexes, (MLINT)2,
                                         &mesh_edge),
             "find edge");
        n++;
      }
    }
    printf("%d surface edges\n", n);
    ref_edge_free(ref_edge);
  }

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_constrain(REF_GRID ref_grid, REF_INT node) {
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT type, item, geom, gref_geom;
  REF_BOOL is_node, is_edge, is_face;
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  ProjectionDataObj projection_data = NULL;
  GeometryGroupObj geom_group = NULL;
  MLVector3D point;
  MLVector3D projected_point;
  MLVector2D uv;
  char entity_name[256];
  MLINT gref;

  RNS(ref_geom->meshlink, "meshlink NULL");
  mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &is_node), "node");
  if (is_node) return REF_SUCCESS; /* can't move geom node */

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &is_edge), "edge");
  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &is_face), "face");
  RAS(is_edge || is_face, "no edge or face for contraint");
  if (is_edge) {
    type = REF_GEOM_EDGE;
  } else {
    type = REF_GEOM_FACE;
  }
  gref_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      gref_geom = geom;
      break;
    }
  }
  RUS(REF_EMPTY, gref_geom, "can't find geom");
  gref = (MLINT)ref_geom_id(ref_geom, gref_geom);

  point[0] = ref_node_xyz(ref_node, 0, node);
  point[1] = ref_node_xyz(ref_node, 1, node);
  point[2] = ref_node_xyz(ref_node, 2, node);

  REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
  REIS(0, ML_createProjectionDataObj(geom_kernel, &projection_data), "make");
  REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");
  REIS(0, ML_projectPoint(geom_kernel, geom_group, point, projection_data),
       "prj");
  REIS(0,
       ML_getProjectionInfo(geom_kernel, projection_data, projected_point, uv,
                            entity_name, 256),
       "info");
  ML_freeProjectionDataObj(&projection_data);

  ref_node_xyz(ref_node, 0, node) = projected_point[0];
  ref_node_xyz(ref_node, 1, node) = projected_point[1];
  ref_node_xyz(ref_node, 2, node) = projected_point[2];

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  SUPRESS_UNUSED_COMPILER_WARNING(node);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_tri_norm_deviation(REF_GRID ref_grid, REF_INT *nodes,
                                           REF_DBL *dot_product) {
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  ProjectionDataObj projection_data = NULL;
  GeometryGroupObj geom_group = NULL;
  MLVector3D center_point;
  MLVector3D normal;
  REF_INT i, id;
  MLINT gref;
  REF_STATUS status;
  REF_DBL tri_normal[3];
  REF_DBL area_sign = 1.0;
  MLVector3D projected_point;
  MLVector2D uv;
  char entity_name[256];
  MLVector3D eval_point;
  MLVector3D dXYZdU;    /* First partial derivative */
  MLVector3D dXYZdV;    /* First partial derivative */
  MLVector3D d2XYZdU2;  /* Second partial derivative */
  MLVector3D d2XYZdUdV; /* Second partial derivative */
  MLVector3D d2XYZdV2;  /* Second partial derivative */
  MLVector3D principalV;
  MLREAL minCurvature;
  MLREAL maxCurvature;
  MLREAL avg;
  MLREAL gauss;
  MLORIENT orientation;

  id = nodes[3];
  RSS(ref_node_tri_normal(ref_grid_node(ref_grid), nodes, tri_normal),
      "tri normal");
  /* collapse attempts could create zero area, reject the step with -2.0 */
  status = ref_math_normalize(tri_normal);
  if (REF_DIV_ZERO == status) return REF_SUCCESS;
  RSS(status, "normalize");

  RNS(ref_geom->uv_area_sign, "uv_area_sign NULL");
  RNS(ref_geom->meshlink, "meshlink NULL");
  mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);

  for (i = 0; i < 3; i++) {
    center_point[i] = (1.0 / 3.0) * (ref_node_xyz(ref_node, i, nodes[0]) +
                                     ref_node_xyz(ref_node, i, nodes[1]) +
                                     ref_node_xyz(ref_node, i, nodes[2]));
  }
  gref = (MLINT)(id);

  REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
  REIS(0, ML_createProjectionDataObj(geom_kernel, &projection_data), "make");
  REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");

  REIS(0,
       ML_projectPoint(geom_kernel, geom_group, center_point, projection_data),
       "prj");
  REIS(0,
       ML_getProjectionInfo(geom_kernel, projection_data, projected_point, uv,
                            entity_name, 256),
       "info");
  if (REF_FALSE) printf(" pre to %f %f of %s\n", uv[0], uv[1], entity_name);
  REIS(0,
       ML_evalCurvatureOnSurface(geom_kernel, uv, entity_name, eval_point,
                                 dXYZdU, dXYZdV, d2XYZdU2, d2XYZdUdV, d2XYZdV2,
                                 normal, principalV, &minCurvature,
                                 &maxCurvature, &avg, &gauss, &orientation),
       "eval");

  ML_freeProjectionDataObj(&projection_data);

  area_sign = ref_geom->uv_area_sign[id - 1];

  *dot_product = area_sign * ref_math_dot(normal, tri_normal);

  if (REF_FALSE) {
    printf("surf %.3f %.3f %.3f disc %.3f %.3f %.3f\n", normal[0], normal[1],
           normal[2], tri_normal[0], tri_normal[1], tri_normal[2]);
  }

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  SUPRESS_UNUSED_COMPILER_WARNING(nodes);
  *dot_product = -2.0;
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_close(REF_GRID ref_grid) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);

  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_infer_orientation(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT min_id, max_id, id;
  REF_DBL normdev, min_normdev, max_normdev;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  RSS(ref_cell_id_range(ref_cell, ref_mpi, &min_id, &max_id), "id range");
  RAS(min_id > 0, "expected min_id greater then zero");
  ref_geom->nface = max_id;
  ref_malloc_init(ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 1.0);

  for (id = min_id; id <= max_id; id++) {
    min_normdev = 2.0;
    max_normdev = -2.0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (id != nodes[ref_cell_id_index(ref_cell)]) continue;
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "norm dev");
      min_normdev = MIN(min_normdev, normdev);
      max_normdev = MAX(max_normdev, normdev);
    }
    normdev = min_normdev;
    RSS(ref_mpi_min(ref_mpi, &normdev, &min_normdev, REF_DBL_TYPE), "mpi max");
    RSS(ref_mpi_bcast(ref_mpi, &min_normdev, 1, REF_DBL_TYPE), "min");
    normdev = max_normdev;
    RSS(ref_mpi_max(ref_mpi, &normdev, &max_normdev, REF_DBL_TYPE), "mpi max");
    RSS(ref_mpi_bcast(ref_mpi, &max_normdev, 1, REF_DBL_TYPE), "max");
    if (min_normdev > 1.5 && max_normdev < -1.5) {
      ref_geom->uv_area_sign[id - 1] = 0.0;
    } else {
      if (min_normdev < -max_normdev) ref_geom->uv_area_sign[id - 1] = -1.0;
      if (ref_mpi_once(ref_mpi))
        printf("gref %3d orientation%6.2f inferred from %6.2f %6.2f\n", id,
               ref_geom->uv_area_sign[id - 1], min_normdev, max_normdev);
    }
  }

  return REF_SUCCESS;
}
