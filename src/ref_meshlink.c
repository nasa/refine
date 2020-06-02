
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

#include "ref_cell.h"
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
  ProjectionDataObj projection_data = NULL;

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

  REIS(0, ML_createProjectionDataObj(geom_kernel, &projection_data), "make");
  ref_geom->meshlink_projection = (void *)projection_data;

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
          REF_INCREASE_LIMIT, "count edgeids");
      if (nedge > 1) {
        id = (REF_INT)ref_node_global(ref_grid_node(ref_grid), node);
        RSS(ref_geom_add(ref_geom, node, REF_GEOM_NODE, id, param), "node");
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_link(REF_GRID ref_grid, const char *block_name) {
  if (NULL == block_name) return REF_SUCCESS;
  printf("extracting mesh_model %s\n", block_name);
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  MeshAssociativityObj mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  MeshModelObj mesh_model;
  MLINT nsheet, msheet;
  REF_INT isheet;
  MeshSheetObj *sheet;
  MLINT nface, mface;
  MeshTopoObj *face;
  REF_INT iface;
  MLINT f2n[4], fn;
  MLINT sheet_gref;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  ParamVertexConstObj vert[3];
  MLINT nvert, vert_gref, mid;
  MLVector2D uv;
  REF_INT i, geom;
  REF_DBL param[2];
  char vref[REF_MESHLINK_MAX_STRING_SIZE];

  MLINT nstring, mstring;
  MeshTopoObj *string;
  REF_INT istring;
  MLINT string_gref;
  MLINT nedge, medge;
  MeshTopoObj *edge;
  REF_INT iedge;
  MLINT e2n[2], en;

  REIS(0, ML_getMeshModelByName(mesh_assoc, block_name, &mesh_model),
       "Error creating Mesh Model Object");

  nsheet = ML_getNumMeshSheets(mesh_model);
  printf("nsheet %" MLINT_FORMAT "\n", nsheet);
  ref_malloc(sheet, nsheet, MeshSheetObj);
  REIS(0, ML_getMeshSheets(mesh_model, sheet, nsheet, &msheet),
       "Error getting array of Mesh Sheets");
  REIS(nsheet, msheet, "sheet miscount");

  for (isheet = 0; isheet < nsheet; isheet++) {
    REIS(0, ML_getMeshTopoGref(sheet[isheet], &sheet_gref),
         "Error getting Mesh Sheet gref");
    nface = ML_getNumSheetMeshFaces(sheet[isheet]);
    printf("nface %" MLINT_FORMAT " for sheet %d with gref %" MLINT_FORMAT "\n",
           nface, isheet, sheet_gref);
    ref_malloc(face, nface, MeshTopoObj);
    REIS(0, ML_getSheetMeshFaces(sheet[isheet], face, nface, &mface),
         "Error getting array of Mesh Sheet Mesh Faces");
    REIS(nface, mface, "face miscount");
    for (iface = 0; iface < nface; iface++) {
      REIS(0, ML_getFaceInds(face[iface], f2n, &fn), "Error Mesh Face Ind");
      REIS(3, fn, "face ind miscount");
      RSS(ref_node_local(ref_node, f2n[0] - 1, &(nodes[0])), "g2l");
      RSS(ref_node_local(ref_node, f2n[1] - 1, &(nodes[1])), "g2l");
      RSS(ref_node_local(ref_node, f2n[2] - 1, &(nodes[2])), "g2l");
      nodes[3] = (REF_INT)sheet_gref;
      RSS(ref_cell_with(ref_grid_tri(ref_grid), nodes, &cell),
          "tri for sheet missing");
      ref_cell_c2n(ref_grid_tri(ref_grid), 3, cell) = nodes[3];
      REIS(0, ML_getParamVerts(face[iface], vert, 3, &nvert),
           "Error Face Vert");
      for (i = 0; i < 3; i++) {
        REIS(0,
             ML_getParamVertInfo(vert[i], vref, REF_MESHLINK_MAX_STRING_SIZE,
                                 &vert_gref, &mid, uv),
             "Error Face Vert");
        param[0] = uv[0];
        param[1] = uv[1];
        RSS(ref_geom_add(ref_geom, nodes[i], REF_GEOM_FACE, nodes[3], param),
            "face uv");
        RSS(ref_geom_find(ref_geom, nodes[i], REF_GEOM_FACE, nodes[3], &geom),
            "find");
        ref_geom_gref(ref_geom, geom) = (REF_INT)vert_gref;
      }
    }
    ref_free(face);
  }
  ref_free(sheet);

  nstring = ML_getNumMeshStrings(mesh_model);
  printf("nstring %" MLINT_FORMAT "\n", nstring);
  ref_malloc(string, nstring, MeshTopoObj);
  REIS(0, ML_getMeshStrings(mesh_model, string, nstring, &mstring),
       "Error getting array of Mesh Strings");
  REIS(nstring, mstring, "string miscount");
  for (istring = 0; istring < nstring; istring++) {
    REIS(0, ML_getMeshTopoGref(string[istring], &string_gref),
         "Error getting Mesh String gref");
    nedge = ML_getNumStringMeshEdges(string[istring]);
    printf("nedge %" MLINT_FORMAT " for string %d with gref %" MLINT_FORMAT
           "\n",
           nedge, istring, string_gref);
    ref_malloc(edge, nedge, MeshTopoObj);
    REIS(0, ML_getStringMeshEdges(string[istring], edge, nedge, &medge),
         "Error getting array of Mesh String Mesh Edges");
    REIS(nedge, medge, "edge miscount");
    for (iedge = 0; iedge < nedge; iedge++) {
      REIS(0, ML_getEdgeInds(edge[iedge], e2n, &en), "Error Mesh Edge Ind");
      REIS(2, en, "edge ind miscount");
      RSS(ref_node_local(ref_node, e2n[0] - 1, &(nodes[0])), "g2l");
      RSS(ref_node_local(ref_node, e2n[1] - 1, &(nodes[1])), "g2l");
      nodes[2] = (REF_INT)string_gref;
      RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &cell),
          "add edg for string");
      REIS(0, ML_getParamVerts(edge[iedge], vert, 2, &nvert),
           "Error Face Vert");
      for (i = 0; i < 2; i++) {
        REIS(0,
             ML_getParamVertInfo(vert[i], vref, REF_MESHLINK_MAX_STRING_SIZE,
                                 &vert_gref, &mid, uv),
             "Error Face Vert");
        param[0] = uv[0];
        param[1] = uv[1];
        RSS(ref_geom_add(ref_geom, nodes[i], REF_GEOM_EDGE, nodes[2], param),
            "face uv");
        RSS(ref_geom_find(ref_geom, nodes[i], REF_GEOM_EDGE, nodes[2], &geom),
            "find");
        ref_geom_gref(ref_geom, geom) = (REF_INT)vert_gref;
      }
    }
    ref_free(edge);
  }
  ref_free(string);

  { /* mark cad nodes */
    REF_INT node;
    REF_INT edges[2];
    REF_INT n, id;
    each_ref_node_valid_node(ref_node, node) {
      RXS(ref_cell_id_list_around(ref_grid_edg(ref_grid), node, 2, &n, edges),
          REF_INCREASE_LIMIT, "count edgeids");
      if (n > 1) {
        id = (REF_INT)ref_node_global(ref_grid_node(ref_grid), node);
        RSS(ref_geom_add(ref_geom, node, REF_GEOM_NODE, id, param), "node");
      }
    }
  }

  RSS(ref_geom_constrain_all(ref_grid), "constrain");

#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_constrain(REF_GRID ref_grid, REF_INT node) {
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT type, item, geom;
  REF_BOOL is_node, is_edge, is_face;
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  ProjectionDataObj projection_data = NULL;
  GeometryGroupObj geom_group = NULL;
  MLVector3D point;
  MLVector3D projected_point;
  MLVector2D uv;
  char entity_name[REF_MESHLINK_MAX_STRING_SIZE];
  MLINT gref;
  REF_DBL node_xyz[3];

  RNS(ref_geom->meshlink, "meshlink NULL");
  mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  projection_data = (ProjectionDataObj)(ref_geom->meshlink_projection);

  RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &is_node), "node");
  if (is_node) {
    node_xyz[0] = ref_node_xyz(ref_node, 0, node);
    node_xyz[1] = ref_node_xyz(ref_node, 1, node);
    node_xyz[2] = ref_node_xyz(ref_node, 2, node);
  }

  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
      gref = (MLINT)ref_geom_id(ref_geom, geom);
      point[0] = ref_node_xyz(ref_node, 0, node);
      point[1] = ref_node_xyz(ref_node, 1, node);
      point[2] = ref_node_xyz(ref_node, 2, node);

      REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
      REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");
      REIS(0, ML_projectPoint(geom_kernel, geom_group, point, projection_data),
           "prj");
      REIS(0,
           ML_getProjectionInfo(geom_kernel, projection_data, projected_point,
                                uv, entity_name, REF_MESHLINK_MAX_STRING_SIZE),
           "info");

      ref_node_xyz(ref_node, 0, node) = projected_point[0];
      ref_node_xyz(ref_node, 1, node) = projected_point[1];
      ref_node_xyz(ref_node, 2, node) = projected_point[2];
    }
  }

  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
      gref = (MLINT)ref_geom_id(ref_geom, geom);
      point[0] = ref_node_xyz(ref_node, 0, node);
      point[1] = ref_node_xyz(ref_node, 1, node);
      point[2] = ref_node_xyz(ref_node, 2, node);

      REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
      REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");
      REIS(0, ML_projectPoint(geom_kernel, geom_group, point, projection_data),
           "prj");
      REIS(0,
           ML_getProjectionInfo(geom_kernel, projection_data, projected_point,
                                uv, entity_name, REF_MESHLINK_MAX_STRING_SIZE),
           "info");

      ref_node_xyz(ref_node, 0, node) = projected_point[0];
      ref_node_xyz(ref_node, 1, node) = projected_point[1];
      ref_node_xyz(ref_node, 2, node) = projected_point[2];
    }
  }

  if (is_node) {
    ref_node_xyz(ref_node, 0, node) = node_xyz[0];
    ref_node_xyz(ref_node, 1, node) = node_xyz[1];
    ref_node_xyz(ref_node, 2, node) = node_xyz[2];
  }

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
  char entity_name[REF_MESHLINK_MAX_STRING_SIZE];
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

  RNS(ref_geom->meshlink, "meshlink NULL");
  mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  projection_data = (ProjectionDataObj)(ref_geom->meshlink_projection);

  for (i = 0; i < 3; i++) {
    center_point[i] = (1.0 / 3.0) * (ref_node_xyz(ref_node, i, nodes[0]) +
                                     ref_node_xyz(ref_node, i, nodes[1]) +
                                     ref_node_xyz(ref_node, i, nodes[2]));
  }
  gref = (MLINT)(id);

  REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
  REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");

  REIS(0,
       ML_projectPoint(geom_kernel, geom_group, center_point, projection_data),
       "prj");
  REIS(0,
       ML_getProjectionInfo(geom_kernel, projection_data, projected_point, uv,
                            entity_name, REF_MESHLINK_MAX_STRING_SIZE),
       "info");
  if (REF_FALSE) printf(" pre to %f %f of %s\n", uv[0], uv[1], entity_name);
  REIS(0,
       ML_evalCurvatureOnSurface(geom_kernel, uv, entity_name, eval_point,
                                 dXYZdU, dXYZdV, d2XYZdU2, d2XYZdUdV, d2XYZdV2,
                                 normal, principalV, &minCurvature,
                                 &maxCurvature, &avg, &gauss, &orientation),
       "eval");

  area_sign = 1.0;
  if (ML_ORIENT_SAME == orientation) area_sign = -1.0;

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
REF_STATUS ref_meshlink_face_curvature(REF_GEOM ref_geom, REF_INT geom,
                                       REF_DBL *kr, REF_DBL *r, REF_DBL *ks,
                                       REF_DBL *s) {
#ifdef HAVE_MESHLINK
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  ProjectionDataObj projection_data = NULL;
  GeometryGroupObj geom_group = NULL;
  MLINT gref, n_entity;
  MLVector3D projected_point;
  MLVector2D uv;
  char entity_name[REF_MESHLINK_MAX_STRING_SIZE];
  MLVector3D eval_point;
  MLVector3D dXYZdU;    /* First partial derivative */
  MLVector3D dXYZdV;    /* First partial derivative */
  MLVector3D d2XYZdU2;  /* Second partial derivative */
  MLVector3D d2XYZdUdV; /* Second partial derivative */
  MLVector3D d2XYZdV2;  /* Second partial derivative */
  MLVector3D normal;
  MLVector3D principalV;
  MLREAL minCurvature;
  MLREAL maxCurvature;
  MLREAL avg;
  MLREAL gauss;
  MLORIENT orientation;

  REIS(REF_GEOM_FACE, ref_geom_type(ref_geom, geom), "face geom expected");
  RNS(ref_geom->meshlink, "meshlink NULL");
  mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  projection_data = (ProjectionDataObj)(ref_geom->meshlink_projection);
  gref = (MLINT)ref_geom_gref(ref_geom, geom);

  REIS(0, ML_getActiveGeometryKernel(mesh_assoc, &geom_kernel), "kern");
  REIS(0, ML_getGeometryGroupByID(mesh_assoc, gref, &geom_group), "grp");
  REIB(0,
       ML_getEntityNames(geom_group, entity_name, 1,
                         REF_MESHLINK_MAX_STRING_SIZE, &n_entity),
       "grp", { printf("gref %" MLINT_FORMAT "\n", gref); });
  REIS(1, n_entity, "single entity expected");
  uv[0] = ref_geom_param(ref_geom, 0, geom);
  uv[1] = ref_geom_param(ref_geom, 1, geom);
  REIS(0,
       ML_evalCurvatureOnSurface(geom_kernel, uv, entity_name, eval_point,
                                 dXYZdU, dXYZdV, d2XYZdU2, d2XYZdUdV, d2XYZdV2,
                                 normal, principalV, &minCurvature,
                                 &maxCurvature, &avg, &gauss, &orientation),
       "eval");
  r[0] = principalV[0];
  r[1] = principalV[1];
  r[2] = principalV[2];
  ref_math_cross_product(normal, r, s);
  *kr = minCurvature;
  *ks = maxCurvature;

  return REF_SUCCESS;
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  SUPRESS_UNUSED_COMPILER_WARNING(geom);
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

REF_STATUS ref_meshlink_close(REF_GRID ref_grid) {
#ifdef HAVE_MESHLINK
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  ProjectionDataObj projection_data = NULL;
  if (NULL != ref_geom->meshlink_projection) {
    projection_data = (ProjectionDataObj)(ref_geom->meshlink_projection);
    ML_freeProjectionDataObj(&projection_data);
  }
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif
  return REF_SUCCESS;
}
