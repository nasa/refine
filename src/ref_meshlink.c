
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
#endif

#include "ref_edge.h"

REF_STATUS ref_meshlink_open(REF_GRID ref_grid, const char *xml_filename) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  if (NULL == xml_filename) return REF_SUCCESS;
#ifdef HAVE_MESHLINK
#define MAX_STRING_SIZE 256
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  MeshAssociativityObj mesh_assoc;
  GeometryKernelObj geom_kernel = NULL;
  MLINT iFile;
  MLINT numGeomFiles;
  MeshLinkFileConstObj geom_file;
  char geom_fname[MAX_STRING_SIZE];

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
    REIS(0, ML_getFilename(geom_file, geom_fname, MAX_STRING_SIZE),
         "Error getting Geometry File Name");
    printf("geom file %" MLINT_FORMAT " %s\n", iFile, geom_fname);
    REIS(0, ML_readGeomFile(geom_kernel, geom_fname),
         "Error reading Geometry File");
  }

#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_cache(REF_GRID ref_grid, const char *block_name) {
  if (NULL == block_name) return REF_SUCCESS;
  printf("extracting mesh_model %s\n", block_name);
#ifdef HAVE_MESHLINK
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  MeshAssociativityObj mesh_assoc = (MeshAssociativityObj)(ref_geom->meshlink);
  MeshModelObj mesh_model;
  REF_INT node;

  REIS(0, ML_getMeshModelByName(mesh_assoc, block_name, &mesh_model),
       "Error creating Mesh Model Object");

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_cell_node_empty(ref_grid_tri(ref_grid), node)) {
      MeshPointObj edge_mesh_point = NULL;
      MeshPointObj face_mesh_point = NULL;
      char ref[256];
      char name[256];
      MLINT edge_gref, face_gref;
      MLINT mid;
      MLINT attIDs[24];
      MLINT numAttIDs;
      ParamVertexConstObj paramVert;
      REF_DBL param[2] = {0.0, 0.0};
      REF_INT id;
      REIS(0,
           ML_findLowestTopoPointByInd(mesh_model, node + 1, &edge_mesh_point),
           "low/edge");
      REIS(0,
           ML_findHighestTopoPointByInd(mesh_model, node + 1, &face_mesh_point),
           "high/face");
      REIS(0,
	   ML_getMeshPointInfo(mesh_assoc, edge_mesh_point, ref, 256, name, 256, &edge_gref,
			       &mid, attIDs, 24, &numAttIDs, &paramVert),
       "bad point info");
      REIS(0,
	   ML_getMeshPointInfo(mesh_assoc, face_mesh_point, ref, 256, name, 256, &face_gref,
			       &mid, attIDs, 24, &numAttIDs, &paramVert),
       "bad point info");
      id = face_gref;
      RSS(ref_geom_add(ref_geom, node, REF_GEOM_FACE, id, param),
          "face uv");
      if (edge_gref != face_gref) {
	id = edge_gref;
      RSS(ref_geom_add(ref_geom, node, REF_GEOM_EDGE, id, param),
          "edge t");
      }

    }
  }

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

REF_STATUS ref_meshlink_close(REF_GRID ref_grid) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);

  return REF_SUCCESS;
}
