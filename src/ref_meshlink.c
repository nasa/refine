
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

#ifdef HAVE_MESHLINK
static REF_STATUS ref_meshlink_tattle_point(REF_INT node,
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

  if (0 != ML_getMeshPointInfo(mesh_assoc, mesh_point, ref, 1024, name, 1024,
                               &gref, &mid, attIDs, sizeAttIDs, &numAttIDs,
                               &paramVert)) {
    printf("evaluateParamPoint: bad point info\n");
  } else {
    printf("node %d geom ref %" MLINT_FORMAT " ref %s name %s\n", node, gref,
           ref, name);
  }
  return REF_SUCCESS;
}
#endif

REF_STATUS ref_meshlink_open(REF_GRID ref_grid, const char *xml_filename,
                             const char *block_name) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  if (NULL == xml_filename || NULL == block_name) return REF_SUCCESS;
#ifdef HAVE_MESHLINK
#define MAX_STRING_SIZE 256
  MeshAssociativityObj mesh_assoc;
  MeshModelObj mesh_model;
  GeometryKernelObj geom_kernel = NULL;
  MLINT iFile;
  MLINT numGeomFiles;
  MeshLinkFileConstObj geom_file;
  char geom_fname[MAX_STRING_SIZE];

  REIS(0, ML_createMeshAssociativityObj(&mesh_assoc),
       "Error creating Mesh Associativity Object");
  printf("have mesh_assoc\n");
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

  printf("extracting mesh_model %s\n", block_name);
  REIS(0, ML_getMeshModelByName(mesh_assoc, block_name, &mesh_model),
       "Error creating Mesh Model Object");
  printf("have mesh model\n");

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

  {
    MeshPointObj mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findLowestTopoPointByInd(mesh_model, i, &mesh_point)) {
        n = i;
        break;
      } else {
        RSS(ref_meshlink_tattle_point(i - 1, mesh_assoc, mesh_point), "tattle");
      }
    }
    printf("%d numpoints\n", n);
  }

  {
    MeshPointObj mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findHighestTopoPointByInd(mesh_model, i, &mesh_point)) {
        n = i;
        break;
      } else {
        RSS(ref_meshlink_tattle_point(i - 1, mesh_assoc, mesh_point), "tattle");
      }
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
    REF_NODE ref_node = ref_grid_node(ref_grid);
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

#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_close(REF_GRID ref_grid) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);

  return REF_SUCCESS;
}
