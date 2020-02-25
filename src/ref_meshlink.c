
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
#include "MeshAssociativity.h" /* for debug info, not c */
#include "MeshAssociativity_c.h"
#include "MeshLinkParser_xerces_c.h"
#endif

#include "ref_edge.h"

REF_STATUS ref_meshlink_open(REF_GRID ref_grid, const char *xml_filename,
                             const char *block_name) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  if (NULL == xml_filename || NULL == block_name) return REF_SUCCESS;
#ifdef HAVE_MESHLINK
  MeshAssociativityObj mesh_assoc;
  MeshModelObj mesh_model;
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
  printf("extracted mesh_model numedges %d numfaces %d\n",
         ((MeshModel *)mesh_model)->getNumEdges(),
         ((MeshModel *)mesh_model)->getNumFaces());

  {
    MeshPointObj mesh_point = NULL;
    REF_INT i, n;
    n = 0;
    for (i = 1; i < 1000; i++) {
      if (0 != ML_findLowestTopoPointByInd(mesh_model, i, &mesh_point)) {
        n = i;
        break;
      } else {
        char ref[1024];
        char name[1024];
        MLINT gref;
        MLINT mid;
        const MLINT sizeAttIDs = 24;
        MLINT attIDs[24];
        MLINT numAttIDs;
        ParamVertexConstObj paramVert;

        if (0 != ML_getMeshPointInfo(mesh_assoc, mesh_point, ref, 1024, name,
                                     1024, &gref, &mid, attIDs, sizeAttIDs,
                                     &numAttIDs, &paramVert)) {
          printf("evaluateParamPoint: bad point info\n");
        } else {
          printf("name %s\n", name);
        }
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
    REF_BOOL tri_side;
    RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");
    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      node0 = ref_edge_e2n(ref_edge, 0, edge);
      node1 = ref_edge_e2n(ref_edge, 1, edge);
      RSS(ref_cell_has_side(ref_cell, node0, node1, &tri_side), "is tri side");
      if (tri_side) {
        MLINT edge_indexes[2];
        MeshEdgeObj mesh_edge = NULL;
        edge_indexes[0] = ref_node_global(ref_node, node0);
        edge_indexes[1] = ref_node_global(ref_node, node1);
        REIS(0,
             ML_findLowestTopoEdgeByInds(mesh_model, edge_indexes, (MLINT)2,
                                         &mesh_edge),
             "find edge");
      }
    }
    ref_edge_free(ref_edge);
  }

#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_close(REF_GRID ref_grid) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);

  return REF_SUCCESS;
}
