
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
#include "MeshAssociativity_c.h"
#include "MeshLinkParser_xerces_c.h"
#endif

REF_STATUS ref_meshlink_open(REF_GEOM ref_geom, const char *filename) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
  if (NULL == filename) return REF_SUCCESS;
  printf("MeskLink with %s\n", filename);
#ifdef HAVE_MESHLINK
  MeshAssociativityObj mesh_assoc;
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
    REIS(0, ML_parserValidateFile(parser, filename, schema_filename),
         "validate");
    REIS(0, ML_parserReadMeshLinkFile(parser, filename, mesh_assoc), "parse");
    ML_freeMeshLinkParserXercesObj(&parser);
  }

#endif
  return REF_SUCCESS;
}

REF_STATUS ref_meshlink_close(REF_GEOM ref_geom) {
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);

  return REF_SUCCESS;
}
