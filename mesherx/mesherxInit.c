
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

#include <stdlib.h>
#include "mesherxInit.h"
#include "gridfiller.h"
#ifdef HAVE_SDK
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#include "MeatLib/UGPatch.h"
#else
#include "FAKEGeom.h"
#endif

Layer *mesherxInit(int vol, int maxNodes)
{
  Grid *grid;
  Layer *layer;

  grid = gridFillFromPart( vol, maxNodes );
  if (NULL == grid){
    printf("ERROR: MesherX_DiscretizeVolume: could not fill grid with part.\n");
    return NULL;
  }

  layer = layerFormAdvancingLayerWithCADGeomBCS( vol, grid );
  layerFindParentGeomEdges(layer);
  gridThawAll(grid);

  return layer;
}

Layer *layerFormAdvancingLayerWithCADGeomBCS( int vol, Grid *grid )
{
  int nFrontFaces, frontFaces[10000];
  int face;
  UGPatchPtr upp;

  int    loop,edge,current;
  int    nloop, edgeindex;
  int    *nedge;
  int    *edges;
  double uv[4];
  int    self,other;
  int    left,right;

  Layer *layer;

  layer = layerCreate( grid );
  if (NULL == layer) printf("ERROR layerCreate failed: %s: %d\n",
			    __FILE__,__LINE__);
  nFrontFaces =0;
  for (face=1;face<=gridNGeomFace(grid);face++){
    upp = CADGeom_FaceGrid(vol,face);
    if (NULL == upp) printf("ERROR CADGeom_FaceGrid(%d,%d) failed: %s: %d\n",
			    vol,face,__FILE__,__LINE__);
    if ( NULL != UGPatch_BC(upp) &&
	 BC_NOSLIP == GeoBC_GenericType(UGPatch_BC(upp))){
      printf("face %d is added to front.\n",face);
      frontFaces[nFrontFaces] = face;
      nFrontFaces++;
    }
  }

  layerPopulateAdvancingFront(layer, nFrontFaces, frontFaces);
  
  for (face=1;face<=gridNGeomFace(grid);face++){
    upp = CADGeom_FaceGrid(vol,face);
    if ( NULL == UGPatch_BC(upp) ||
	 BC_NOSLIP != GeoBC_GenericType(UGPatch_BC(upp))){
      if( !CADGeom_GetFace(vol,face,uv,&nloop,&nedge,&edges) ) {/* Face Info */
        ErrMgr_Set(__FILE__,__LINE__,
		   "%s\nCould NOT get Volume %d, Face %d Info",
		   ErrMgr_GetErrStr(),vol,face);
        return( NULL );
      }

      edgeindex = 0;
      for( loop=0; loop<nloop; loop++ ) {                 /* Each Loop */
        for( current=0; current<nedge[loop]; current++) { /* Each Edge */
	  self = other = -1;
	  if( edges[edgeindex*2+1] > 0 )
            CADTopo_EdgeFaces(vol,edges[edgeindex*2],&self,&other);
          else
            CADTopo_EdgeFaces(vol,edges[edgeindex*2],&other,&self);
	  upp = CADGeom_FaceGrid(vol,other);
	  if ( NULL != UGPatch_BC(upp) &&
	       BC_NOSLIP == GeoBC_GenericType(UGPatch_BC(upp))){
	    layerConstrainNormal(layer,self);
	    printf("face %d is used to constrain normals.\n",face);
	  }
	  edgeindex++;
        }
      }
    }
  }

  for (edge=1;edge<=gridNGeomEdge(grid);edge++){
    CADTopo_EdgeFaces(vol,edge,&left,&right);
    if ( layerConstrainingGeometry(layer,left) &&
	 layerConstrainingGeometry(layer,right) ){
      layerConstrainNormal(layer,-edge);
      printf("edge %d is used to constrain normals.\n",edge);
    }
  }

  layerFeasibleNormals(layer,-1.0,-1.0);
  layerVisibleNormals(layer,-1.0,-1.0);

  return layer;
}
