

/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include "mesherxInit.h"
#include "gridfiller.h"
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#include "MeatLib/UGPatch.h"

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

  layerVisibleNormals(layer,-1.0,-1.0);

  return layer;
}
