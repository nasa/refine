
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>         /* Needed in some systems for DBL_MAX definition */
#include <float.h>
#include "mesherx.h"
#include "mesherxInit.h"
#include "mesherxRebuild.h"
#include "grid.h"
#include "gridfiller.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#include "Goolache/MeshMgr.h"
#include "Goolache/UGPatch.h"
#include "MeatLib/ErrMgr.h"
#include "MeatLib/GeoBC.h"

int MesherX_DiscretizeVolume( int maxNodes, double scale, char *project,
			      bool mixedElement,
			      bool blendElement,
			      bool qualityImprovement,
			      bool copyGridY,
			      bool bil )
{
  char outputProject[256];
  int vol=1;
  Grid *grid;
  Layer *layer;
  int i;
  double h;
  double rate;
  int nLayer;
  int face;
  double gapHeight;
  double origin[3] = {0.0, -1.5, 0.0};
  double direction[3] = {0, 1, 0};

  layer = mesherxInit(vol, maxNodes);
  if (NULL == layer) return 0; 
  grid = layerGrid(layer);

  /* case dependant */

  if (bil) {
    nLayer = (int)(20.0/scale);
    rate = exp(scale*log(1.05));
  }else{
    nLayer = (int)(100.0/scale);
    rate = exp(scale*log(1.20));
  }

  printf("rate is set to %10.5f for %d layers\n",rate,nLayer);

  if (mixedElement) layerToggleMixedElementMode(layer);

  if (bil ) { 
    layerAssignPolynomialNormalHeight(layer, 0.002, 0.01, 2.0, 
				      origin, direction );
  }else{
  }

  origin[0] = -0.0001;
  origin[1] = 0.0;
  origin[2] = 0.0;
  direction[0] = 1.0;
  direction[1] = 0.0;
  direction[2] = 0.0;
  layerSetPolynomialMaxHeight(layer, 0.71, 0.05, 1.0, 
			      origin, direction );
  origin[0] = 1.2;
  layerSetPolynomialMaxHeight(layer, 0.07, 0.1, 1.0, 
			      origin, direction );
  origin[0] = 1.7;
  layerSetPolynomialMaxHeight(layer, 0.12, 0.04, 1.0, 
			      origin, direction );
  origin[0] = -0.0001;
  layerAssignPolynomialNormalHeight(layer, 1.1e-3, 2.0e-3, 2.0,
                                    origin, direction );
  origin[0] = 1.1;
  layerAssignPolynomialNormalHeight(layer, 5.0e-5, 0.0, 1.0,
                                    origin, direction );
  origin[0] = 1.5;
  layerAssignPolynomialNormalHeight(layer, 5.0e-5, 1.0e-3, 2.0,
                                    origin, direction );
  origin[0] = 2.5;
  layerAssignPolynomialNormalHeight(layer, 0.00105, 1.5e-4, 1.0,
                                    origin, direction );
  layerScaleNormalHeight(layer,scale);
  layerSaveInitialNormalHeight(layer);
  layerComputeNormalRateWithBGSpacing(layer,1.0);

  i=0;
  while (i<nLayer &&
	 layerNNormal(layer)>layerTerminateNormalWithLength(layer,1.0)){

    layerSmoothNormalDirection(layer);
    layerSetNormalHeightWithMaxRate(layer,rate);
    layerTerminateCollidingFronts(layer);
    layerAdvance(layer);
    layerWriteTecplotFrontGeometry(layer);
    printf("advance layer %d rate %f\n",i,rate);
    i++;
  }

  /* case dep */

  if ( layer != layerRebuildInterior(layer,vol) ) return 0;

  if ( qualityImprovement ){
    printf(" -- QUALITY IMPROVEMENT\n");
    layerThaw(layer);
    printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));
    for (i=0;i<3;i++){
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));
    }
  }

  printf("total grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  printf(" -- DUMP PART\n");

  if ( copyGridY ) {
    printf("copy grid about y=0.\n");
    gridCopyAboutY0(grid);
    printf("total grid size: %d nodes %d faces %d cells.\n",
	   gridNNode(grid),gridNFace(grid),gridNCell(grid));
  }

  if ( NULL != project ) {

    if (mixedElement) {
      sprintf(outputProject,"%s_MX.ugrid",project);
      printf("writing output AFL3R file %s\n",outputProject);
      gridExportAFLR3( grid, outputProject  );
    }else{
      sprintf(outputProject,"%s_MX.fgrid",project);
      printf("writing output FAST file %s\n",outputProject);
      gridExportFAST( grid, outputProject  );
    }
  }

  if (!bil && !copyGridY ) {
    if ( project == NULL ) {
      gridSavePart( grid, NULL );
    }else{
      sprintf(outputProject,"%s_MX",project);
      printf("writing output GridEx/CADGeom/CAPRI project %s\n",outputProject);
      gridSavePart( grid, outputProject );
    }
  }else{
    printf("skip save grid.\n");
  }

  printf("MesherX Done.\n");
  return 1;
}

Layer *layerComputeNormalRateWithBGSpacing(Layer *layer, double finalRatio)
{
  int normal, root;
  double xyz[3], length, normalDirection[3];
  double initialDelta, finalDelta, rate;
  double spacing[3], direction[9];

  for(normal=0;normal<layerNNormal(layer);normal++){
    root = layerNormalRoot(layer, normal );
    gridNodeXYZ(layerGrid(layer),root,xyz);
    layerNormalDirection(layer,normal,normalDirection);
    length = layerNormalMaxLength(layer,normal);
    xyz[0] += (length * normalDirection[0]);
    xyz[1] += (length * normalDirection[1]);
    xyz[2] += (length * normalDirection[2]);
    MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
    finalDelta = spacing[0]*finalRatio;
    initialDelta = layerNormalInitialHeight(layer,normal);
    rate = (finalDelta - initialDelta) / length;
    gridNodeXYZ(layerGrid(layer),root,xyz);
    layerSetNormalRate(layer, normal, rate);
  }
  return layer;
}

Layer *layerComputeInitialCellHeightWithBGSpacing(Layer *layer, double finalRatio)
{
  int normal, root;
  double xyz[3], length, normalDirection[3];
  double initialDelta, finalDelta, rate;
  double spacing[3], direction[9];

  for(normal=0;normal<layerNNormal(layer);normal++){
    root = layerNormalRoot(layer, normal );
    gridNodeXYZ(layerGrid(layer),root,xyz);
    layerNormalDirection(layer,normal,normalDirection);
    length = layerNormalMaxLength(layer,normal);
    xyz[0] += (length * normalDirection[0]);
    xyz[1] += (length * normalDirection[1]);
    xyz[2] += (length * normalDirection[2]);
    MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
    finalDelta = spacing[0]*finalRatio;
    rate = layerNormalRate(layer,normal);
    initialDelta = finalDelta / pow(rate,length*10);
    printf("init %15.10f %15.10f %15.10f %15.10f\n",rate, initialDelta,finalDelta,length);
    layerSetNormalInitialHeight(layer, normal, initialDelta);
  }
  return layer;
}

Layer *layerCreateWakeWithBGSpacing(Layer *layer, 
				    double *origin, double *direction, 
				    double length )
{
  int i;
  double xyz[3];
  double spacing[3], map[9];
  double extrude[3], extrusion;

  i=0;
  extrusion = 0.0;
  xyz[0] = origin[0];  xyz[1] = origin[1];  xyz[2] = origin[2];
  while (extrusion < length){
    i++;
    MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,map);
    printf("wake%5d length%8.3f spacing%10.5f xyz%8.3f%8.3f%8.3f\n",
	   i,extrusion,spacing[0],xyz[0],xyz[1],xyz[2]);
    extrude[0] = direction[0]*spacing[0];
    extrude[1] = direction[1]*spacing[0];
    extrude[2] = direction[2]*spacing[0];
    layerExtrudeBlend(layer,extrude[0],extrude[1],extrude[2]); 
    xyz[0] += extrude[0];  xyz[1] += extrude[1];  xyz[2] += extrude[2];
    extrusion += spacing[0];
  }

}

int layerTerminateNormalWithBGSpacing(Layer *layer, 
				      double normalRatio, double edgeRatio)
{
  int normal, root;
  double xyz[3];
  double spacing[3];
  double direction[9];
  double height;
  int triangle, normals[3];
  double edgeLength, center[3];
  int totalterm;

  if (layerNNormal(layer) == 0 ) return EMPTY;

  for (normal=0;normal<layerNNormal(layer);normal++){
    layerGetNormalHeight(layer,normal,&height);

    root = layerNormalRoot(layer, normal );
    gridNodeXYZ(layerGrid(layer),root,xyz);
    MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);

    if (height > normalRatio*spacing[0]) {     /* Assume Isotropic for now */
      layerTerminateNormal(layer, normal);
    }
  }

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleMaxEdgeLength(layer,triangle,&edgeLength );
    layerTriangleCenter(layer,triangle,center);

    MeshMgr_GetSpacing(&(center[0]),&(center[1]),&(center[2]),
		       spacing,direction);

    if ( edgeLength > edgeRatio*spacing[0]) { /* Assume Isotropic for now */
      layerTriangleNormals(layer, triangle, normals);
      layerTerminateNormal(layer, normals[0]);
      layerTerminateNormal(layer, normals[1]);
      layerTerminateNormal(layer, normals[2]);
    }
  }

  totalterm = layerNNormal(layer)-layerNActiveNormal(layer);
  printf("%d of %d normals terminted.\n",
	 totalterm,layerNNormal(layer) );
  return totalterm;
}


