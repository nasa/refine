
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <math.h>
#include <limits.h>         /* Needed in some systems for DBL_MAX definition */
#include <float.h>
#include "mesherx.h"
#include "mesherx_hla.h"
#include "mesherxInit.h"
#include "mesherxRebuild.h"
#include "grid.h"
#include "gridfiller.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#include "Goolache/MeshMgr.h"
#include "MeatLib/UGPatch.h"
#include "MeatLib/ErrMgr.h"
#include "MeatLib/GeoBC.h"

#include "geometricStretch.h"
#include "layerNormalExtra.h"

// a file parameter that has file scope

int MesherX_DiscretizeVolumeHLA(int maxNodes,
				double scale,
				double scalev,
				char *project,
				char *outputName,
				GridBool mixedElement,
				GridBool blendElement,
				GridBool qualityImprovement,
				GridBool bil )
{
  char outputProject[256];
  int vol=1;
  Grid *grid;
  Layer *layer;
  int i;
  double rate;
  int nLayer;
  double origin[3] = {0.0, -1.5, 0.0};
  double direction[3] = {0, 1, 0};

  int startRateIncrease;
  double rateOfRate;
  double rateMax;
  double rateMin;
  double hFirst;
  double hLayer;
  double trSpacing;
  double bgSpacingJet;
  double bgSpacingTunnel;

  double rate1,rate2;
  int nLayer1,nLayer2;
  int numSmooIter;
  int numSmooIterNP[4];
  double maxRate;

  layer = mesherxInit(vol, maxNodes);
  if( NULL == layer ) return 0;
  grid = layerGrid(layer);

// parameters controlling rate (in general)


  startRateIncrease=10000*scalev;
  rateOfRate = 1.02;
  rateMax = 1.3;
  rateMin = 1.1;
  hFirst = 0.005*scalev;
  hLayer = 150.0*hFirst;
  trSpacing = 0.9;			// spacing matching criterion
  bgSpacingJet = 0.25*scale*trSpacing;
  bgSpacingTunnel = 2.4*scale*trSpacing;

  if (bil) {
    nLayer = (int)(20.0/scalev);
    rate = exp(scalev*log(1.05));
  }else{
//    nLayer = (int)(20.0/scalev);
//    rate = exp(scalev*log(1.10));

// for SyJet

// inner region
    rate1 = MAX( 1.01, RateOfGeometricStretch( hLayer, hFirst, bgSpacingJet));
    nLayer1 = NptsOfGeometricStretch( hLayer, hFirst, rate1 );
    printf( "NumPts and rate in cavity: %d %f\n", nLayer1, rate1 );

// inner region
    rate2 = MAX( 1.01, RateOfGeometricStretch( 24.0, hFirst, bgSpacingTunnel ));
    nLayer2 = NptsOfGeometricStretch( 24.0, hFirst, rate2 );
    printf( "NumPts and rate in BL: %d %f\n", nLayer2, rate2 );

    rate = MIN( rate1, rate2 );
    rate = MIN( rateMax, MAX( rateMin, rate ));

    nLayer1 = MAX( nLayer1, NptsOfGeometricStretch( hLayer, hFirst, rate ));
    nLayer2 = MAX( nLayer2, NptsOfGeometricStretch( 24.0, hFirst, rate ));
    nLayer = MAX( nLayer1, nLayer2 );


    if( nLayer < 3 ) exit(0);
  }

  printf("rate is set to %10.5f for %d layers\n",rate,nLayer);

//  nLayer /= 2;
//  printf("number of layers is reduced to %d.\n",nLayer);


  if (mixedElement) layerToggleMixedElementMode(layer);

// streamwise direction

  direction[0] = 1.0;
  direction[1] = 0.0;
  direction[2] = 0.0;

  layerAssignPolarGrowthHeight(layer, scale, scalev,
				 direction );


  numSmooIterNP[0] = 0.0;
  numSmooIterNP[1] = (int)(0.5+10.0/scalev);
  numSmooIterNP[2] = numSmooIterNP[1];
  numSmooIterNP[3] = 0.0;
  layerSmoothNormalProperty( layer, numSmooIterNP, 0.5, TRUE );

  layerSaveInitialNormalHeight(layer);
  layerComputeNormalRateWithBGSpacing2(layer,trSpacing);

  numSmooIter = (int)(0.5+10.0/scalev);
  layerSmoothRate(layer, numSmooIter, 0.5, TRUE );


  /* only needed for formAdvancingFront freeze distant volume nodes */
  gridThawAll(grid);
  layerFindParentGeomEdges(layer);
  if (bil) {
    layerAssignPolynomialNormalHeight(layer, 0.002, 0.01, 2.0, 
				      origin, direction );
  }

  if (blendElement) {
    printf( "inserting blends...\n");
    layerBlend(layer,-1.0); 

/*
    printf( "extrude blends...\n");
    origin[0] = 1.0;
    origin[1] = 0.0;
    origin[2] = 0.0;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    layerCreateWakeWithBGSpacing(layer, origin, direction, 1.0 );

    origin[0] = -0.01;
    origin[1] = 0.0;
    origin[2] = 0.0;
    direction[0] = 1.0;
    direction[1] = 0.0;
    direction[2] = 0.0;
    layerAssignPolynomialNormalHeight(layer, 1.5e-4, 1.0e-3, 1.0, 
				      origin, direction );
    origin[0] = 1.0;
    layerAssignPolynomialNormalHeight(layer, 1.15e-3, 1.0e-2, 2.0, 
				      origin, direction );

    layerScaleNormalHeight(layer,scalev);
*/

  }

  i=0;
  while (i<nLayer
	&& layerNNormal(layer)>layerTerminateNormalWithBGSpacing(layer,trSpacing,2.9)
	&& layerNNormal(layer)>layerTerminateNormalWithLength(layer,0.90) 
		 ) {

    if( i > 5.0/scalev ) layerRelaxNormalDirection(layer, 1, 0.5);
    layerWriteTecplotFrontWithData(layer, i);
    layerAdvance(layer, FALSE);

    maxRate = 1.0;
    if( i > startRateIncrease ) {
        maxRate = rateMax;
        printf(" increasing rate towards %f ", maxRate );
    }
    layerSetNormalHeightWithRateAcceleration(layer,maxRate);


    printf("advance layer %d rate %f\n",i,rate);
    i++;
  }
  layerWriteTecplotFrontWithData(layer, i);

  if( layer != layerRebuildInterior(layer,vol) ) return 0;

  if ( qualityImprovement ){
    printf( " -- QUALITY IMPROVEMENT\n");
    layerThaw(layer);
    printf( "minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));
    for (i=0;i<3;i++){
      printf( "edge swapping grid...\n");gridSwap(grid);
      printf( "minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));
    }
  }

  printf( "total grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  printf( " -- DUMP PART\n");

  if ( NULL != outputName ) {

    if (mixedElement) {
      sprintf( outputProject,"%s.ugrid",outputName);
      printf( "writing output AFL3R file %s\n",outputProject);
      gridExportAFLR3( grid, outputProject  );
    }else{
      sprintf( outputProject,"%s.fgrid",outputName);
      printf( "writing output FAST file %s\n",outputProject);
      gridExportFAST( grid, outputProject  );
    }
  }

  if (!bil) {
    if ( outputName == NULL ) {
      gridSavePart( grid, NULL );
    }else{
      sprintf( outputProject,"%s",outputName);
      printf( "writing output GridEx/CADGeom/CAPRI project %s\n",outputProject);
      gridSavePart( grid, outputProject );
    }
  }

  return 1;
}
