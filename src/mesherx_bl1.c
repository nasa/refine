
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
				bool mixedElement,
				bool blendElement,
				bool qualityImprovement,
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

  double MMtoFt = 12.0*2.54;

  int startRateIncrease=100/scalev;
  double rateOfRate = 1.02;
  double rateMax = 1.0;
  double hFirst = 0.005*scalev;
  double hLayer = 150.0*hFirst;
  double trSpacing = 1.0;			// spacing matching criterion
  double bgSpacingJet = 0.25*scale*trSpacing;
  double bgSpacingTunnel = 2.4*scale*trSpacing;

  layer = mesherxInit(vol, maxNodes);
  if( NULL == layer ) return 0;
  grid = layerGrid(layer);
  if (mixedElement) layerToggleMixedElementMode(layer);

// parameters controlling rate (in general)
// parameters set for a coarse grid for domain in units of ft

  hFirst = hFirst/MMtoFt;
  rate = pow(1.3 , scalev);
  nLayer = NptsOfGeometricStretch( 0.9, hFirst, rate );
  printf( "NumPts and rate in cavity: %d %f\n", nLayer, rate );

  startRateIncrease= NptsOfGeometricStretch( 1.0/12.0, hFirst, rate );
  rateMax = 1.3;

  printf("rate is set to %10.5f for %d layers\n", rate, nLayer);
  printf("increasing rate towards %f after %d layers\n", rateMax, startRateIncrease );

// setting distribution of rate

  double power = 6.0/7.0;
  double constant = 0.4*hFirst;
  double slope = (hFirst-constant)/pow( 5.5-origin[0], power );
  double hh = hFirst0;
  double rr = rate;
  double yy = 0.0;

  double hFirst0 = constant + slope*pow(-origin[0],power);
  double hMax0   = 0.1 * 0.7;

  double L0 = LengthFromDyMax( hMax0/2.0, hFirst0, rate );
  int i=0;

  origin[0] = -0.2;
  origin[1] = 0.0;
  origin[2] = 0.0;
  direction[0] = 1.0;
  direction[1] = 0.0;
  direction[2] = 0.0;

/*
  {
    double ff = exp(-5.0*pow(L0, 0.5));
    double h00 = hMax0/2.0;
    printf( " i=%d,  ff=%f,  h00=%f,  L0=%f\n", i, ff, h00, L0);
    while( i<4 ) {
      ff = exp(-5.0*sqrt(L0));
      h00 = 0.025*ff + hMax0*(1.0-ff);
      L0 = LengthFromDyMax( h00, hFirst0, rate );
      i++;
      printf( " i=%d,  ff=%f,  h00=%f,  L0=%f\n", i, ff, h00, L0);
    }
  }
*/

  double L1 = MIN( 0.5, L0 );
  int nL0 = NptsOfGeometricStretch( L1, hFirst0, rate );

/*
  {
    int i=0;
    while( yy < hMax0 ) {
      yy += hh;
      if( i > startRateIncrease ) rr *= 1.02;
      hh *= rr;
      i++;
    }
    nLayer = MIN( nLayer, MIN( i, nL0 ));
  }
*/
  nLayer = MAX( nLayer, nL0 );
	
    
  printf( "rate = %f, h0 = %f, H0 = %f, L0 = %f, i = %d, nL0 = %d\n", 
           rate, hFirst0, hMax0, L0, i, nL0 );
  printf( "Changing number of layers to %d \n", nLayer );

  layerAssignPolynomialNormalHeight(layer, constant, slope, power, 
				      origin, direction );

// a drop dead upper bound on height

  constant = 0.8*(8.0/12.0);
  slope = (8.0/12.0 - constant)/pow( 5.5-origin[0], power );
  layerSetPolynomialMaxHeight(layer, constant, slope, power,
                                      origin, direction );

  layerSetPolynomialMaxHeight(layer, 0.8, 0.0, 1.0,
                                      origin, direction );

  layerSetAllNormalRate( layer, rate );
  layerSaveInitialNormalHeight( layer );

//  layerSaveInitialNormalHeight(layer);
//  layerComputeNormalRateWithBGSpacing(layer, trSpacing);

// omega = zero just to generate tecplot files

  layerSmoothRate(layer, 1, 0.0, TRUE );


  if (blendElement) {
    printf( "inserting blends...\n");
    layerBlend(layer); 

/*
// old code kept as an example
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

    printf( "split blends...\n");
    layerSplitBlend(layer); 

/*
    printf( "split blends...\n");
    layerSplitBlend(layer); 

    printf( "split blends...\n");
    layerSplitBlend(layer); 
    if (scalev < 0.55) layerSplitBlend(layer);
*/
  }

  i=0;
  while (i<nLayer
	&& layerNNormal(layer)>layerTerminateNormalWithBGSpacing(layer,trSpacing,1.9)
	&& layerNNormal(layer)>layerTerminateNormalWithLength(layer,1.0) 
//	&& layerNNormal(layer)>layerTerminateNormalIfInX(layer,-1, 0.5 ) 
		 ) {

    if( i > 5.0/scalev ) layerRelaxNormalDirection(layer, 1, 0.5);
    layerWriteTecplotFront(layer,i);
    layerAdvance(layer);

    double maxRate = 1.0;
    if( i > startRateIncrease ) {
        maxRate = rateMax;
        printf(" increasing rate towards %f ", maxRate );
    }
    layerSetNormalHeightWithRateAcceleration(layer, maxRate);


    printf("advance layer %d rate %f\n",i,rate);
    i++;
  }
  layerWriteTecplotFront(layer,i);

  if ( layer != layerRebuildInterior(layer,vol) ) return 0;

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

    if ( outputName == NULL ) {
      gridSavePart( grid, NULL );
    }else{
      sprintf( outputProject,"%s",outputName);
      printf( "writing output GridEx/CADGeom/CAPRI project %s\n",outputProject);
      gridSavePart( grid, outputProject );
    }

  return 1;
}



Layer *layerComputeNormalRateWithBGSpacing2(Layer *layer, double finalRatio)
{
  int normal, root, maxNpts, npts;
  
  double xyz[3], xyz0[3], length, normalDirection[3];
  double initialDelta, finalDelta, rate, maxRate;
  double spacing[3], direction[9];

  maxNpts = 0;
  maxRate = 0;

  FILE *tecPlot, *tecPlot2;
  tecPlot = fopen("CNRwBGS.plt","w");
  fprintf( tecPlot, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE\n");

//  tecPlot2 = fopen("CNRwBGS2.plt","w");
//  fprintf( tecPlot2, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE\n");

  for(normal=0;normal<layerNNormal(layer);normal++){
    double sp1, sp2, z, r;
    double minSpacing, lengthOfMinSpacing;
    root = layerNormalRoot(layer, normal );
    gridNodeXYZ(layerGrid(layer),root,xyz);
    MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
    sp1 = spacing[0];
    z = -xyz[2];
    r = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    layerNormalDirection(layer,normal,normalDirection);
    length = layerNormalMaxLength(layer,normal);
    xyz0[0] = xyz[0] + (length * normalDirection[0]);
    xyz0[1] = xyz[1] + (length * normalDirection[1]);
    xyz0[2] = xyz[2] + (length * normalDirection[2]);
    MeshMgr_GetSpacing(&(xyz0[0]),&(xyz0[1]),&(xyz0[2]),spacing,direction);
    sp2 = spacing[0];
    

// sample the background spacing at several points
// in the direction of the normal

    minSpacing = MAX( sp1, sp2 );
    lengthOfMinSpacing = length;

/*
    int sample;
    for( sample=1; sample<10; sample++ ) {
       double len = sample*length/10.0;
       xyz0[0] = xyz[0] + (len * normalDirection[0]);
       xyz0[1] = xyz[1] + (len * normalDirection[1]);
       xyz0[2] = xyz[2] + (len * normalDirection[2]);
       MeshMgr_GetSpacing(&(xyz0[0]),&(xyz0[1]),&(xyz0[2]),spacing,direction);
	if( spacing[0] < minSpacing ) {
	   minSpacing = spacing[0];
	   lengthOfMinSpacing = len;
        }
    }
*/

    finalDelta = minSpacing*finalRatio;


    initialDelta = layerNormalInitialHeight(layer,normal);

// check some stuff

    if( lengthOfMinSpacing < initialDelta ) {
	printf("Node %d : layer thickness < dy-Min :  %f  %f ",
                normal, lengthOfMinSpacing, initialDelta );
	printf("at (%f, %f).     Resetting dy-Min\n", r, z );

        initialDelta = 0.99*lengthOfMinSpacing;
        layerSetNormalHeight( layer, normal, initialDelta );
    }

    if( finalDelta*0.75 < initialDelta ) {
/*
	printf(" node %d : dy-Max  <  0.75 * dy-Min :  %f   %f",
		 normal, finalDelta, initialDelta );
	printf("  at (%f, %f)  : TERMINATING NORMAL.\n",
		 r, z );
*/

// Setting height to 0.0 generates "CAPrI Projection errors???" errors
	initialDelta = finalDelta*0.75; //0.0;
        layerSetNormalHeight( layer, normal, initialDelta );
	layerTerminateNormal(layer,normal);
	rate = 1.001;
        layerSetNormalRate(layer, normal, rate );
        npts = 0;
    } else {

       if( lengthOfMinSpacing < finalDelta || initialDelta > finalDelta ) {
	   printf(" node %d : layer thickness  <  dy-Max :  %f  %f",
		    normal, lengthOfMinSpacing, finalDelta );
	   printf("  at (%f, %f)\n", r, z );

	   finalDelta = (2.0*lengthOfMinSpacing+initialDelta)/3.0;
       }

       rate = RateOfGeometricStretch( lengthOfMinSpacing, initialDelta, finalDelta );
       rate = MIN( 1.3, MAX( 1.01, rate ));
       npts = NptsOfGeometricStretch( lengthOfMinSpacing, initialDelta, rate );
       maxNpts = MAX( maxNpts, npts );
       maxRate = MAX( maxRate, rate );


    }

    layerSetNormalRate(layer, normal, rate);

//    if( z >= -1.0e-5 && r > 3.17501 )
       fprintf(tecPlot, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, lengthOfMinSpacing, initialDelta,
               finalDelta, rate, npts);
//    else
//       fprintf(tecPlot2, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
//               normal, xyz[0], xyz[1], z, r, lengthOfMinSpacing, initialDelta,
//               finalDelta, rate, npts);
  }

  fclose(tecPlot);
//  fclose(tecPlot2);
  printf(" Max rate and thickness layer : %f  %d \n", maxRate, maxNpts );
//  exit(0);
  return layer;
}

