
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include <string.h>
#include "refine_defs.h"
#include "layer.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "grid.h"
#include "near.h"
#include "geometricStretch.h"
#include "layerNormalExtra.h"

//
// Functions:
//        layerAssignPolarGrowthHeight		cut from mesherx_SyJet.c
//        layerComputeNormalRateWithBGSpacing2  cut from layerNormalExtra.c
//
//  primary methods controlling spacing for SyJet
//        

// Set the desired layer thickness and height of first element
//
// Inside of jet,
//     make the first spacing a function of r that
//          constant for for r < rTrans0
//          grows like the sqrt(r) for rTrans0 < r < rTrans1
//          and then grows exponentially for r > rTrans1: dy0 = c^(r-rTrans1)
//          (must keep an eye on c relative to the spacing of surface mesh)
//
//     The layer thickness is C * sqrt(r) but bounded to half the height
//     of the jet cavity
//
// In the tunnel, the b.l. has a constant first spacing and total height
// based on input from experiment.
//

Layer *layerAssignPolarGrowthHeight(Layer *layer, 
                                         double scale,
					 double scalev,
					 double *referenceDirection)
{

  int normal;
  double cosAngle, arcAngle, distance, height, layerThickness;
  double normalDirection[3], xyzRoot[3];

// "r" dimensions of jet

  double rMin     = 3.175;
  double rMax     = 70.0/rMin;

// parameters describing (controlling) grid inside jet

  double rTrans0  = 2.0;
  double rTrans1  = 2.0;
  double pwrInner = 0.5;                 // growth for Laminar B.L.
  double hMinTunnel = 0.005*scalev;
  double hMinCavity = 0.003*scalev;
  double hMax     = 1.0*scale;
  int nLayer      = 70;
  double hLmax    = 3.0;

  double hTrans   = hMinCavity * pow( rTrans1/rTrans0, pwrInner );
  double sRate    = pow( hMax/hTrans, 1.0/(rMax-rTrans1) );

  hMax     = MIN ( hLmax, hMax );

  printf( " Parameters describing first spacing and layer thickness\n");
  printf( " rMin 	= %f\n", rMin);
  printf( " rMax 	= %f\n", rMax);
  printf( " rTrans1	= %f\n", rTrans1);
  printf( " pwrInner	= %f\n", pwrInner);
  printf( " hTrans	= %f\n", hTrans);
  printf( " hMinTunnel	= %f\n", hMinTunnel);
  printf( " hMinCavity  = %f\n", hMinCavity);
  printf( " hMax	= %f\n", hMax);
  printf( " sRate	= %f\n", sRate);
  
  for(normal=0;normal<layerNNormal(layer);normal++){
    double r;
    gridNodeXYZ(layerGrid(layer), layerNormalRoot(layer,normal), xyzRoot );
    layerNormalDirection(layer, normal, normalDirection);

    r = sqrt( xyzRoot[0]*xyzRoot[0] + xyzRoot[1]*xyzRoot[1] );

    if( xyzRoot[2] < -1.0e-5 || r < 3.175 + 1.0e-5) { // lip is treated as inside cavity
      double h1, xx;
          distance = sqrt( xyzRoot[0]*xyzRoot[0] +
		        xyzRoot[1]*xyzRoot[1] +
		        xyzRoot[2]*xyzRoot[2] )/rMin;
   

          if( distance < rTrans0 )  {
             height = hMinCavity;
          } else if( distance < rTrans1 )  {
             height = hMinCavity * pow( distance/rTrans0, pwrInner );
          } else {
	     xx = (distance-rTrans1);
	     height = hTrans * pow( sRate, xx );
          }

          layerThickness = LengthFromN( nLayer, height, 1.15 );

          layerThickness = MIN( 1.0 + 2.0*r/20.0, layerThickness );
          layerThickness = MIN( hLmax, layerThickness );

    } else {	// assume node is in the crossflow
		// use normal boundary layer profile
           layerThickness = 30.0;
           height = hMinTunnel;

    }

    height = MIN( height, 0.99*layerThickness );

    layerSetNormalMaxLength( layer, normal, layerThickness );
    layerSetNormalHeight( layer, normal, height);

  }
  return layer;	
}



Layer *layerComputeNormalRateWithBGSpacing2(Layer *layer, double finalRatio)
{
  int normal, root, maxNpts, npts;
  
  double xyz[3], xyz0[3], length, normalDirection[3];
  double initialDelta, finalDelta, rate, maxRate;
  double spacing[3], direction[9];
  FILE *tecPlot, *tecPlot2, *messages;

  maxNpts = 0;
  maxRate = 0;

  messages = fopen("CNRwGBS.msg","w");

  tecPlot  = fopen("CNRwBGS.plt","w");
  tecPlot2 = fopen("CNRwBGS2.plt","w");
  fprintf( tecPlot, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In tunnel after computing rate\"\n");
  fprintf( tecPlot2, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In cavity after computing rate\"\n");

  for(normal=0;normal<layerNNormal(layer);normal++){
    double sp1, sp2, z, r;
    double minSpacing, lengthOfMinSpacing;

    root = layerNormalRoot(layer, normal );
    gridNodeXYZ(layerGrid(layer),root,xyz);
    UG_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
    sp1 = spacing[0];
    z = xyz[2];
    r = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1]);
    layerNormalDirection(layer,normal,normalDirection);
    length = layerNormalMaxLength(layer,normal);
    xyz0[0] = xyz[0] + (length * normalDirection[0]);
    xyz0[1] = xyz[1] + (length * normalDirection[1]);
    xyz0[2] = xyz[2] + (length * normalDirection[2]);
    UG_GetSpacing(&(xyz0[0]),&(xyz0[1]),&(xyz0[2]),spacing,direction);
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
       UG_GetSpacing(&(xyz0[0]),&(xyz0[1]),&(xyz0[2]),spacing,direction);
	if( spacing[0] < minSpacing ) {
	   minSpacing = spacing[0];
	   lengthOfMinSpacing = len;
        }
    }
*/

    finalDelta = minSpacing*finalRatio;

// override and hardcode if in tunnel
//    if( length > 5.0 ) finalDelta = (finalDelta + 2.4)/2.0;

    initialDelta = layerNormalHeight(layer,normal);

// check some stuff

    if( lengthOfMinSpacing < initialDelta ) {
	fprintf(messages,"Node %d : layer thickness < dy-Min :  %f  %f ",
                normal, lengthOfMinSpacing, initialDelta );
	fprintf(messages,"at (%f, %f).     Resetting dy-Min\n", r, z );

        initialDelta = 0.99*lengthOfMinSpacing;
        layerSetNormalHeight( layer, normal, initialDelta );
    }

    if( finalDelta*0.75 < initialDelta ) {
	fprintf(messages," node %d : dy-Max  <  0.75 * dy-Min :  %f   %f",
		 normal, finalDelta, initialDelta );
	fprintf(messages,"  at (%f, %f)  : TERMINATING NORMAL.\n",
		 r, z );

// Setting height to 0.0 generates "CAPrI Projection errors???" errors
	initialDelta = finalDelta*0.75; //0.0;
        layerSetNormalHeight( layer, normal, initialDelta );
	layerTerminateNormal(layer,normal);
	rate = 1.001;
        layerSetNormalRate(layer, normal, rate );
        npts = 0;
    } else {

       if( lengthOfMinSpacing < finalDelta || initialDelta > finalDelta ) {
	   fprintf(messages," node %d : layer thickness  <  dy-Max :  %f  %f",
		    normal, lengthOfMinSpacing, finalDelta );
	   fprintf(messages,"  at (%f, %f)\n", r, z );

	   finalDelta = (2.0*lengthOfMinSpacing+initialDelta)/3.0;
       }

       rate = RateOfGeometricStretch( lengthOfMinSpacing, initialDelta, finalDelta );
       rate = MIN( 1.3, MAX( 1.15, rate ));
       npts = NptsOfGeometricStretch( lengthOfMinSpacing, initialDelta, rate );
       maxNpts = MAX( maxNpts, npts );
       maxRate = MAX( maxRate, rate );


    }

    layerSetNormalRate(layer, normal, rate);

    if( z >= -1.0e-5 && r > 3.17501 )
       fprintf(tecPlot, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, lengthOfMinSpacing, initialDelta,
               finalDelta, rate, npts);
    else
       fprintf(tecPlot2, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, lengthOfMinSpacing, initialDelta,
               finalDelta, rate, npts);
  }

  fclose(tecPlot);
  fclose(tecPlot2);

  printf(" Max rate and thickness layer : %f  %d \n", maxRate, maxNpts );
  return layer;
}

