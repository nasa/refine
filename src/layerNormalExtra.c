
// taken from layer.c

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include <string.h>
#include "master_header.h"
#include "layer.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "grid.h"
#include "near.h"
#include "geometricStretch.h"
#include "layerNormalExtra.h"

#include "layerStruct.h"


// adjust the height (and rate) between each layer.

// compiled and linked so as to overide the version in mesherx.c
// (had to comment version in mesherx.c)

Layer *layerSetNormalHeightWithRateAcceleration(Layer *layer, double maxRate)
{
  int normal;
  double height;

  for(normal=0;normal<layerNNormal(layer);normal++){

    double rate = layerNormalRate(layer, normal);

// accelerate the rate towards maxRate
// "fast" normals are unchanged

    if( maxRate > rate ) {
       rate = MIN( maxRate, 1.02*rate );
       layerSetNormalRate( layer, normal, rate );
    }

    height = layerNormalHeight( layer, normal );
    layerSetNormalHeight( layer, normal, height*rate);
  }

  return layer;
}


Layer *layerSmoothRate(Layer *layer, int itMax, double omega, bool iprt)
{
  int normal, iter, triangle, normals[3], total, i;
  double norm[3], avgdir[3], rate, denom;
  AdjIterator it;
  int minTriangle, lastTriangle;

  FILE *tecPlot, *tecPlot2;
  if( iprt ) {
    tecPlot = fopen("SmooRt.plt","w");
    tecPlot2 = fopen("SmooRt2.plt","w");
    fprintf( tecPlot, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In tunnel after smoothing rate\"\n");
    fprintf( tecPlot2, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In cavity after smoothing rate\"\n");
  }


  if (layerNNormal(layer) == 0 ) return NULL;

  for (iter=0;iter<itMax;iter++){

    for (normal=0;normal<layerNNormal(layer);normal++){
      double aveRate, rr;
      double rate = layerNormalRate(layer,normal);

        total = 0;
	aveRate = 0.0;
        for ( it = adjFirst(layer->adj,normal);
              adjValid(it);
              it = adjNext(it) ){
          triangle = adjItem(it);
          layerTriangleNormals(layer,triangle,normals);
          for (i=0;i<3;i++){
            if ( !layerNormalTerminated( layer,normals[i]) ){
	      rr = layerNormalRate(layer,normals[i]);
              aveRate += rr;
              total++;
            }
          }
        }
	if( total > 0 ) {
          aveRate /= (double)total;
	  rate += omega*(aveRate-rate);
          layerSetNormalRate( layer, normal, rate );
	}

// write a tecplot file on last iteration

      if( iter == itMax-1 && !layerNormalTerminated( layer,normal) && iprt ) {
	double sp1, z, r, length, height;
	int npts;
         int root = layerNormalRoot(layer, normal );
         double xyz[3], spacing[3], direction[9];
         gridNodeXYZ(layerGrid(layer),root,xyz);
         MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
         sp1 = spacing[0];
         z = xyz[2];
         r = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1]);
         length = layerNormalMaxLength(layer,normal);
         height = layerNormalHeight(layer,normal);
         npts = NptsOfGeometricStretch( length, height, rate );

         if( z >= -1.0e-5 && r > 3.17501 )
            fprintf(tecPlot, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, length, height,
               sp1, rate, npts);
         else
            fprintf(tecPlot2, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, length, height,
               sp1, rate, npts);
      }
    }
  }

  if( iprt ) {
    fclose(tecPlot);
    fclose(tecPlot2);
  }
  return layer;
}


// A generalized smoothing routine that
// will smooth any of several properties
// Initially smooths rate, height, and length
// with plans to add DyMax at a later date.
//
// NOTE: the 4 properties are not independent,
// so only 3 of the 4 should be smoothed at any given time.
// Rate, length and height are stiffly coupled to each other while
// height can be halfed or doubled with little effect on the others.
//
// DyMax is not directly stored in the normal but
// is computed as needed from rate, length, and height.
// If DyMax is smoothed, then either rate or length should be over ridden.
// If the user inadvertently request to smooth all 4 properties,
// then length is over ridden, but is constrained to be less than a
// specified value (intended to avoid colisions with geometry).
// For now, the constraint on length is 1.5 * its initial user specified
// value which is temporally stored in "normal.initialheight".

Layer *layerSmoothNormalProperty(Layer *layer, int itMax[4], double omega, bool iprt)
{
  int normal, iter, triangle, normals[3], i, ItMax;
  int totalR, totalH, totalL, totalY;
  double norm[3], avgdir[3], rate, length, height, denom;
  AdjIterator it;
  int minTriangle, lastTriangle;

  FILE *tecPlot, *tecPlot2;
  if( iprt ) {
    tecPlot = fopen("SmooNp.plt","w");
    tecPlot2 = fopen("SmooNp2.plt","w");
    fprintf( tecPlot , " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In tunnel after smoothing h and L\"\n");
    fprintf( tecPlot2, " VARIABLES = \"n\" \"x\" \"y\" \"z\" \"r\" \"L\" \"dY_0\" \"dY_M_a_x\" \"rate\" \"npts\" \n ZONE t=\"In cavity after smoothing h and L\"\n");
  }


  if (layerNNormal(layer) == 0 ) return NULL;

  ItMax = MAX( itMax[0], MAX( itMax[1], itMax[2]));
  for (iter=0;iter<ItMax;iter++){

    for (normal=0;normal<layerNNormal(layer);normal++){
      double aveRate   = 0.0;
      double aveHeight = 0.0;
      double aveLength = 0.0;
      double aveDyMax  = 0.0;

// temporally store a gross upper bound on length in the "initialHeight"

      length = layerNormalMaxLength(layer,normal);
//      height = layerNormalHeight(layer,normal);
      if( iter == 0 ) layerSetNormalInitialHeight(layer,normal, length);

        totalR = 0;
        totalH = 0;
        totalL = 0;
	totalY = 0;
        for ( it = adjFirst(layer->adj,normal);
              adjValid(it);
              it = adjNext(it) ){
          triangle = adjItem(it);
          layerTriangleNormals(layer,triangle,normals);
          for (i=0;i<3;i++){
            if ( !layerNormalTerminated( layer,normals[i]) ){

	      if( iter < itMax[0] ) {
                 aveRate += layerNormalRate(layer,normals[i]);
                 totalR++;
	      }
	      if( iter < itMax[1] ) {
                 aveHeight+= layerNormalHeight(layer,normals[i]);
                 totalH++;
	      }
	      if( iter < itMax[2] ) {
                 aveLength += layerNormalMaxLength(layer,normals[i]);
                 totalL++;
	      }
/*
	      if( iter < itMax[3] ) {
                 aveDyMax += layerNormalMaxDy(layer,normals[i]);
                 totalY++;
	      }
*/

            }
          }
        }
	if( totalR > 0 ) {
          aveRate /= (double)totalR;
          rate = layerNormalRate(layer,normal);
	  rate += omega*(aveRate-rate);
          layerSetNormalRate( layer, normal, rate );
	}
	if( totalH > 0 ) {
          aveHeight /= (double)totalH;
          height = layerNormalHeight(layer,normal);
	  height += omega*(aveHeight-height);
          layerSetNormalHeight( layer, normal, height);
	}
	if( totalL > 0 ) {
          aveLength /= (double)totalL;
          length = layerNormalMaxLength(layer,normal);
	  length += omega*(aveLength-length);
          layerSetNormalMaxLengthConstrained( layer, normal, length );
	}
/*
	if( totalY > 0 ) {
          aveDyMax /= (double)totalY;
          DyMax = layerNormalMaxDy(layer,normal);
	  DyMax += omega*(aveDyMax-DyMax);
          layerSetNormalMaxDy( layer, normal, DyMax );
	}
*/

// write a tecplot file on last iteration

      if( iter == ItMax-1 && !layerNormalTerminated( layer,normal) && iprt ) {
	double sp1, z, r;
	int npts;
         int root = layerNormalRoot(layer, normal );
         double xyz[3], spacing[3], direction[9];
         gridNodeXYZ(layerGrid(layer),root,xyz);
         MeshMgr_GetSpacing(&(xyz[0]),&(xyz[1]),&(xyz[2]),spacing,direction);
         sp1 = spacing[0];
         z = xyz[2];
         r = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	 rate = layerNormalRate(layer,normal);
         length = layerNormalMaxLength(layer,normal);
         height = layerNormalHeight(layer,normal);
         npts = NptsOfGeometricStretch( length, height, rate );
	 layerSetNormalInitialHeight( layer, normal, height );

         if( z >= -1.0e-5 && r > 3.17501 )
            fprintf(tecPlot, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, length, height,
               sp1, rate, npts);
         else
            fprintf(tecPlot2, " %d  %f  %f  %f  %f  %f  %f  %f  %f  %d \n",
               normal, xyz[0], xyz[1], z, r, length, height,
               sp1, rate, npts);
      }
    }
  }

  if( iprt ) {
    fclose(tecPlot);
    fclose(tecPlot2);
  }
  return layer;
}

double layerNormalHeight( Layer *layer, int i ) {
   return layer->normal[i].height;
}

double layerNormalMaxDy( Layer *layer, int i ) {

   double height = layer->normal[i].height;
   double npts = RptsOfGeometricStretch(layer->normal[i].length,
					height,
					layer->normal[i].rate );
   return pow(height,npts);

}

Layer *layerSetNormalMaxDy( Layer *layer, int normal, double DyMax ) {
   double height = layer->normal[normal].height;
   double rate   = layer->normal[normal].rate;
   layerSetNormalMaxLength(layer,normal,  (rate*DyMax-height)/(rate-1.0) );
}

// note: in the following function, initialHeight has been loaded
// with a gross upper bound on the length

Layer *layerSetNormalMaxLengthConstrained( Layer *layer, int normal, double length ) {
      double ih = layer->normal[normal].initialheight;
      layerSetNormalMaxLength(layer,normal, MIN( length, ih ) );
}

// function to find out why some points are terminating

void WriteTerminationMessage(Layer* layer, int normal, char *message ) {
  double z, r;
     int root = layerNormalRoot(layer, normal );
     double xyz[3], spacing[3], direction[9];
     gridNodeXYZ(layerGrid(layer),root,xyz);
     z = xyz[2];
     r = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1]);
//     if( z >= -1.0e-5 && r > 3.17501 && r < 3.3 ) {
	printf(" terminating normal %d because %s\n", normal, message );
//    }
}

Layer *layerRelaxNormalDirection(Layer *layer, int itMax, double omega )
{
  int normal, iter, triangle, normals[3], total, i;
  double norm[3], avgdir[3], denom;
  AdjIterator it;
  int minTriangle, lastTriangle;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (layerNBlend(layer) != 0 ) return NULL;

  layerProjectNormalsToConstraints(layer);

  for (iter=0; iter<itMax; iter++){
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 < layerConstrained(layer,normal) ){
        total = 0;
        avgdir[0]=0.0;
        avgdir[1]=0.0;
        avgdir[2]=0.0;
        for ( it = adjFirst(layer->adj,normal);
              adjValid(it);
              it = adjNext(it) ){
          triangle = adjItem(it);
          layerTriangleNormals(layer,triangle,normals);
          for (i=0;i<3;i++){
            if (normal != normals[i] && 0 != layerConstrained(layer,normal) ){
              layerNormalDirection(layer,normals[i],norm);
              avgdir[0] += norm[0];
              avgdir[1] += norm[1];
              avgdir[2] += norm[2];
              total++;
            }
          }
        }
	if( total > 0 ) {
           denom = 1.0 / (double)total;
	   avgdir[0] *= denom;
	   avgdir[1] *= denom;
	   avgdir[2] *= denom;
	   avgdir[0] -= layer->normal[normal].direction[0];
	   avgdir[1] -= layer->normal[normal].direction[1];
	   avgdir[2] -= layer->normal[normal].direction[2];
           layer->normal[normal].direction[0] += avgdir[0] * omega;
           layer->normal[normal].direction[1] += avgdir[1] * omega;
           layer->normal[normal].direction[2] += avgdir[2] * omega;
	}
      }
    }
    layerVisibleNormals(layer,0.5,1.0e-5);
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 == layerConstrained(layer,normal) ){
        total = 0;
        avgdir[0]=0.0;
        avgdir[1]=0.0;
        avgdir[2]=0.0;
        for ( it = adjFirst(layer->adj,normal);
              adjValid(it);
              it = adjNext(it) ){
          triangle = adjItem(it);
          layerTriangleNormals(layer,triangle,normals);
          for (i=0;i<3;i++){
            if (normal != normals[i]){
              layerNormalDirection(layer,normals[i],norm);
              avgdir[0] += norm[0];
              avgdir[1] += norm[1];
              avgdir[2] += norm[2];
              total++;
            }
          }
        }
	if( total > 0 ) {
           denom = 1.0 / (double)total;
	   avgdir[0] *= denom;
	   avgdir[1] *= denom;
	   avgdir[2] *= denom;
	   avgdir[0] -= layer->normal[normal].direction[0];
	   avgdir[1] -= layer->normal[normal].direction[1];
	   avgdir[2] -= layer->normal[normal].direction[2];
           layer->normal[normal].direction[0] += avgdir[0] * omega;
           layer->normal[normal].direction[1] += avgdir[1] * omega;
           layer->normal[normal].direction[2] += avgdir[2] * omega;
	}
      }
    }
    layerVisibleNormals(layer,0.5,1.0e-5);
  }
  return layer;
}


int layerTerminateNormalIfInX(Layer *layer, int s, double x)
{
  double xyz[3];
  int normal;
  int iflag = 1;
  int totalterm;

// terminate all remaining normals if any normal is a region has
// been terminated

  if (layerNNormal(layer) == 0 ) return 0;


  for (normal=0;normal<layerNNormal(layer)&&iflag;normal++){
    gridNodeXYZ(layer->grid, layer->normal[normal].root, xyz);
    if ( s*xyz[0] > s*x && layerNormalTerminated(layer,normal) ) iflag = 0;
  }

  if( !iflag ) {
     for (normal=0;normal<layerNNormal(layer);normal++)
	layerTerminateNormal(layer,normal);
     totalterm = layerNNormal(layer)-layerNActiveNormal(layer);
     printf("%d of %d normals terminted by region criterion.\n",
         totalterm,layerNNormal(layer) );
     return totalterm;
  }
  return 0;

}
