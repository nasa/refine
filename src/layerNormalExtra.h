#ifndef EXTRAS_H
#define EXTRAS_H


Layer *layerAssignPolarGrowthHeight(Layer*, double, double, double* );
Layer *layerSetNormalHeightWithRateAcceleration(Layer*, double );

Layer *layerSmoothRate(Layer*, int, double, bool);
Layer *layerSmoothNormalProperties(Layer*, int*, double, bool);
Layer *layerSetNormalMaxLengthConstrained(Layer*, int, double );
Layer *layerSetNormalMaxDy(Layer*, int, double );
double layerNormalMaxDy(Layer*, int );
double layerNormalHeight(Layer*, int );

Layer *layerRelaxNormalDirection(Layer*, int, double);

void WriteTerminationMessage(Layer*,int,char*);
int layerTerminateNormalIfInX(Layer*, int, double);

//Layer *layerGetTriangleRates( Layer*, int, double* );
//int layerComputeNormalRateWithBGSpacing(Layer *layer, double finalRatio);

#endif
