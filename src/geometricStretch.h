//
//  some general helper functions
//
//
// functions related to geometric progressions
//
//  dy(n+1) = rate * dy(n);
//
//  y(n) = sum dy(m), m=0..n   for n < N
//
//  5 general parameters:
//
//      dy(0)	->	the first spacing
//	dy(N)	->	the last and max spacing
//	y(N)	->	max y
//	N	->	number of points
//	rate	->	growth rate
//
// given three out of 5 of the basic properties: return one of the other two

#ifndef GEOMETRICSTRETCH_H
#define GEOMETRICSTRETCH_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

double LengthFromN( int n, double dy0, double rate );

double LengthFromDyMax( double DyMax, double dy0, double rate );

double RateOfGeometricStretch( double yMax, double dyMin, double dyMax );

int NptsOfGeometricStretch( double yMax, double dyMin, double rate );

double RptsOfGeometricStretch( double yMax, double dyMin, double rate );

double DyMaxFromN( int n, double rate, double dy0 );

double DyMaxFromYmax( double Ymax, double rate, double dy0 );

END_C_DECLORATION

#endif
