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


double LengthFromN( int n, double dy0, double rate ) {
   return dy0*(pow(rate,n)-1.0)/(rate-1.0);
}

double LengthFromDyMax( double DyMax, double dy0, double rate ) {
   return (rate*DyMax - dy0)/(rate-1.0);
}

double RateOfGeometricStretch( double yMax, double dyMin, double dyMax )
{
       return ( yMax - dyMin ) / ( yMax - dyMax );
}

int NptsOfGeometricStretch( double yMax, double dyMin, double rate )
{
       if( rate <= 1.0001 ) return (int)(yMax/dyMin+0.5);
       else return (int)(0.5 + log( (rate-1.0)*yMax/dyMin + 1.0) / log(rate));
}

double RptsOfGeometricStretch( double yMax, double dyMin, double rate )
{
       if( rate <= 1.0001 ) return yMax/dyMin;
       else return log( (rate-1.0)*yMax/dyMin + 1.0) / log(rate) ;
}

double DyMaxFromN( int n, double rate, double dy0 ) {
	return dy0 * pow(rate,n);
}

double DyMaxFromYmax( double Ymax, double rate, double dy0 ) {
	return (Ymax*(rate - 1.0) + dy0 )/rate;
}

