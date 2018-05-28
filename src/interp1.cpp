/*! \file interp1.cpp
  \brief Interpolation function
*/

#include "interp1.h"

double 
interp1( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
{
   bool	found = false;
   int	size = xData.size(),
		i = 0,
		begin = 0,
		end = xData.size() - 1;
		
   if ( x >= xData[size - 2] )                                                 	// find left end of interval for interpolation
   {																			// special case: beyond right end 
      i = size - 2;
   }
   else
   {
	  while( begin <= end && found == false  ){
		
		i = ( begin + end ) / 2;

		if( x == xData[i] ){
			found = true;
		}
		else{
			if( x < xData[i] ){
				end = i-1;
			}
			else{
				begin = i+1;
			}
		}
	}
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends) 
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating 
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }
   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient 
   return yL + dydx * ( x - xL );                                              // linear interpolation
}