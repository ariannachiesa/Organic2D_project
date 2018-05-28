/*! \file interp1.h
  \brief Interpolation function
*/

#include<vector>

double 
interp1( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate );