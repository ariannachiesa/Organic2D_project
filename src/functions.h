/*! \file functions.h
  \brief auxiliar functions for classes Newton and NLPoisson
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "probl.h"

#include<vector>

///	Computes the density of charge n and its derivative w.r.t the electric potential
///	by approximating them with Gauss-Hermite quadrature formula
void 
org_gaussian_charge_n( std::vector<double>& V,
						std::vector<double>& rhon, std::vector<double>& drhon_dV);
							
std::vector<double>
n_approx(  std::vector<double>& V, Probl& P);

std::vector<double> 
dn_dV_approx(  std::vector<double>& V, Probl& P);

/// Compute the (Inf,L2,H1) norm of a piecewise linear function.
void
bim2a_norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type);
					
#endif