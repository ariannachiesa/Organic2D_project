/*! \file functions.cpp
  \brief auxiliar functions for classes Newton and NLPoisson
*/

#include "functions.h"

void org_gaussian_charge_n( std::vector<double>& V,
							// output
							std::vector<double>& rhon, std::vector<double>& drhon_dV){

	double	q = _q;
    std::vector<double> n = n_approx(V,P);
	
	rhon = n;
    for(unsigned i=0; i<n.size(); i++){
		rhon[i] *= -q;
    }

    std::vector<double> dn_dV = dn_dV_approx(V,P);
	
	drhon_dV = dn_dV;
    for(unsigned i=0; i<dn_dV.size(); i++){
		drhon_dV[i] *= -q;
    }
};

std::vector<double> n_approx(  std::vector<double>& V, Probl& P){

    std::vector<double> coeff(V.size(),0),
						n(V.size(),0);
    double	q = P._q,
			kT = P._Kb*P._T0,
			sigman = P._sigman,
			N0 = P._N0,
			denom;

	std::vector<double>	gx = P._gx;
	std::vector<double>	gw = P._gw;	
	
    for(unsigned i=0; i<gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){		
            coeff[j] = (sqrt(2) * sigman * gx[i] - q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            n[j] += N0 / sqrt(M_PI) * gw[i] / denom;
        }
    }
	gx.clear();
	gw.clear();
	coeff.clear();	
	
	return n;
};

std::vector<double> dn_dV_approx(  std::vector<double>& V, Probl& P){

    std::vector<double> coeff(V.size(),0),
						dn_dV(V.size(),0);
    double	kT = P._Kb * P._T0,
			sigman = P._sigman,
			N0 = P._N0,
			q = P._q;
	double	denom;
	std::vector<double>	gx = P._gx;
	std::vector<double>	gw = P._gw;

    for(unsigned i=0; i<gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){
            coeff[j] = (sqrt(2) * sigman * gx[i] - q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            dn_dV[j] += - q * N0 / sigman * sqrt(2/M_PI) * gw[i]*gx[i] / denom;
		}
    }
	gx.clear();
	gw.clear();
	coeff.clear();
	return dn_dV;
};

void
bim2a_norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type)
{
  if (type == Inf)
    {
      norm = 0.0;
      for (unsigned int i = 0; i < v.size (); ++i)
        {
          double temp = std::fabs (v[i]);
          if (norm < temp)
            norm = temp;
        }
    }
  else if (type == L2 || type == H1)
    {
      norm = 0.0;
      std::vector<double> ecoeff (msh.num_local_quadrants (), 1.0);
      std::vector<double> ncoeff (msh.num_local_nodes (), 1.0);
	  std::vector<double> psi (msh.num_local_nodes (),0);
      sparse_matrix M;
	  M.resize(msh.num_local_nodes ());
      bim2a_reaction (msh, ecoeff, ncoeff, M);
      if (type == H1)
        bim2a_advection_diffusion (msh, ecoeff, psi, M);
      std::vector<double> temp;
      temp = M * v;
      for (unsigned int i = 0; i < v.size (); ++i)
        norm += v[i] * temp[i];
      norm = sqrt (norm);
    }
};