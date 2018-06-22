/*! \file newton12.h
  \brief Class Newton: 1st + 2nd eq. only
*/

#ifndef NEWTON_H
#define NEWTON_H

#include "bim_sparse.h"
#include "interp1.h"
#include "probl.h"

class Newton{

	protected:
	sparse_matrix		_jac;
	std::vector<double>	_res;

	std::vector<double>	_V;
	std::vector<double>	_n;
	std::vector<double>	_tout;
	

	public:
        Newton() = delete; // default constructor
        Newton(	Probl& P,
				std::vector<double>& Vin, std::vector<double>& nin, std::vector<double>& tspan, double freq);	// constructor

	
	///Assemble residual vector
	std::vector<double>							
	org_secs2d_newton_residual(	Probl& P,										
								std::vector<double>& V, std::vector<double>& n, 
								std::vector<double>& V0, std::vector<double>& n0,
								double deltat, ordering& ordV, ordering& ordn);
	
	///Assemble jacobian matrix
	void
	org_secs2d_newton_jacobian(	Probl& P, std::vector<double>& V, std::vector<double>& n, double deltat, 
								ordering& ordV, ordering& ordn,	sparse_matrix& jacobian);
						
	void 
	compute_residual_norm (	double& resnrm, int& whichone, std::vector<double>& resall, std::vector<double>& res,
							int nnodes, ordering& idxV, ordering& idxn);
	
	void
	compute_variation (	std::vector<double>& Va, std::vector<double>& na,
						std::vector<double>& Vb, std::vector<double>& nb,
						Probl& P, double clamping, double& incrV, double& incrn);
	
	void
	CONV_MSG (int tstep, int Nstep, int mNstep, double t, std::string reason, int field, double incr, std::vector<double>& res);
	
	void
	MAXINCR_MSG (int tstep, double t, int Nstep, int field, double incr, std::vector<double>& res, Probl& P);
	
	void
	DIV_MSG (	int tstep, double t, int Nstep, int field, std::vector<double>& incrhist, double incr, std::vector<double>& res, int nsteps_check);

	void
	DIV_MN_MSG (int tstep, double t, int Nstep, int mNstep, int field, std::vector<double>& incrhist, double incr, std::vector<double>& res, int nsteps_check);
							
	void
	org_secs_state_predict (Probl& P, std::vector<double>& Vold, std::vector<double>& nold,
							std::vector<double>& Voldold, std::vector<double>& noldold,
							int& tstep, std::vector<double>& tout, std::vector<double>& V0, std::vector<double>& n0);
	
	void
	org_secs_safe_increment (	std::vector<double>& V0, std::vector<double>& n0,
								std::vector<double>& dV, std::vector<double>& dn,
								Probl& P,
								std::vector<double>& V, std::vector<double>& n,
								double& clamp, double& tauk);
								
	double
	infnorm(std::vector<double>& in);
	
	bool
	any(std::vector<int>& v);

	bool
	any(std::vector<double>& v);
								
	void
	org_physical_models2d (std::vector<double>& n, Probl& P, std::vector<double>& mobilityn, std::vector<double>& alphan, std::vector<double>& der_dalpha_n);
	
	void
	org_physical_models2d (std::vector<double>& n, Probl& P, std::vector<double>& mobilityn, std::vector<double>& alphan);
	
	std::vector<double>
	org_secs_mobility_EGDM (double mu0, std::vector<double>& c, double C0, double sigma, Probl& P);
	
	void
	org_gaussian_charge_n(std::vector<double>& V, Probl& P,std::vector<double>& rhon, std::vector<double>& drhon_dV);
	
	void
	org_gaussian_charge_n(std::vector<double>& V, Probl& P,std::vector<double>& rhon);
	
	double
	org_gaussian_charge_n(double V, Probl& P);
	
	double
	n_approx(double V, Probl& P);
	
	std::vector<double> 
	n_approx(std::vector<double>& V, Probl& P);
	
	std::vector<double>
	dn_dV_approx(std::vector<double>& V, Probl& P);
	
	void
	org_einstein_n (std::vector<double>& n,Probl& P,std::vector<double>& alphan, std::vector<double>& dalphan_dn);
	
	std::vector<double>
	alphan_fun (std::vector<double>& phi, Probl& P);
	
	std::vector<double>
	dn_dphi_approx (std::vector<double>& phi, Probl& P);
	
	std::vector<double>
	d2n_dphi2_approx (std::vector<double>& phi, Probl& P);
	
	std::vector<double>
	dalphan_dn_fun (std::vector<double>& phi, std::vector<double>& alphan, Probl& P);
				
	/// It saves Newton output at the end of each time step
	void
	saveNEWT(std::vector<double>& Vold, std::vector<double>& nold, double told, 
			std::vector<double>& V, std::vector<double>& n, std::vector<double>& res,
			double t, double dt, int nsaves, int newton_solves, int modified_newton_solves, double freq);
			
	void
	saveJAC (int nrows, int ncols, std::vector<double>& vals);
	
	void
	saveVn(std::vector<double>& V, std::vector<double>& n, const char* FileName);
	
	void
	saveRES(std::vector<double>& res, const char* FileName);
};

#endif	// NEWTON_H