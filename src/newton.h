/*! \file newton.h
  \brief Class Newton
*/

#ifndef NEWTON_H
#define NEWTON_H

#include "bcs_circuit.h"
#include "bim_sparse.h"
#include "probl.h"

class Newton{

	protected:
	sparse_matrix		_jac;
	std::vector<double>	_res;
	
	// std::vector< std::vector<double> >	_V;
	// std::vector< std::vector<double> >	_n;
	// std::vector< std::vector<double> >	_F;
	// std::vector< std::vector<double> >	_I;
	sparse_matrix	_V;
	sparse_matrix	_n;
	sparse_matrix	_F;
	sparse_matrix	_I;
	std::vector<double>	_tout;
	

	public:
        Newton() = delete; // default constructor
        Newton(	Probl& P,
				std::vector<double>& Vin, std::vector<double>& nin, std::vector<double>& tspan, std::vector<double>& Fin,
				std::vector<double>& Iin, BCS_CIRC& bcs, double freq);	// constructor

	
	///Assemble residual vector
	// std::vector<double>							
	// org_secs2d_newton_residual(	Probl& P,										
								// std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I,		
								// std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0,	
								// double deltat, BCS_CIRC& bcs, std::vector<int>& indexingV, std::vector<int>& indexingn, std::vector<int>& indexingF, 
								// std::vector<int>& indexingI);
	
	///Assemble jacobian matrix
	// void
	// org_secs2d_newton_jacobian(	Probl& P,
								// std::vector<double>& V, std::vector<double>& n, std::vector<double>& F,			
								// double deltat, BCS_CIRC& bcs, std::vector<int>& indexingV, std::vector<int>& indexingn,		
								// std::vector<int>& indexingF, std::vector<int>& indexingI, sparse_matrix& jacobian);
						
	void 
	compute_residual_norm (	double& resnrm, int& whichone, std::vector<double>& resall, std::vector<double>& res,
							std::vector<int>& idxV, std::vector<int>& idxn, std::vector<int>& idxF, std::vector<int>& idxI);
	
	// void
	// compute_variation (std::vector<double>& Va, std::vector<double>& na, std::vector<double>& Fa, std::vector<double>& Ia, 
                        // std::vector<double>& Vb, std::vector<double>& nb, std::vector<double>& Fb, std::vector<double>& Ib,
                        // Algor* alg, Device* dev, Material* mat, Constants* c, double clamping,
                        // double& incrV, double& incrn, double& incrF, double& incrI);
	
	void
	CONV_MSG (int tstep, int Nstep, int mNstep, double t, std::string reason, int field, double incr, std::vector<double>& res);
	
	void
	MAXINCR_MSG (int tstep, double t, int Nstep, int field, double incr, std::vector<double>& res, Probl& P);
	
	void
	DIV_MSG (	int tstep, double t, int Nstep, int field, std::vector<double>& incrhist, double incr, 
				std::vector<double>& res, int nsteps_check);

	void
	DIV_MN_MSG (int tstep, int t, int Nstep, int mNstep, int field, std::vector<double>& incrhist, double incr, 
				std::vector<double>& res, int nsteps_check);
				
	void
	org_secs_state_predict (Probl& P, sparse_matrix& V, sparse_matrix& n, sparse_matrix& F, sparse_matrix& I, int& tstep, std::vector<double>& tout,
							std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0);
	void
	diff (sparse_matrix& in, double a, double b, std::vector<double>& out);
	
	void
	org_secs_safe_increment (	std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0,
								std::vector<double>& dV, std::vector<double>& dn, std::vector<double>& dF, std::vector<double>& dI,
								Probl& P,
								std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I,
								int& clamp, double& tauk);
				
	/// It saves Newton output at each time step
	// void
	// saveNEWT(	std::vector<double>& Vold, std::vector<double>& nold, std::vector<double>& Fold, std::vector<double>& Iold, double told, 
				// std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I, std::vector<double>& res,
				// double t, double dt, int nsaves, int newton_solves, int modified_newton_solves, double freq);

};

#endif	// NEWTON_H