/*! \file probl.h
  \brief Class Problem
*/

#ifndef PROBL_H
#define PROBL_H

#include "bim_sparse.h"
//#include "functions.h"
#include "mumps_class.h"
#include "quad_operators.h"
#include "sandia_rules.hpp"
#include "tmesh.h"

#include<array>
#include<vector>

/**
 * class Probl :
 * It groups all other parameters classes (Constants, Material, Quad, Algor, Device) 
 *	as it contains pointers to them + stores an interpolation table needed when performing Newton's algorithm
 *	Input parameter for its constructor: number of refinement cycles for the mesh.
 *	N.B : if n.cycles = 0 then the mesh will be constitute by a single row of elements on the y-axis
 *	otherwise it will be refined n.cycles times and will have 2*n.cycles elements on both the x-axis and the y-axis
 */
class Probl
{	

	protected:
		
	std::vector<double>	_data_phi_lumo;		/**< interpolation table */
	std::vector<double>	_data_n;			/**< interpolated values */	
		
    public:
	
	/// Class to store physical constants
	/// class	Constants
	double _T0;				/**< absolute temperature */
	double _Kb;				/**< Boltzmann's constant */
	double _q;				/**< electron charge */
	double _eps0;			/**< vacuum permittivity */
	double _Vth;			/**< threshold voltage */
	
	/// Class Material: stores material characteristic prameters
	/// class	Material
	double _eps_semic_r;			/**< relative permittivity of semiconductor */
	double _eps_ins_r;				/**< relative permittivity of insulator */
	double _eps_semic;
	double _eps_ins;
	double _PhiB;					/**< Metal to semiconductor barrier. [eV] */
	double _N0;						/**< Total number of available states (per unit volume) [m^{-3}] */
	double _Egap;					/**< LUMO-HOMO energy gap. [eV] */
	double _sigman;					/**< [J]	the disorder parameter */
	double _sigman_kT;
	double _mu0n;					/**< [m^2 V^{-1} s^{-1}] low-field and low-density charge mobility for electrons */
	double _ni;
		
	/// Stores nodes and weigths of Gauss-Hermite quadrature rules
	/// class	Quad
	std::vector<double> _gx;			/**< vectors of nodes for quad. rules */
	std::vector<double> _gw;			/**< vectors of weigths for quad. rules */	
	
	/// Stores parameters for loops
	/// class	Algor
	int _pmaxit;
	int _maxit;
	int _maxit_mnewton;
	int _nsteps_check;
	double _maxnpincr;
	double _ptoll;
	double _toll;
	double _dt0;
	double _dtcut;
	double _dtmax;
	double _dtmin;
	double _maxdtincr;

	bool _clampOnOff;
	bool _savedata;
	
	std::vector<double> _colscaling;
	std::vector<double> _rowscaling;
	std::vector<int> _clamping;
	
	///	Class which stores the geomterical settings of the device,
	///	generate the mesh and assemble vectors useful to understand
	///	which nodes are in the semiconductor and which in the insulator,
	///	which elements are in the semiconductor and which in the insulator
	///	and which are the boundary nodes
	/// class	Device	
	bool _ins;									/**< true = (default) if insulator is present , false = if there's only semiconductor */
	
	double _section;							/**< device section [m^2] */
	double _Efield;								/**< device source-drain electric field [V / m] */
	
	double _Vshift;								/**< [V] shift potential sensed at the gate terminal */
	double _Csb;								/**< [F] value of the capacitor of the external control circuit */
	double _t_semic;							/**< [m] thickness of semiconductor, along y-axis */
	double _t_ins;
	double _L;

	std::vector<int> _scnodes;					/**< vector of 1s and 0s; 1 = node in the semiconductor, 0 = node in the insulator */
	std::vector<int> _insulator;				/**< vector of 1s and 0s; 0 = element in the semiconductor, 1 = element in the insulator */
	std::array<int,2> _pins;					/**< [Bulk, gate] : region 1 = semiconductor , region 0 = insulator */
	std::array<int,2> _contacts;				/**< number of the geometrical border containing the side edges where the contacts are :
												*   edge 0 of tree 0, edge 1 of tree 1
												*/
	std::vector<int> _alldnodes;				/**< vector of all boundary nodes on bulk and gate */
	std::vector< std::vector<int> > _dnodes;	/**< vector with gate nodes + vector with bulk nodes */
	
	tmesh	_msh;
	
	std::vector<double> Vin;		/**< Initial guess for the potential */
	std::vector<double> nin;		/**< Initial guess for the electron density */
	std::vector<double> res;		/**< residual vector */
	int niter;						/**< n. iterations required to compute the initial guess */
		
	Probl(	int maxcycle,
			double T0 = 300,																												// Constants
			double PhiB = 0.54, double sigman = 2.6, double mu0n = 4.29110133911508e-6,														// Material
			int nq = 101,																													// Quad
			int pmaxit = 1000, int maxit = 5, int maxit_mnewton = 30, int nsteps_check = 3, double maxnpincr = 1e-3, double ptoll = 1e-10, 
			double toll = 1e-4, double dt0 = 1e-10, double dtcut = 0.25, double dtmax = 1, double dtmin = 1e-12, double maxdtincr = 2,		// Algor
			double Vshift = 1.79738, double Csb = 1.16183675549126e-11, double t_semic = 3.49436549222355e-8, double t_ins = 4.41e-7, 
			double L = 1e-5, bool ins = true,																								// Device
			std::array<int,2> pins = {1, 0}, std::array<int,2> contacts = {0, 1}, double section = 0.00000081, double Vdrain = 5);
			// constructor

    /// METHODS

	void Constants(double T0);
	void Material(double PhiB, double sigman, double mu0n);
	void Quad(int n);
	void Algor(	int pmaxit, int maxit, int maxit_mnewton, int nsteps_check, double maxnpincr, double ptoll, 
				double toll, double dt0, double dtcut, double dtmax, double dtmin, double maxdtincr);
	void Device(double Vshift, double Csb, double t_semic, double t_ins, double L, bool ins, 
				std::array<int,2>& pins, std::array<int,2>& contacts, double section, double Vdrain, int maxcycle);
				
	void LinearPoisson(	std::vector<double>& phi0);
	
	/// It saves Poisson output
	void savePoisson(std::vector<double>& V, std::vector<double>& n, double niter, std::vector<double>& resnrm, const char* FileName);
	
	// void org_gaussian_charge_n( std::vector<double>& V, std::vector<double>& rhon, std::vector<double>& drhon_dV);
							
	// std::vector<double> n_approx(  std::vector<double>& V);
	// std::vector<double> dn_dV_approx(  std::vector<double>& V);

	/// Compute the (Inf,L2,H1) norm of a piecewise linear function.
	void bim2a_norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type);
	
	/// get methods
	std::vector<double>& get_data_phi_lumo();
	std::vector<double> get_data_n();
	
	/// set methods
	void set_T0(double T0);
	
	void set_PhiB(double PhiB);
	void set_sigman(double sigman);
	void set_sigmankT(double sigmankT);
	void set_mu0n(double mu0n);
	
	void set_pmaxit(int pmaxit);
	void set_maxit(int maxit);
	void set_maxit_mnewton(int maxit_mnewton);
	void set_nsteps_check(int nsteps_check);
	void set_maxnpincr(double maxnpincr);
	void set_ptoll(double ptoll);
	void set_toll(double toll);
	void set_dt0(double dt0);
	void set_dtcut(double dtcut);
	void set_dtmax(double dtmax);
	void set_dtmin(double dtmin);
	void set_maxdtincr(double maxdtincr);
	
	void set_Vshift(double Vshift);
	void set_Csb(double Csb);
	void set_Vdrain(double Vdrain);
	void set_section(double section);
	
};

	// /// class NLPoisson
	// class NLPoisson : public Probl
	// {
		// std::vector<double> _Vin;		/**< Initial guess for the potential */
		// std::vector<double> _nin;		/**< Initial guess for the electron density */
		// std::vector<double> _res;		/**< residual vector */
		// int _niter;						/**< n. iterations required to compute the initial guess */
		
		// public:
		// NLPoisson() = delete;
		// NLPoisson(	std::vector<double>& phi0,
					// int maxcycle,
			// double T0 = 300,	// Constants
			// double PhiB = 0.54, double sigman = 2.6, double mu0n = 4.29110133911508e-6,		// Material
			// int nq = 101,		// Quad
			// int pmaxit = 1000, int maxit = 5, int maxit_mnewton = 30, int nsteps_check = 3, double maxnpincr = 1e-3, double ptoll = 1e-10, 
			// double toll = 1e-4, double dt0 = 1e-10, double dtcut = 0.25, double dtmax = 1, double dtmin = 1e-12, double maxdtincr = 2,	// Algor
			// double Vshift = 1.79738, double Csb = 1.16183675549126e-11, double t_semic = 3.49436549222355e-8, double t_ins = 4.41e-7, 
			// double L = 1e-5, bool ins = true,		// Device
			// std::array<int,2> pins = {1, 0}, std::array<int,2> contacts = {0, 1}, double section = 0.00000081, double Vdrain = 5);//: Probl(int maxcycle);
			// // constructor
		
		// // ///	Computes the density of charge n and its derivative w.r.t the electric potential
		// // ///	by approximating them with Gauss-Hermite quadrature formula
		// // void 
		// // org_gaussian_charge_n( std::vector<double>& V,
								// // std::vector<double>& rhon, std::vector<double>& drhon_dV);
							
		// // std::vector<double>
		// // n_approx(  std::vector<double>& V, Probl& P);

		// // std::vector<double> 
		// // dn_dV_approx(  std::vector<double>& V, Probl& P);

		// // /// Compute the (Inf,L2,H1) norm of a piecewise linear function.
		// // void
		// // bim2a_norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type);
	// };

#endif /* PROBL_H */