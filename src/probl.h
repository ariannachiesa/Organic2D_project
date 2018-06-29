/*! \file probl.h
  \brief Class Problem
*/

#ifndef PROBL_H
#define PROBL_H

#include "bim_sparse.h"
#include "mumps_class.h"
#include "nonlinear_solver.h"
#include "octave_file_io.h"
#include "quad_operators.h"
#include "sandia_rules.hpp"
#include "tmesh.h"

#include<algorithm>
#include<array>
#include<fstream>
#include<vector>

/**
 * class Probl :
 * It groups all other parameters classes (Constants, Material, Quad, Algor, Device) 
 *	Among the input parameters for its constructor: number of refinement cycles for the mesh.
 *	N.B : if n.cycles = 0 then the mesh will be constitute by a single row of elements on the y-axis
 *	otherwise it will be refined n.cycles times and will have 2*n.cycles elements on both the x-axis and the y-axis
 */
class Probl
{
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
	double _eps_semic;				/**< permittivity of semiconductor */
	double _eps_ins;				/**< permittivity of insulator */
	double _PhiB;					/**< Metal to semiconductor barrier. [eV] */
	double _N0;						/**< Total number of available states (per unit volume) [m^{-3}] */
	double _Egap;					/**< LUMO-HOMO energy gap. [eV] */
	double _sigman;					/**< the disorder parameter [J] */
	double _sigman_kT;
	double _mu0n;					/**< low-field and low-density charge mobility for electrons [m^2 V^{-1} s^{-1}] */
	double _ni;
		
	/// Stores nodes and weigths of Gauss-Hermite quadrature rules
	/// class	Quad
	std::vector<double> _gx;			/**< vectors of nodes for quad. rules */
	std::vector<double> _gw;			/**< vectors of weigths for quad. rules */	
	
	/// Stores parameters for loops
	/// class	Algor
	int _pmaxit;
	double _ptoll;
	
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
	double _VG;									/**< [V] voltage applied at the gate terminal */
	double _VB;									/**< [V] voltage applied at the bulk terminal */
	double _Csb;								/**< [F] value of the capacitor of the external control circuit */
	double _t_semic;							/**< [m] thickness of semiconductor, along y-axis */
	double _t_ins;								/**< [m] thickness of insulator, along y-axis */
	double _L;									/**< [m] width of device, along x-axis */
	int _nTrees;								/**< number of Trees in the mesh */

	std::vector<int> _scnodes;					/**< vector of 1s and 0s; 1 = node in the semiconductor, 0 = node in the insulator */
	std::vector<int> _insulator;				/**< vector of 1s and 0s; 0 = element in the semiconductor, 1 = element in the insulator */
	std::array<int,2> _pins;					/**< [Bulk, gate] : region 1 = semiconductor , region 0 = insulator */
	std::array<int,2> _contacts;				/**< number of the geometrical border containing the side edges where the contacts are :
												*   edge 2 of tree 0, edge 3 of last tree
												*/
	std::vector< std::vector<int> > _dnodes;	/**< vector with gate nodes + vector with bulk nodes */
	
	tmesh	_msh;					/**< quadrangular mesh */
	
	std::vector<double> Vin;		/**< Output potential */
	std::vector<double> nin;		/**< Output electron density */
	std::vector<double> resnrm;		/**< residual norm */
	int niter;						/**< n. iterations required until convergence */
	
	Probl(	int maxcycle,
			double T0 = 300,																												// Constants
			double PhiB = 0.54, double sigman = 2.6, double mu0n = 4.29110133911508e-6,														// Material
			int nq = 101,																													// Quad
			int pmaxit = 1000, double ptoll = 1e-10,																						// Algor
			double Vshift = 1.79738, double Csb = 1.16183675549126e-11, double t_semic = 3.49436549222355e-8, double t_ins = 4.41e-7, 
			double L = 1.4e-3, bool ins = true,																								// Device
			std::array<int,2> pins = {1, 0}, std::array<int,2> contacts = {2, 3}, double section = 0.00000081, double Vdrain = 5, 
			double VG = 0.0, double VB = 0.0);			// constructor

    /// METHODS

	void Constants(double T0);
	void Material(double PhiB, double sigman, double mu0n);
	void Quad(int n);
	void Algor(	int pmaxit, double ptoll);				
	void Device(double Vshift, double Csb, double t_semic, double t_ins, double L, bool ins, 
				std::array<int,2>& pins, std::array<int,2>& contacts, double section, double Vdrain, double VG, double VB, int maxcycle);
				
	void Laplace();
	void LinearPoisson();
	void NonLinearPoisson(std::vector<double>& phi0);
	
	void org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& rhon, std::vector<double>& drhon_dV);
	void org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& drhon_dV);
							
	std::vector<double> n_approx(  std::vector<double>& V);
	std::vector<double> dn_dV_approx(  std::vector<double>& V);
	
	double CVcurve(std::vector<double>& phi);
	
	/// save methods
	void savePoisson(std::vector<double>& V, std::vector<double>& n, double niter, std::vector<double>& resnrm, const char* FileName);
	void saveCV(std::vector<double>& V, std::vector<double>& C, const char* FileName);

	/// Compute the (Inf,L2,H1) norm of a piecewise linear function.
	void Norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type);
	
	/// set methods
	void set_T0(double T0);

};

#endif /* PROBL_H */