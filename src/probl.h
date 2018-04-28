/*! \file probl.h
  \brief Class Problem
*/

#ifndef PROBL_H
#define PROBL_H

#include "sandia_rules.hpp"
#include "tmesh.h"

#include<array>
#include<vector>

/**
 * It groups all other parameters classes (Constants, Material, Quad, Algor, Device) 
 *	as it contains pointers to them + stores an interpolation table needed when performing Newton's algorithm
 */
class Probl
{	
    public:
	
		/// Class to store physical constants
		class	Constants;
		
		/// Class Material: stores material characteristic prameters
		class	Material{
			    protected:
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
	
				public:
				Material() = delete; // default constructor
				Material(Constants c, double PhiB, double sigman, double mu0n); // constructor
				
				friend class Probl;
	
				/*double get_PhiB();
				double get_eps_ins();
				double get_eps_ins_r();
				double get_eps_semic();
				double get_sigman();
				double get_sigman_kT();
				double get_N0();
				double get_mu0n();
				double get_ni();
				void assign_ni(double ni);
				void assign_eps_ins(double x);
				void assign_eps_ins_r(double x);*/
		};
		
		class	Constants{
			protected:
			double _Kb;				/**< Boltzmann's constant */
			double _q;				/**< electron charge */
			double _eps0;			/**< vacuum permittivity */
			double _Vth;			/**< threshold voltage */
			double _T0;				/**< absolute temperature */
	
			public:
			Constants() = delete;
			Constants(double T0); // constructor

			/*double get_Kb();
			double get_q();
			double get_eps0();
			double get_T0();
			double get_Vth();*/
	
			friend class Probl::Material;
			friend class Probl;
		};
		
		/// Stores nodes and weigths of Gauss-Hermite quadrature rules
		class	Quad{
			protected:
			std::vector<double> _gx;			/**< vectors of nodes for quad. rules */
			std::vector<double> _gw;			/**< vectors of weigths for quad. rules */
	
			public:
			Quad() = delete; 	// default constructor
			Quad(int n); 		// constructor

			/*std::vector<double> get_gx();
			std::vector<double> get_gw();
	
			double* get_pgx();
			std::vector<double>* get_pgw();*/
			
			friend class Probl;
		};
		
		/// Stores parameters for loops
		class	Algor{
			protected:
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
	
			public:
			Algor() = delete; // default constructor
			Algor(	int pmaxit, int maxit, int maxit_mnewton, int nsteps_check, double maxnpincr, double ptoll, 
					double toll, double dt0, double dtcut, double dtmax, double dtmin, double maxdtincr); // constructor
	
			/*int get_pmaxit();
			int get_maxit();
			int get_maxit_mnewton();
			double get_ptoll();
			double get_toll();
			double get_dt0();
			double get_dtmax();
			double get_dtmin();
			double get_dtcut();
			double get_maxnpincr();
			double get_maxdtincr();
			int get_nsteps_check();
			bool get_plotsOff();
			bool get_clampOnOff();
			bool get_savedata();
			std::vector<double> get_rowscaling();
			std::vector<double> get_colscaling();
			std::vector<int> get_clamping();
			void assign_row(std::vector<double>& r);*/
			
			friend class Probl;
		};
		
		
		///	Class which stores the geomterical settings of the device,
		///	generate the mesh and assemble vectors useful to understand
		///	which nodes are in the semiconductor and which in the insulator,
		///	which elements are in the semiconductor and which in the insulator
		///	and which are the boundary nodes
		class	Device{
			protected:
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
	
			public:
			Device() = delete; // default constructor
			Device(	double Vshift, double Csb, double t_semic, double t_ins, double L, bool ins, 
					std::array<int,2>& pins, std::array<int,2>& contacts, double section, double Vdrain, int maxcycle); // constructor
	
			bool get_ins();
			double get_section();
			double get_Efield();
			double get_Vshift();
			double get_t_semic();
			double get_t_ins();
			double get_L();
			double get_Csb();
			void assign_Csb(double x);

			double get_mshp(tmesh::idx_t i, tmesh::idx_t j);
	
			tmesh* get_msh();
			tmesh::idx_t get_msh_gnodes();
			tmesh::idx_t get_msh_lnodes();
			tmesh::idx_t get_msh_gelem();
			tmesh::idx_t get_msh_lelem();
	
			std::array<int,2> get_pins();

			std::vector<int> get_alldnodes();
			std::vector<int> get_insulator();
			std::vector<int> get_scnodes();
			std::vector< std::vector<int> > get_dnodes();
			
			friend class Probl;
		};
		
        //Probl() = delete;				
		Probl(	int maxcycle,
				double T0 = 300,	// Constants
				double PhiB = 0.54, double sigman = 2.6, double mu0n = 4.29110133911508e-6,		// Material
				int nq = 101,		// Quad
				int pmaxit = 1000, int maxit = 5, int maxit_mnewton = 30, int nsteps_check = 3, double maxnpincr = 1e-3, double ptoll = 1e-10, 
				double toll = 1e-4, double dt0 = 1e-10, double dtcut = 0.25, double dtmax = 1, double dtmin = 1e-12, double maxdtincr = 2,	// Algor
				double Vshift = 1.79738, double Csb = 1.16183675549126e-11, double t_semic = 3.49436549222355e-8, double t_ins = 4.41e-7, 
				double L = 1e-5, bool ins = true,		// Device
				std::array<int,2> pins = {1, 0}, std::array<int,2> contacts = {0, 1}, double section = 0.00000081, double Vdrain = 5);
				// constructor

    /// METHODS

	std::vector<double>& get_data_phi_lumo();
	std::vector<double> get_data_n();
	
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

	protected:
		
		std::vector<double>	_data_phi_lumo;		/**< interpolation table */
		std::vector<double>	_data_n;			/**< interpolated values */
		Constants*			_cnst;
		Material*			_mat;
		Quad*				_quad;
		Algor*				_alg;
		Device*				_dev;
	
};

#endif /* PROBL_H */