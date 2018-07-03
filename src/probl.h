/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file probl.h
*  \brief Class Problem
* It stores informations about the device, material parameters, physical constans, nodes and wiegths for Gaussian quadrature rules,
* tolerances for algorithms and methods which solve physical problems. 
* Among the input parameters for its constructor: number of refinement cycles for the mesh.
* N.B : if n.cycles = 0 then the mesh will be constitute by a single row of elements on the y-axis
* otherwise it will be refined n.cycles-times and will have 2*n.cycles elements on both the x-axis and the y-axis
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
 * It stores informations about the device, material parameters, physical constans, nodes and wiegths for Gaussian quadrature rules,
 * tolerances for algorithms and methods which solve physical problems. 
 * Among the input parameters for its constructor: number of refinement cycles for the mesh.
 * N.B : if n.cycles = 0 then the mesh will be constitute by a single row of elements on the y-axis
 * otherwise it will be refined n.cycles-times and will have 2*n.cycles elements on both the x-axis and the y-axis
 */
class Probl
{
public:
	
  double _T0;			///< Absolute temperature
  double _Kb;			///< Boltzmann's constant 
  double _q;			///< Electron charge 
  double _eps0;			///< Vacuum permittivity 
  double _Vth;			///< Threshold voltage 

  
  double _eps_semic_r;		///< Relative permittivity of semiconductor 
  double _eps_ins_r;	        ///< Relative permittivity of insulator 
  double _eps_semic;		///< Permittivity of semiconductor 
  double _eps_ins;		///< Permittivity of insulator 
  double _PhiB;		        ///< Metal to semiconductor barrier. [eV] 
  double _N0;		        ///< Total number of available states (per unit volume) [m^{-3}] 
  double _Egap;			///< Energy gap. [eV] 
  double _sigman;		///< Disorder parameter [J] 
  double _sigman_kT;            ///< sigman / (Kb * T0)
  double _mu0n;		        ///< Low-field and low-density charge mobility for electrons [m^2 V^{-1} s^{-1}] 
  double _ni;                   ///< Intrinsic concentration of electrons
		

  std::vector<double> _gx;	///< Vector of nodes for quadrature rules 
  std::vector<double> _gw;	///< Vector of weigths for quadrature rules
	

  int _pmaxit;                  ///< Maximum allowed number of iterations
  double _ptoll;                ///< Tolerances on th residual norm of the increment
	

  bool _ins;		        ///< true = (default) if insulator is present , false = if there's only semiconductor	
  double _section;		///< Device section [m^2]
  double _Efield;	        ///< Device source-drain electric field [V / m]	
  double _Vshift;		///< Shift potential sensed at the gate terminal [V]
  double _VG;			///< Voltage applied at the gate terminal [V]
  double _VB;			///< Voltage applied at the bulk terminal [V]
  double _t_semic;		///< Thickness of semiconductor, along y-axis [m]
  double _t_ins;		///< Thickness of insulator, along y-axis [m]
  double _L;			///< Width of device, along x-axis [m]
  int _nTrees;			///< Number of Trees in the mesh
  std::vector<int> _scnodes;	                ///< Boolean vector; 1 = node in the semiconductor, 0 = node in the insulator
  std::vector<int> _insulator;	                ///< Boolean vector; 0 = element in the semiconductor, 1 = element in the insulator
  std::array<int,2> _pins;	                ///< [Bulk, gate] : region 1 = semiconductor , region 0 = insulator
  std::array<int,2> _contacts;	                ///< Number of the geometrical border containing the side edges where the contacts are: edge 2 of tree 0, edge 3 of last tree
  std::vector< std::vector<int> > _dnodes;	///< Vector with gate nodes + vector with bulk nodes
  tmesh	_msh;					///< Quadrangular mesh

  
  std::vector<double> Vin;		       ///< Output potential 
  std::vector<double> nin;		       ///< Output electron density 
  std::vector<double> resnrm;		       ///< Residual norm of the increment
  int niter;			               ///< Number of iterations required until convergence 

  
  /// Constructor.
  Probl(int maxcycle,
	double T0 = 300,						                 
	double PhiB = 0.54, double sigman = 2.6, double mu0n = 4.29110133911508e-6,             
	int nq = 101,								                 
	int pmaxit = 1000, double ptoll = 1e-10,				               
	double Vshift = 1.79738, double t_semic = 3.49436549222355e-8, double t_ins = 4.41e-7, double L = 1.4e-3, bool ins = true, std::array<int,2> pins = {1, 0},
	std::array<int,2> contacts = {2, 3}, double section = 0.00000081, double Vdrain = 5, double VG = 0.0, double VB = 0.0);

  /// Set physical constants.
  void Constants(double T0);

  /// Set material parameters.
  void Material(double PhiB, double sigman, double mu0n);

  /// Computes nodes and wigths for Gaussian quadrature rules using Sandia library.
  void Quad(int n);

  /// Set tolerances for Newton's method.
  void Algor(int pmaxit, double ptoll);

  /// Set geometrical/physical parameters of the device.
  void Device(double Vshift, double t_semic, double t_ins, double L, bool ins, std::array<int,2>& pins,
	      std::array<int,2>& contacts, double section, double Vdrain, double VG, double VB, int maxcycle);

  /// Solve Laplace problem.
  void Laplace();

  /// Solve Linear Poisson problem.
  void LinearPoisson();

  /// Solve Non-Linear Poisson problem.
  void NonLinearPoisson(std::vector<double>& phi0);

  /// Computes electron charge density and its derivative depending on the potential.
  void org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& rhon, std::vector<double>& drhon_dV);
  void org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& drhon_dV);

  /// Computes the approximation of n with Gaussian formulas.
  std::vector<double> n_approx(std::vector<double>& V);

  /// Computes the approximation of the derivative of n with Gaussian formulas.
  std::vector<double> dn_dV_approx(std::vector<double>& V);

  /// Computes capacitance of the MIS capacitor.
  double CVcurve(std::vector<double>& phi);
	
  /// Save output of the Poisson methods.
  void savePoisson(std::vector<double>& V, std::vector<double>& n, double niter, std::vector<double>& resnrm, const char* FileName);

  /// Save output of the CVcurve method.
  void saveCV(std::vector<double>& V, std::vector<double>& C, const char* FileName);

  /// Compute the (Inf,L2,H1) norm of a piecewise linear function.
  void Norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type);
	
  /// Set T0 method.
  /** It sets T0 value and update all the other parameters that depend on the temperature:
      Vth and sigman_kT
  */
  void set_T0(double T0);

};

#endif /* PROBL_H */
