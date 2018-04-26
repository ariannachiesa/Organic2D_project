/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file mis_main2d_CF.cpp
  \brief Simulation of 2D Organic Thin Film Transistors
*/

#include<algorithm>
#include<ctime>
#include<complex>
#include<iostream>
#include<fstream>
#include <math.h>
#include <mpi.h>
#include<sstream>
#include<stdlib.h>
#include<string>
#include<vector>

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);

	/// Input Parameters to be specified by the user

	/// for Constants class:
    double    T0;// = 295;	// absolute temperature
	
	/// for Device class:
	/// Device dimensions
	bool	ins = true;			// true = (default) if insulator is present , false = if there's only semiconductor 
	double  t_ins,		// [m]		thickness of insulator, along y-axis 
			L, 		// [m]		width of the device, along x-axis 
			t_semic,		// [m]		thickness of semiconductor, along y-axis 
			section,		// [m^2]		device section 
			Vdrain = 0,			// [V]		Drain potential 
			Vshift = 0,			// [V]		shift potential sensed at the gate terminal 
			Csb = 0;			// [F]		value of the capacitor of the external control circuit 
	
	std::array<int,2>	pins,		// [Bulk, gate] , default: region 1 = semiconductor , region 0 = insulator 
						contacts;	// contacts between the device and the external control circuit 
						
	int	cycle = 0;				// number of cycles of mesh refinement; default 2 cycles = 45 nodes, 32 elements 
	
	/// for Quad class:
	int	n;			// Num. nodes/weigths for Quad Rules */
	
	/// for Algor class:
	/// Parameters for the algorithm
	int pmaxit,
		maxit,
		maxit_mnewton,
		nsteps_check;
		
	double	maxnpincr,
			ptoll,
			toll,
			dt0,
			dtcut,
			dtmax,
			dtmin,
			maxdtincr;
	
	/// for Material class:
	double  PhiB = 0.54,		// [eV]	Metal to semiconductor barrier. 
            sigman = 2.6,		// [J]	the disorder parameter 
			mu0n;				// [m^2 V^{-1} s^{-1}]	low-field and low-density charge mobility for electrons 
	
	clock_t tstart_probl;
	tstart_probl = clock();
	
	/// class Probl
	Probl 		P(	T0, PhiB, sigman, mu0n, n, pmaxit, maxit, maxit_mnewton, nsteps_check, maxnpincr, ptoll, toll, dt0, dtcut,
					dtmax, dtmin, maxdtincr, Vshift, Csb, t_semic, t_ins, L, ins, pins, contacts, section, Vdrain, cycle);

	tstart_probl = clock() - tstart_probl;
	std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};