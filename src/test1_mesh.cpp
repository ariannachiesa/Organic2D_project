/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file mis_main2d_CF.cpp
  \brief Simulation of 2D Organic Thin Film Transistors
*/

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);
	
	clock_t tstart_probl;
	tstart_probl = clock();
	
	/// class Probl :
	/// input parameter to be specified by the user: number of refinement cycles for the mesh.
	/// N.B : if n.cycles = 0 then the mesh will be constitute by a single row of elements on the y-axis
	/// otherwise it will be refined n.cycles times and will have 2*n.cycles elements on both the x-axis and the y-axis
	Probl P(0);
	
	P.set_T0(295);
	
	P.set_Vdrain(0);
	P.set_Csb(0);	
	P.set_Vshift(0);	

	tstart_probl = clock() - tstart_probl;
	std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};