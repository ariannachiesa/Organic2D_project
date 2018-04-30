/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test3_mesh.cpp
  \brief Test 3: mesh creation
*/

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);
	
	clock_t tstart_probl;
	tstart_probl = clock();
	
	/// n. refinement cycles : 5
	/// very dense mesh
	Probl P(5);
	
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