/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test_NLPoisson.cpp
  \brief Test: Non Linear Poisson with applied V at gate
*/

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);
	
	clock_t tstart_probl, tstart_p;
	tstart_probl = clock();
	
	/// n. refinement cycles : 1
	/// the mesh is constituted only by two rows of quadrants along the y-axis
	Probl P(1);	
	P.set_T0(295);

	tstart_probl = clock() - tstart_probl;
	std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	int	nnodes = P._msh.num_global_nodes();
	std::vector<double>	Vguess(nnodes, P._PhiB);
	
	/// Export nodal field Vguess to a octbin.gz file for visualization.
	P._msh.octbin_export ("Vguess_visualization", Vguess);
	
	/// Solve Non Linear Poisson
	tstart_p = clock();

	P.set_VG(35);
	P.NonLinearPoisson(Vguess);
	P.savePoisson(P.Vin, P.nin, P.niter, P.resnrm, "NLPoisson_output.gz");
	
	tstart_p = clock() - tstart_p;
	std::cout<<"Non-Linear Poisson run time: "<<tstart_p<<" , ("<<((float)tstart_p)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};