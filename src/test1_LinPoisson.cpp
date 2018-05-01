/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test1_LinPoisson.cpp
  \brief Test 1: Linear Poisson
*/

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);
	
	clock_t tstart_probl, tstart_nlp;
	tstart_probl = clock();
	
	/// n. refinement cycles : 1
	/// the mesh is constituted only by two rows of quadrants along the y-axis
	Probl P(1);	

	tstart_probl = clock() - tstart_probl;
	std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	int	nnodes = P._msh.num_global_nodes();
	std::vector<double>	Vguess(nnodes, P._PhiB);
	
	/// Export nodal field Vguess to a octbin.gz file for visualization.
	P._msh.octbin_export ("Vguess_visualization", Vguess);
	
	/// Solve Linear Poisson
	tstart_nlp = clock();
	P.NLPoisson(Vguess);
	P.saveNLP(P.Vin, P.nin, P.niter, P.res, "LinPoisson_output.gz");
	tstart_nlp = clock() - tstart_nlp;
	std::cout<<"NLPoisson run time: "<<tstart_nlp<<" , ("<<((float)tstart_nlp)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};