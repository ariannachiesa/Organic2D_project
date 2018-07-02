/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test2_NLPoisson.cpp
  \brief Test: Non Linear Poisson with applied V at gate
*/

#include "probl.h"

int main(int argc, char** argv){
	
  MPI_Init(&argc,&argv);
	
  clock_t tstart_probl, tstart_p;
  tstart_probl = clock();
	
  /// n. refinement cycles : 1
  /// the mesh is constituted only by two rows of quadrants along the y-axis
  Probl P(1, 295, 0.54, 2.6, 4.29110133911508e-6, 101, 1000, 1e-10, 1.79738, 1.16183675549126e-11,
	  3.49436549222355e-8, 4.41e-7, 1.4e-3, true, {1, 0}, {2, 3}, 0.00000081, 5, 35, 0.0);

  tstart_probl = clock() - tstart_probl;
  std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
  int nnodes = P._msh.num_global_nodes();
  
  std::vector<double> Vguess(nnodes, P._PhiB);
	
  /// Export nodal field Vguess to a octbin.gz file for visualization.
  P._msh.octbin_export ("Vguess_visualization", Vguess, default_ord);
	
  /// Solve Non Linear Poisson problem with custom boundary conditions
  tstart_p = clock();

  P.NonLinearPoisson(Vguess);
  P.savePoisson(P.Vin, P.nin, P.niter, P.resnrm, "NLPoisson_output");
	
  tstart_p = clock() - tstart_p;
  std::cout<<"Non-Linear Poisson run time: "<<tstart_p<<" , ("<<((float)tstart_p)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
  std::cout<<"End of program"<<std::endl;
	
  MPI_Finalize();
	
  return 0;
};
