/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test_CVcurve.cpp
  \brief Test : CV curve
*/

#include "probl.h"

int main(int argc, char** argv){
	
	MPI_Init(&argc,&argv);
	
	clock_t tstart_probl, tstart_cv;
	tstart_probl = clock();
	
	/// n. refinement cycles : 0
	/// the mesh is constituted only by a row of quadrants along the y-axis
	Probl P(0);	

	tstart_probl = clock() - tstart_probl;
	std::cout<<"Construction class Probl run time: "<<tstart_probl<<" , ("<<((float)tstart_probl)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	int	nnodes = P._msh.num_global_nodes();
	std::vector<double>	Vguess(nnodes, P._PhiB);
	
	/// Export nodal field Vguess to a octbin.gz file for visualization.
	P._msh.octbin_export ("Vguess_visualization", Vguess);
	
	double	VGstart = -30,
			VGend = 30,
			dVG = 0.1;

	/// CV curve
	tstart_cv = clock();
	P.CVcurve (Vguess,VGstart,VGend,dVG,"CVcurve.gz");	
	tstart_cv = clock() - tstart_cv;
	std::cout<<"CV curve run time: "<<tstart_cv<<" , ("<<((float)tstart_cv)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};