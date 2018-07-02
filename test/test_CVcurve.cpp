/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test_CVcurve.cpp
  \brief Test: CV curve
*/

#include "probl.h"

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  clock_t tstart_probl, tstart_p;
  tstart_probl = clock();

  /// n. refinement cycles : 1
  /// the mesh is constituted only by two rows of quadrants along the y-axis
  Probl P(1, 295);

  tstart_probl = clock() - tstart_probl;
  std::cout << "Construction class Probl run time: " << tstart_probl << " , (" << ((float)tstart_probl) / CLOCKS_PER_SEC << " seconds)." << std::endl;

  double VG_min = -35.0,
         VG_max = 35.0,
         dV = 1.0,
         c;

  int nnodes = P._msh.num_global_nodes(),
      length = std::abs(VG_max - VG_min) / dV;

  std::vector<double> Vguess(nnodes, P._PhiB),
                      C(length, 0.0),
                      V(length, 0.0),
                      Q;

  tstart_p = clock();
  P._VG = VG_min;
  for (int i = 0; i < length; i++) {

    std::cout << "iter = " << i << std::endl;

    P._VG += dV;
    V[i] = P._VG;

    /// Solve Non Linear Poisson problem
    P.NonLinearPoisson(Vguess);
    P.savePoisson(P.Vin, P.nin, P.niter, P.resnrm, "NLPoisson_output");

    /// Compute capacitance
    c = P.CVcurve(P.Vin);
    C[i] = c;

  }
  P.saveCV(V, C, "CVcurve_output");


  tstart_p = clock() - tstart_p;
  std::cout << "CV curve run time: " << tstart_p << " , (" << ((float)tstart_p) / CLOCKS_PER_SEC << " seconds)." << std::endl;

  std::cout << "End of program" << std::endl;

  MPI_Finalize();

  return 0;
};
