/*! \file csvread.h
  \brief Read and store input data from file csv
*/

#include "bim_sparse.h"

#include <fstream>
#include <istream>

/// It reads input data file
void
csvread(std::ifstream& file,sparse_matrix& CF);