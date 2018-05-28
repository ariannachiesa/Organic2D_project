/*! \file bcs_circuit.h
  \brief Enforce boundary conditions of the external control circuit
*/

#ifndef BCS_CIRC_H
#define BCS_CIRC_H

#include <bim_sparse.h>
#include <math.h>
#include <vector>

///	Compute matrices and vectors to enforce boundary conditions 
///	of the external control circuit attached to the device
class BCS_CIRC
{	
	sparse_matrix _A;
	sparse_matrix _B;
	sparse_matrix _r;
	std::vector<double> _C;	
	std::vector<double> _F;
	double _amp;
	double _freq;
	int _voltage;
	
	
    public:
	
    BCS_CIRC() = delete;
	BCS_CIRC(	double freq, const int voltage, double Csb, double Vshift, std::vector<double>& F,
				const int Rg = 0, const int Rb = 0, double amp = 0.1);
	
	//friend class Newton;
	
	void get_A(sparse_matrix& mat);
	void get_B(sparse_matrix& mat);
	std::vector<double> get_C();
	std::vector<double> get_F();
	//void get_r(std::vector<std::vector<double>>& v);
	void get_r(sparse_matrix& v);
	
	void assign(double t, double Vshift, double Csb, std::vector<double>& F);
};

#endif // BCS_CIRC_H