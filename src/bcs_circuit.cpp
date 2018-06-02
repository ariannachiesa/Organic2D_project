/*! \file bcs_circuit.cpp
  \brief Enforce boundary conditions of the external control circuit
*/

#include "bcs_circuit.h"

BCS_CIRC::BCS_CIRC(	double freq, int voltage, double Csb, double Vshift, std::vector<double>& F, int Rg, int Rb, double amp)
{	
	_A.resize(4);
	_B.resize(4);

	_A[0][2] = 1;
	_B[0][3] = -1;
	_B[1][0] = 1;
	_B[2][1] = 1;
	_B[3][2] = 1;	
	
	_F = F;

	_C.resize(4);
	
	_r.resize(4);
	
	_r[0][0] = 1;
	_r[1][0] = Rg;
	_r[2][1] = Rb;
	
	_amp = amp;
	_freq = freq;
};

void BCS_CIRC::assign(double t, double Vshift, double Csb, std::vector<double>& F){
	std::vector<double>	s(4,0);
	std::vector<double>	C(4,0);
	double	VG;
	bool d1, d2, d3;
	
	d1 = ((t >= -90) && (t < -50));
	d2 = ((t >= -50) && (t <   0));
	d3 = (t >=   0);
  
	VG =(_voltage * (t + 90) / 40) * d1 + _voltage * d2 + (_voltage + _amp * sin (2 * M_PI * _freq * t)) * d3;
	//VG = _voltage * d1 + _voltage * d2 + (_voltage + _amp * sin (2 * M_PI * _freq * t)) * d3;
	VG = VG - Vshift;
	
	s[1] = -VG;
	s[3] = - VG * Csb;
	_F = F;
	C = _B*_F;	
	for(unsigned i=0; i<4; i++){
		C[i] += s[i];
	}
	_C = C;
	
	s.clear();
	C.clear();
};

void BCS_CIRC::get_A(sparse_matrix& mat){

	for(unsigned i=0; i<4; i++){
		for(unsigned j=0; j<4; j++){
			mat[i][j] = _A[i][j];
		}
	}
};

void BCS_CIRC::get_B(sparse_matrix& mat){
	
	for(unsigned i=0; i<4; i++){
		for(unsigned j=0; j<4; j++){
			mat[i][j] = _B[i][j];
		}
	}
};

std::vector<double> BCS_CIRC::get_C(){
	return _C;
};

std::vector<double> BCS_CIRC::get_F(){
	return _F;
};

void BCS_CIRC::get_r(sparse_matrix& v){

	for(unsigned i=0; i<_r.size(); i++){
		for(unsigned j=0; j<_r.size(); j++){	// potrebbero creare problemi le dimensioni; r così è quadr e non rettang
			v[i][j] = _r[i][j];
		}
	}
};