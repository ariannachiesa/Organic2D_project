/*
  Arianna Chiesa
  Project for the course:
  "Advanced Programming for Scientific Computing"
*/

/*! \file test_DD.cpp
  \brief Test: Non Linear Poisson + Drift-Diffusion
*/

#include "bcs_circuit.h"
#include "csvread.h"
#include "interp1.h"
#include "newton.h"
#include "probl.h"

#include <iostream>

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
	std::cout<<"NLPoisson run time: "<<tstart_p<<" , ("<<((float)tstart_p)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	///Other physical parameters
	double	sum, average;
	sum = std::accumulate( P._scnodes.begin(), P._scnodes.end(), 0.0);
	average = std::accumulate( P.nin.begin(), (P.nin.begin()+sum), 0.0)/sum;
	P.set_ni(average);
	
	///  2D simulation
	double  eps_ins0 = P._eps_ins,
			Csb0 = P._Csb;
	
	sparse_matrix CF;
	std::ifstream file;
	file.open("input_data.csv");
	if(!file){
		std::cout<<"error: the file does not exist"<<std::endl;
		exit (EXIT_FAILURE);
	}
	else{
		/// read file 'input_data.csv'
		csvread(file,CF);
	}
	
	/// Assembly of the frequency vector: frequency = [CF(2:end, 1); 2 * logspace(5, 11, 91)(2:end)']; # [Hz]
	double	a, b, step, freq,
			tmin, tmax,
			N = 91;
	int		j;
	std::complex<double> complx(0.0,1.0);
	std::vector<double>	frequency((N-1) + CF.size(),0.0),
						f(CF.size(),0.0),
						eps_ins_r(CF.size(),0.0),
						tspan,
						Fin(4,0.0);

	for(unsigned i=0; i<CF.size(); i++){
		frequency[i] = CF[i][0];
	}
	
	a = pow(10,5);
	b = pow(10,11);
	step = pow( (b/a) , (1/(N-1)) );

	j = CF.size();
    for(auto i=a*step; i<=b; i=i*step){
        frequency[j] = (2*i);
		j++;
    }
	
	for(int i=0; i<61; i++){
		f[i] = frequency[i];
	}
	
	for(unsigned i=0; i<CF.size(); i++){
			eps_ins_r[i] = CF[i][4];
	}
	
	clock_t tstart_loop,
            tstart_simul;

    tstart_loop = clock();	
	for (unsigned i=0; i<frequency.size(); i++){
	
		freq = frequency[i];
		/// Compute frequency-dependent eps_ins and Csb.
		if (freq <= frequency[60]){
			P.set_eps_ins_r( interp1( f, eps_ins_r, freq, true ) );
		}
        else{
			complx *= freq;
			complx *= 1e-2;
			complx += 1;
			P.set_eps_ins_r( 	std::real(2.56+(eps_ins_r[0]-2.56) /  std::pow(complx,0.24) )); //# [~]
        }
		
		P.set_eps_ins( P._eps0 * P._eps_ins_r ); // [F / m]
        P.set_Csb( Csb0 * P._eps_ins / eps_ins0 ); // [F]
	
		/// Time span for simulation
        tmin  = -100;
        tmax  = 5 / freq;
		a = 0;
		b = tmax;
		N = 1000;
		step = (b-a)/(N-1);
		
		tspan.resize(N+3);
		tspan[0] = tmin;	tspan[1] = -90;		tspan[2] = -50;
		j = 3;
		for(auto h=a; h<=b; h+=step){
			tspan[j] = h ;
			j++;
		}		

		/// Circuit boundary conditions and initial condition.
        /// Full circuit.
		
		Fin[0] = - P._Vshift;
        Fin[2] = - P._Vshift * P._Csb;
		
		tstart_simul = clock();
		
		/// Enforcing boundary conditions of the attached control circuit
		BCS_CIRC	bcs(freq, P._VG, P._Csb, P._Vshift, Fin);
		
		/// Newton's algorithm
		// Newton		newt(	P, Vin, nin, tspan, Fin, Iin, bcs, freq );
		
		tstart_simul = clock() - tstart_simul;
        std::cout<<"Simulation run time: "<<tstart_simul<<" , ("<<((float)tstart_simul)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	}

	tstart_loop = clock() - tstart_loop;
    std::cout<<"Loop run time: "<<tstart_loop<<" , ("<<((float)tstart_loop)/CLOCKS_PER_SEC<<" seconds)."<<std::endl;
	
	std::cout<<"End of program"<<std::endl;
	
	MPI_Finalize();
	
	return 0;
};