/*! \file probl.cpp
  \brief Class Problem
*/

#include "probl.h"

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

Probl::Constants::Constants(double T0){
	_Kb   = 1.380648813131e-23;
	_q    = 1.602176620898e-19;		
	_eps0 = 8.854187817620e-12;				
	_T0 = T0;
	_Vth = _Kb * _T0 / _q;
};
	
Probl::Material::Material(Constants c, double PhiB, double sigman, double mu0n){
	_eps_semic_r = 2.90;			
	_eps_ins_r   = 2.82222;			
	_eps_semic = c._eps0 * _eps_semic_r;
	_eps_ins   = c._eps0 * _eps_ins_r;

	_PhiB = -PhiB; 
	_N0 = 1e27; 
	_Egap = 1.55; 

	_sigman = sigman * c._Kb * 300;	
	_sigman_kT = _sigman / (c._Kb * c._T0);
	_mu0n = mu0n;
};

Probl::Quad::Quad(int n){
	/// Quadrature nodes and weights
	double  gx[n],
			gw[n];
	webbur::hermite_compute (n, gx, gw);

	for(int i=0; i<n; i++){
		_gx.push_back( gx[i] );
		_gw.push_back( gw[i] );
	}
};

// std::vector<double> 
// Probl::Quad::get_gx(){
	// return _gx;
// };

// std::vector<double> 
// Probl::Quad::get_gw(){
	// return _gw;
// };

// double* 
// Probl::Quad::get_pgx(){
	// return &_gx[0];
// };

// std::vector<double>* 
// Probl::Quad::get_pgw(){
	// return &_gw;
// };

Probl::Algor::Algor(	int pmaxit, int maxit, int maxit_mnewton, int nsteps_check, double maxnpincr, double ptoll, 
						double toll, double dt0, double dtcut, double dtmax, double dtmin, double maxdtincr){

	_pmaxit = pmaxit;
	_maxit = maxit;
	_maxit_mnewton = maxit_mnewton;
	_nsteps_check = nsteps_check;
	_maxnpincr = maxnpincr;
	_ptoll = ptoll;
	_toll = toll;
	_dt0 = dt0;
	_dtcut = dtcut;
	_dtmax = dtmax;
	_dtmin = dtmin;
	_maxdtincr = maxdtincr;

	// Scaling and clamping coefficients
	for (int i=0; i<4; i++){
		_colscaling.push_back(1);
		_rowscaling.push_back(1);
		_clamping.push_back(1);
	}
	_clampOnOff = true;
	_savedata = true;
};
	
Probl::Device::Device(	double Vshift, double Csb, double t_semic, double t_ins, double L, bool ins, std::array<int,2>& pins, 
						std::array<int,2>& contacts, double section, double Vdrain, int maxcycle){
	
	int	recursive, partforcoarsen;

	if( ins ){
		if(t_ins == 0){
			_t_ins = 0;
			_ins = false;
		}
		else{
			_t_ins = t_ins;
			_ins = ins;
		}
	}
	else{
		_t_ins = 0;
	}

	_Vshift = Vshift;
	_Csb = Csb;
	_t_semic = t_semic;
	_t_ins = t_ins;
	_L = L;
	
	std::vector<int>	row1, row2;
	
	if(maxcycle == 0){ 	// Mesh solo lungo l'asse y

		const int	nNodes = 801;
		std::cout<<"Nnodes = "<<nNodes<<std::endl;

		// Define mesh.
		if(ins){	// if insulator is present: mesh with both the insulator and the semiconductor
   
			int	nNodesSemic = std::floor (0.6 * nNodes);
			int	nNodesIns   = nNodes - nNodesSemic;

			std::cout<<"NnodesSemic = "<<nNodesSemic<<std::endl;
			std::cout<<"NnodesIns = "<<nNodesIns<<std::endl;

			constexpr p4est_topidx_t simple_conn_num_vertices = nNodes*2;
			constexpr p4est_topidx_t simple_conn_num_trees = nNodes-1;
	
			double simple_conn_p[simple_conn_num_vertices*2] = {};
			p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {};

			double	step = t_semic / (nNodesSemic-1);
			int		j = 0;

			for(int i=0; i<(nNodesSemic*4); i=i+4){
				if( i==0 ){
					simple_conn_p[0]	=	0;
					simple_conn_p[1]	=	-t_semic;
					simple_conn_p[2]	=	L;
					simple_conn_p[3]	=	-t_semic;
					j++;
				}
				else{
					simple_conn_p[i]	=	L;
					simple_conn_p[i+1]	=	-t_semic + j*step;
					simple_conn_p[i+2]	=	0;
					simple_conn_p[i+3]	=	-t_semic + j*step;
					j++;
				}
			}
	
			j=1;
			step = t_ins / (nNodesIns);
			for(int i=(nNodesSemic*4); i<(nNodesSemic*4+nNodesIns*4); i=i+4){
				simple_conn_p[i]	=	L;
				simple_conn_p[i+1]	=	j*step;
				simple_conn_p[i+2]	=	0;
				simple_conn_p[i+3]	=	j*step;
				j++;
			}

			for(int i=0; i<(simple_conn_num_trees*5); i=i+5){
				if( i==0 ){
					simple_conn_t[0]	=	1;
					simple_conn_t[1]	=	2;
					simple_conn_t[2]	=	3;
					simple_conn_t[3]	=	4;
					simple_conn_t[4]	=	1;
				}
				else{
					simple_conn_t[i]	=	simple_conn_t[i-2];
					simple_conn_t[i+1]	=	simple_conn_t[i-3];
					simple_conn_t[i+2]	=	simple_conn_t[i] + 1;
					simple_conn_t[i+3]	=	simple_conn_t[i+2] + 1;
					simple_conn_t[i+4]	=	1;
				}
			}
			_msh.read_connectivity (simple_conn_p, simple_conn_num_vertices, simple_conn_t, simple_conn_num_trees);
		}
		else{		// mesh with semiconductor only
			constexpr p4est_topidx_t simple_conn_num_vertices = 4;
			constexpr p4est_topidx_t simple_conn_num_trees = 1;
			const double simple_conn_p[simple_conn_num_vertices*2] = {0, 0, L, 0, L, t_semic, 0, t_semic};
			const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1};	
			_msh.read_connectivity (simple_conn_p, simple_conn_num_vertices, simple_conn_t, simple_conn_num_trees);
		}

		_msh.vtk_export ("mesh semic-ins");

		using idx_t = p4est_gloidx_t;
		idx_t			indexE;
		p4est_locidx_t	indexT;
		for (auto quadrant = _msh.begin_quadrant_sweep ();
			quadrant != _msh.end_quadrant_sweep ();
			++quadrant)
		{
			for(int i=0; i<4; i++){
				indexE = quadrant->e(i);			// index of the boundary which the i-th vertex of the quadrant lies on (in current tree)
				indexT = quadrant->get_tree_idx ();	// index of the current tree

				if( (indexE == 0) && (indexT == 0)){
					_alldnodes.push_back( quadrant->gt(i) );
					row1.push_back( quadrant->gt(i) );
				}
				if( (indexE == 1) && (indexT == (nNodes-2)) ){	
					_alldnodes.push_back( quadrant->gt(i) );
					row2.push_back( quadrant->gt(i) );
				}
			}
		}
	}
	else{	// Mesh rifinita uniforme lungo entrambi gli assi

		// Define mesh.
		if(_ins){	// if insulator is present: mesh with both the insulator and the semiconductor
			constexpr p4est_topidx_t simple_conn_num_vertices = 6;
			constexpr p4est_topidx_t simple_conn_num_trees = 2;
			const double simple_conn_p[simple_conn_num_vertices*2] = {0, -t_semic, L, -t_semic, L, 0, 0, 0, L, t_ins, 0, t_ins};
			const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1, 4, 3, 5, 6, 1};
			_msh.read_connectivity (simple_conn_p, simple_conn_num_vertices, simple_conn_t, simple_conn_num_trees);
		}
		else{		// mesh with semiconductor only
			constexpr p4est_topidx_t simple_conn_num_vertices = 4;
			constexpr p4est_topidx_t simple_conn_num_trees = 1;
			const double simple_conn_p[simple_conn_num_vertices*2] = {0, 0, L, 0, L, t_semic, 0, t_semic};
			const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1};	
			_msh.read_connectivity (simple_conn_p, simple_conn_num_vertices, simple_conn_t, simple_conn_num_trees);
		}

		recursive = 0; partforcoarsen = 1;
		for (int cycle = 0; cycle < maxcycle; ++cycle)	// loop which refines the mesh uniformly
		{
			_msh.set_refine_marker (uniform_refinement);
			_msh.refine (recursive, partforcoarsen);
		}

		_msh.vtk_export ("mesh semic-ins refined");

		using idx_t = p4est_gloidx_t;
		idx_t			indexE;
		p4est_locidx_t	indexT;
	
		for (auto quadrant = _msh.begin_quadrant_sweep ();
			quadrant != _msh.end_quadrant_sweep ();
			++quadrant)
		{
			for(int i=0; i<4; i++){
				indexE = quadrant->e(i);			// index of the boundary which the i-th vertex of the quadrant lies on (in current tree)
				indexT = quadrant->get_tree_idx ();	// index of the current tree
			
				if( (indexE == 0) && (indexT == 0)){
					_alldnodes.push_back( quadrant->gt(i) );
					row1.push_back( quadrant->gt(i) );
				}
				if( (indexE == 1) && (indexT == 1) ){	
					_alldnodes.push_back( quadrant->gt(i) );
					row2.push_back( quadrant->gt(i) );
				}
			}
		}
	}
	
	std::vector<int>::iterator it1, it2, it3;
		
	std::sort(_alldnodes.begin(),_alldnodes.end());
	it1 = std::unique (_alldnodes.begin(), _alldnodes.end());
	_alldnodes.resize( std::distance(_alldnodes.begin(),it1) );
	
	std::sort (row1.begin(), row1.end());				// row1 = vector of the nodes on boundary 0 of the semiconductor
	
	it2 = std::unique (row1.begin(), row1.end());		// removal of the common vertices between quadrants which have been
														//	  inserted more than once in the vector
	row1.resize( std::distance(row1.begin(),it2) );
		
	std::sort (row2.begin(), row2.end());				// row2 = vector of the nodes in the insulator
	it3 = std::unique (row2.begin(), row2.end());
	row2.resize( std::distance(row2.begin(),it3) );
	
	_dnodes.push_back(row1);
	_dnodes.push_back(row2);
		
	row1.clear();
	row2.clear();
		
	/// Contacts
	_pins = {pins[0], pins[1]};
	_contacts = {contacts[0], contacts[1]};
		
	double 	x = 0, y = 0, smc = 0;
	
	if(! _ins){
			
		for(auto quadrant = _msh.begin_quadrant_sweep ();
		quadrant != _msh.end_quadrant_sweep ();
		++quadrant){
			_insulator.push_back(0);
		}
		
		std::vector<int>	scnodes(_msh.num_global_nodes(),1);
		_scnodes = scnodes;
		scnodes.clear();
	}
	else{
	// scorro sui quadranti
	for(auto quadrant = _msh.begin_quadrant_sweep ();
		quadrant != _msh.end_quadrant_sweep ();
		++quadrant)
	{
		smc = 0;
		// per ogni quadrante scorro sui vertici
		for (int ii = 0; ii < 4; ++ii)
		{
			x = quadrant->p(0, ii);
			y = quadrant->p(1, ii);
			
			if(x >= 0 && x <= L && y >= -t_semic && y <= 0){
				smc++;
				row1.push_back( quadrant->gt(ii) );
			}
			else{
				row2.push_back( quadrant->gt(ii) );
			}
		}
		if(smc==4){
			_insulator.push_back(0);	// if all the vertices of the quadrant lie in the semiconductor
										//  --> the quadrant is in the semiconductor : 0
										//  --> otherwise it is in the insulator : 1
		}
		else{
			_insulator.push_back(1);
		}
	}
	
	std::sort (row1.begin(), row1.end());				// row1 = vector of the semiconductor nodes 
		
	it2 = std::unique (row1.begin(), row1.end());		// removal of the common vertices between quadrants which have been
														//	inserted more than once in the vector
	row1.resize( std::distance(row1.begin(),it2) );
	
	std::sort (row2.begin(), row2.end());				// row2 = vector of the insulator nodes 
	
	it3 = std::unique (row2.begin(), row2.end());
	row2.resize( std::distance(row2.begin(),it3) );
		
	// creo il vettore complessivo dei nodi: 1 = semiconductor , 0 = insulator
	std::vector<int>	scnodes(row1.size()+row2.size(),0);
	for(unsigned i=0; i<row1.size(); i++){
		scnodes[i] = 1;
	}
	_scnodes = scnodes;
	scnodes.clear();
	row1.clear();
	row2.clear();
	}
		
	_section = section;
	_Efield = Vdrain / L;
};

Probl::Probl(	double T0,	// Constants
				double PhiB, double sigman, double mu0n,		// Material
				int nq,		// Quad
				int pmaxit, int maxit, int maxit_mnewton, int nsteps_check, double maxnpincr, double ptoll,
				double toll, double dt0, double dtcut, double dtmax, double dtmin, double maxdtincr,	// Algor
				double Vshift, double Csb, double t_semic, double t_ins,
				double L, bool ins,		// Device
				std::array<int,2> pins, std::array<int,2> contacts, double section, double Vdrain, int maxcycle)
{	
	///	Calculation of the interpolation table
	int j = 0;
	double	a, b, N, step;
	a = -10;
	b = 10;
	N = 1e6+1;
	step = (b-a)/(N-1);
	_data_phi_lumo.resize(N);
	for(auto i=a; i<=b; i+=step){
		_data_phi_lumo[j] = ( i );
		j++;
	}

	std::cout<<"mu0n = "<<mu0n<<std::endl;
	
	Constants	Cnst(T0);
	Material	Mat(Cnst, PhiB, sigman, mu0n);
	Quad		Q(nq);
	
	///	Calculation of interpolated n
    std::vector<double> coeff(_data_phi_lumo.size(),0),
						n(_data_phi_lumo.size(),0),
						gx = Q._gx,
						gw = Q._gw;
						
    double	q = Cnst._q,
			kT = Cnst._Kb * Cnst._T0,
			N0 = Mat._N0,
			denom;
	
	_data_n.resize(_data_phi_lumo.size());
    for(unsigned i=0; i<gx.size(); i++){
        for(unsigned j=0; j<_data_phi_lumo.size(); j++){		
            coeff[j] = (sqrt(2) * sigman * gx[i] - q * _data_phi_lumo[j]) / kT ;
            denom = 1+exp(coeff[j]);
			n[j] += N0 / sqrt(M_PI) * gw[i] / denom;
			_data_n[j] = -q*n[j];
        }
    }

	gx.clear();
	gw.clear();
	coeff.clear();

};

std::vector<double>& Probl::get_data_phi_lumo(){
	return _data_phi_lumo;
};

std::vector<double> Probl::get_data_n(){
	return _data_n;
};