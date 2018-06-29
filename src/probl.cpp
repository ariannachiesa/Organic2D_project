/*! \file probl.cpp
  \brief Class Problem
*/

#include "probl.h"

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

/// Method which sets T0 value and update Vth and sigman_kT values
void 
Probl::set_T0(double T0){
	_T0 = T0;
	_Vth = _Kb * _T0 / _q;
	_sigman_kT = _sigman / (_Kb * _T0);
};

void 
Probl::Constants(double T0){
	_Kb   = 1.380648813131e-23;
	_q    = 1.602176620898e-19;		
	_eps0 = 8.854187817620e-12;				
	_T0 = T0;
	_Vth = _Kb * _T0 / _q;
};

void 
Probl::Material(double PhiB, double sigman, double mu0n){
	_eps_semic_r = 2.90;			
	_eps_ins_r   = 2.82222;			
	_eps_semic = _eps0 * _eps_semic_r;
	_eps_ins   = _eps0 * _eps_ins_r;

	_PhiB = -PhiB; 
	_N0 = 1e27; 
	_Egap = 1.55; 

	_sigman = sigman * _Kb * 300;
	_sigman_kT = _sigman / (_Kb * _T0);
	_mu0n = mu0n;
};

void
Probl::Quad(int n){
	/// Quadrature nodes and weights
	double  gx[n],
			gw[n];
	webbur::hermite_compute (n, gx, gw);

	for(int i=0; i<n; i++){
		_gx.push_back( gx[i] );
		_gw.push_back( gw[i] );
	}
};

void
Probl::Algor(int pmaxit, double ptoll)
{	
	_pmaxit = pmaxit;
	_ptoll = ptoll;	
};

void
Probl::Device(	double Vshift, double Csb, double t_semic, double t_ins, double L, bool ins, std::array<int,2>& pins, 
				std::array<int,2>& contacts, double section, double Vdrain, double VG, double VB, int maxcycle)
{
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

	_VG = VG;
	_VB = VB;
	_Vshift = Vshift;
	_Csb = Csb;
	_t_semic = t_semic;
	_t_ins = t_ins;
	_L = L;
	
	std::vector<int>	row1, row2;
	
		const int	nNodes = 401;
		std::cout<<"Nnodes = "<<nNodes<<std::endl;
		
		_nTrees = nNodes-1;

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

		_msh.vtk_export ("Mesh");
		
		recursive = 0; partforcoarsen = 1;
		for (int cycle = 0; cycle < maxcycle; ++cycle)	// loop which refines the mesh uniformly
		{
			_msh.set_refine_marker (uniform_refinement);
			_msh.refine (recursive, partforcoarsen);
		}

		_msh.vtk_export ("Refined Mesh");
		
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

				if( (indexE == 2) && (indexT == 0)){
					row1.push_back( quadrant->gt(i) );
				}
				if( (indexE == 3) && (indexT == (nNodes-2)) ){
					row2.push_back( quadrant->gt(i) );
				}
			}
		}
	
	std::vector<int>::iterator it2, it3;
	
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
	
	if(! _ins){	// if there's not insulator
			
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
			
				if(x >= 0 && x <= L && y >= -t_semic && (y < 0.0 || y == 0.0)){
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
	_Efield = Vdrain / _L;
};


Probl::Probl(	int maxcycle,
				double T0,																						// Constants
				double PhiB, double sigman, double mu0n,														// Material
				int nq,																							// Quad
				int pmaxit, double ptoll,																		// Algor
				double Vshift, double Csb, double t_semic, double t_ins,
				double L, bool ins,																				// Device
				std::array<int,2> pins, std::array<int,2> contacts, double section, double Vdrain, double VG, double VB)
{
	
	Constants(T0);
	Material(PhiB, sigman, mu0n);
	Quad(nq);
	Algor(pmaxit, ptoll);
	Device(Vshift, Csb, t_semic, t_ins, L, ins, pins, contacts, section, Vdrain, VG, VB, maxcycle);

};


void
Probl::Laplace(){
	
	int	nnodes = _msh.num_global_nodes(),
		nelements = _msh.num_global_quadrants();
	
	std::vector<double>	phiout(nnodes,0.0),
						nout(nnodes,0.0),
						f(nnodes,0.0),
						epsilon(nelements,_eps_semic),
						psi(nnodes,0.0);
						
	sparse_matrix   A;
	A.resize(nnodes);
	
	///Assemble system matrices.	
	if(_ins){
        for(unsigned i=0; i<_insulator.size(); i++){
            if(_insulator[i]==1){
                epsilon[i] = _eps_ins;
            }
        }
    }
	
	bim2a_advection_diffusion (_msh, epsilon, psi, A, default_ord);
	
	int indexT = _nTrees-1;
	std::tuple<int, int, func_quad>	tupla1(0,2,[&](tmesh::quadrant_iterator quad, tmesh::idx_t i){return _PhiB + _VB;}),
									tupla2(indexT,3,[&](tmesh::quadrant_iterator quad, tmesh::idx_t i){return _Vshift + _VG;});
												
	dirichlet_bcs_quad	bcs;
	bcs.push_back(tupla1);
	bcs.push_back(tupla2);
	
	/// BCs Dirichlet type: phi(-t_semic) = PhiB ; phi(t_ins) = Vgate + Vshift;
	bim2a_dirichlet_bc (_msh,bcs,A,f, default_ord);
	
	/// Assembling rhs term:
	phiout = f;
	
	mumps mumps_solver;
      
	std::vector<double> vals(nnodes,0.0);
	std::vector<int> 	irow(nnodes,0),
						jcol(nnodes,0);
	  
	A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
	  
	mumps_solver.set_lhs_structure (A.rows(), irow, jcol);
		
	mumps_solver.analyze ();
	mumps_solver.set_lhs_data (vals);
      
	mumps_solver.factorize ();

	mumps_solver.set_rhs (phiout);
      
	mumps_solver.solve ();
	mumps_solver.cleanup ();
	
	/// Post-processing.	
	Vin = phiout;
	nin = nout;
	niter = 1;
};


void
Probl::LinearPoisson(){
	
	int	nnodes = _msh.num_global_nodes(),
		nelements = _msh.num_global_quadrants();
		
	double	areaS = _L * std::abs(_t_semic),
			Nd = 1e+19;
	
	std::vector<double>	phiout(nnodes,0.0),
						nout(nnodes,0.0),
						f(nnodes,0.0),
						rho(nnodes,_q/areaS*Nd),
						epsilon(nelements,_eps_semic),
						ecoeff(nelements,0.0),
						ncoeff(nnodes,1.0),
						psi(nnodes,0.0);
						
	sparse_matrix   A, M;
	A.resize(nnodes);
	M.resize(nnodes);
	
	///Assemble system matrices.	
	if(_ins){
        for(unsigned i=0; i<_insulator.size(); i++){
            if(_insulator[i]==1){
                epsilon[i] = _eps_ins;
            }
        }
    }
	
	bim2a_advection_diffusion (_msh, epsilon, psi, A, default_ord);
	
	for (unsigned i=0; i<_insulator.size(); i++){
		if(_insulator[i]==0){
			ecoeff[i] = 1;
		}
	}
	bim2a_reaction (_msh, ecoeff, ncoeff, M, default_ord);
	
	int indexT = _nTrees-1;	
	std::tuple<int, int, func_quad>	tupla1(0,2,[&](tmesh::quadrant_iterator quad, tmesh::idx_t i){return _PhiB + _VB;}),
									tupla2(indexT,3,[&](tmesh::quadrant_iterator quad, tmesh::idx_t i){return _Vshift + _VG;});
												
	dirichlet_bcs_quad	bcs;
	bcs.push_back(tupla1);
	bcs.push_back(tupla2);

	for(unsigned i=0; i<_scnodes.size(); i++){
		if(_scnodes[i] == 0){
			rho[i] = 0;
		}
	}
	
	f = M*rho;
	for(unsigned i=0; i<f.size(); i++){
		f[i] *= (-1);
	}
	
	/// BCs Dirichlet type: phi(-t_semic) = PhiB ; phi(t_ins) = Vgate + Vshift;
	bim2a_dirichlet_bc (_msh,bcs,A,f, default_ord);
	
	/// Assembling rhs term:
	phiout = f;
	
	mumps mumps_solver;
      
	std::vector<double> vals(nnodes,0);
	std::vector<int> 	irow(nnodes,0),
						jcol(nnodes,0);
	  
	A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
	  
	mumps_solver.set_lhs_structure (A.rows(), irow, jcol);
		
	mumps_solver.analyze ();
	mumps_solver.set_lhs_data (vals);
      
	mumps_solver.factorize ();

	mumps_solver.set_rhs (phiout);
      
	mumps_solver.solve ();
	mumps_solver.cleanup ();
	
	/// Post-processing.	
	for(unsigned i=0; i<_scnodes.size(); i++){
		if(_scnodes[i]==1){
			nout[i] = - rho[i]/_q;
		}
	}
	
	Vin = phiout;
	nin = nout;
	niter = 1;
};


void
Probl::NonLinearPoisson(std::vector<double>& phi0){

	int	nnodes = _msh.num_global_nodes(),
		nelements = _msh.num_global_quadrants(),
		iter;

	std::vector<double>	phiout(nnodes,0.0),
						phi(nnodes,0.0),
						dphi(nnodes,0.0),
						res(nnodes,0.0),
						res1(nnodes,0.0),
						res2(nnodes,0.0),
						rho(phiout.size(),0.0),
						drho(phiout.size(),0.0),
						nout(nnodes,0.0),
						epsilon(nelements,_eps_semic),
						ecoeff(nelements,0.0),
						ncoeff(nnodes,1.0),
						psi(nnodes,0.0);
						
	sparse_matrix   jac, A, M;
	
	jac.resize(nnodes);
	A.resize(nnodes);
	M.resize(nnodes);
	resnrm.resize(_pmaxit);
	
	///Assemble system matrices.	
	if(_ins){
        for(unsigned i=0; i<_insulator.size(); i++){
            if(_insulator[i]==1){
                epsilon[i] = _eps_ins;
            }
        }
    }
	bim2a_advection_diffusion (_msh, epsilon, psi, A, default_ord, default_ord);
	
	for (unsigned i=0; i<_insulator.size(); i++){
		if(_insulator[i]==0){
			ecoeff[i] = 1;
		}
	}
	bim2a_reaction (_msh, ecoeff, ncoeff, M, default_ord, default_ord);	
	
	int indexT = _nTrees-1;
	double	Vshift = _Vshift,
			Vg = _VG;
	std::tuple<int, int, func_quad>	tupla1(0,2,[](tmesh::quadrant_iterator quad, tmesh::idx_t i){return 0.0;}),
									tupla2(indexT,3,[&Vg,&Vshift,&phiout](tmesh::quadrant_iterator quad, tmesh::idx_t i)
													{return (Vg+Vshift-phiout[quad->gt(i)]);});
												
	dirichlet_bcs_quad	bcs;
	bcs.push_back(tupla1);
	bcs.push_back(tupla2);
	
	phi = phi0;
	phiout = phi0;
	
	for(iter=1; iter<=_pmaxit; iter++){

		phiout = phi;
		
		org_gaussian_charge_n(phiout, rho,drho);
		
		res1 = A*phiout;
		res2 = M*rho;
		for(unsigned i=0; i<res.size(); i++){
			res[i] = res2[i] - res1[i];
		}
		res1.clear();
		res2.clear();
		
		jac = A;
		for(int i=0; i<nnodes; i++){
			jac[i][i] -= M[i][i]*drho[i];
		}

		
		/// BCs Dirichlet type: phi(-t_semic) = PhiB ; phi(t_ins) = Vgate + Vshift;
		bim2a_dirichlet_bc (_msh,bcs,jac,res, default_ord);
	
		/// Assembling rhs term:
		dphi = res;
		
		mumps mumps_solver;
      
		std::vector<double> vals(nnodes,0);
		std::vector<int> 	irow(nnodes,0),
							jcol(nnodes,0);
	  
		jac.aij(vals, irow, jcol, mumps_solver.get_index_base ());
	  
		mumps_solver.set_lhs_structure (jac.rows(), irow, jcol);
		
		mumps_solver.analyze ();
		mumps_solver.set_lhs_data (vals);
      
		mumps_solver.factorize ();

		mumps_solver.set_rhs (dphi);
      
		mumps_solver.solve ();
		mumps_solver.cleanup ();
		
		for(unsigned i=0; i<phi.size(); i++){
			phi[i] += dphi[i];
		}
		
		double norm;
		Norm (_msh,dphi,norm,Inf);
		resnrm[iter-1] = norm;
		std::cout<<"norm = "<<norm<<std::endl;
		
		if(resnrm[iter-1]<_ptoll){
			std::cout<<"Poisson: resnrm < ptoll"<<std::endl;
            break;
		}
	}
	
	phiout = phi;
	
	/// Post-processing.
	
	org_gaussian_charge_n(phiout, rho,drho);
	drho.clear();

	for(unsigned i=0; i<_scnodes.size(); i++){
		if(_scnodes[i]==1){
			nout[i] = - rho[i]/_q;
		}
	}
	rho.clear();
	
	Vin = phiout;
	nin = nout;
	niter = iter;
	resnrm.resize(iter);

};


void
Probl::savePoisson(std::vector<double>& V, std::vector<double>& n, double niter, std::vector<double>& resnrm, const char* FileName)
{
  ColumnVector oct_V (V.size (), 0.0);
  ColumnVector oct_n (n.size (), 0.0);
  ColumnVector oct_res (resnrm.size (), 0.0);

  std::copy_n (V.begin (), V.size (), oct_V.fortran_vec ());
  std::copy_n (n.begin (), n.size (), oct_n.fortran_vec ());
  std::copy_n (resnrm.begin (), resnrm.size (), oct_res.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("niter", niter);
  the_map.assign ("V", oct_V);
  the_map.assign ("n", oct_n);
  the_map.assign ("res", oct_res);
  
  octave_io_mode m = gz_write_mode;
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);
  assert (octave_save ("Poisson", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};


void
Probl::saveCV(std::vector<double>& V, std::vector<double>& C, const char* FileName)
{
  ColumnVector oct_V (V.size (), 0.0);
  ColumnVector oct_C (C.size (), 0.0);

  std::copy_n (V.begin (), V.size (), oct_V.fortran_vec ());
  std::copy_n (C.begin (), C.size (), oct_C.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("V", oct_V);
  the_map.assign ("C", oct_C);
  
  octave_io_mode m = gz_write_mode;
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);
  assert (octave_save ("CVcurve", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};


void
Probl::Norm (tmesh& msh, const std::vector<double>& v, double& norm, norm_type type)
{
  if (type == Inf)
    {
      norm = 0.0;
      for (unsigned int i = 0; i < v.size (); ++i)
        {
          double temp = std::fabs (v[i]);
          if (norm < temp)
            norm = temp;
        }
    }
  else if (type == L2 || type == H1)
    {
      norm = 0.0;
      std::vector<double> ecoeff (msh.num_local_quadrants (), 1.0);
      std::vector<double> ncoeff (msh.num_local_nodes (), 1.0);
	  std::vector<double> psi (msh.num_local_nodes (),0);
      sparse_matrix M;
	  M.resize(msh.num_local_nodes ());
      bim2a_reaction (msh, ecoeff, ncoeff, M);
      if (type == H1)
        bim2a_advection_diffusion (msh, ecoeff, psi, M);
      std::vector<double> temp;
      temp = M * v;
      for (unsigned int i = 0; i < v.size (); ++i)
        norm += v[i] * temp[i];
      norm = sqrt (norm);
    }
};

void
Probl::org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& rhon, std::vector<double>& drhon_dV)
{
    std::vector<double> n = n_approx(V);
	
	rhon = n;
    for(unsigned i=0; i<n.size(); i++){
		rhon[i] *= -_q;
    }

    std::vector<double> dn_dV = dn_dV_approx(V);
	
	drhon_dV = dn_dV;
    for(unsigned i=0; i<dn_dV.size(); i++){
		drhon_dV[i] *= -_q;
    }
};

void
Probl::org_gaussian_charge_n(std::vector<double>& V, std::vector<double>& drhon_dV)
{
    std::vector<double> dn_dV = dn_dV_approx(V);
	
	drhon_dV = dn_dV;
    for(unsigned i=0; i<dn_dV.size(); i++){
		drhon_dV[i] *= -_q;
    }
};

std::vector<double>
Probl::n_approx(std::vector<double>& V)
{
    std::vector<double> coeff(V.size(),0),
						n(V.size(),0);
    double	kT = _Kb * _T0,
			denom;
	
    for(unsigned i=0; i<_gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){		
            coeff[j] = (sqrt(2) * _sigman * _gx[i] - _q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            n[j] += _N0 / sqrt(M_PI) * _gw[i] / denom;
        }
    }
	coeff.clear();		
	return n;
};

std::vector<double>
Probl::dn_dV_approx(std::vector<double>& V)
{
    std::vector<double> coeff(V.size(),0.0),
						dn_dV(V.size(),0.0);
						
    double	kT = _Kb * _T0,
			denom;

    for(unsigned i=0; i<_gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){
            coeff[j] = (sqrt(2) * _sigman * _gx[i] - _q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            dn_dV[j] += - _q * _N0 / _sigman * sqrt(2/M_PI) * _gw[i]*_gx[i] / denom;
		}
    }
	coeff.clear();
	return dn_dV;
};


double
Probl::CVcurve(std::vector<double>& phi){
    
    int nnodes = _msh.num_global_nodes(),
        nelements = _msh.num_global_quadrants();
		
	double dQ, sum, C;
    
    std::vector<double>	epsilon(nelements,_eps_semic),
                        dphi(nnodes,0.0),
                        rhon(nnodes,0.0),
                        drhon_dV(nnodes,0.0),
						psi(nnodes,0.0),
						ecoeff(nelements,0.0),
						ncoeff(nnodes,1.0),
						f(nnodes,0.0),
						S((_dnodes.size()*_dnodes[1].size()),0.0);
						
	sparse_matrix::col_iterator J;
    
    sparse_matrix   A1, M1, A2, M2, jac;
   
	A1.resize(nnodes);
	A2.resize(nnodes);
	M1.resize(nnodes);
	M2.resize(nnodes);
	jac.resize(nnodes);
    
    /// Assemble system matrices
    if(_ins){
        for(unsigned i=0; i<_insulator.size(); i++){
            if(_insulator[i]==1){
                epsilon[i] = _eps_ins;
            }
        }       
    }
    bim2a_advection_diffusion (_msh, epsilon, psi, A1);
	
	for (unsigned i=0; i<_insulator.size(); i++){
		if(_insulator[i]==0){
			ecoeff[i] = 1;
		}
	}
	bim2a_reaction (_msh, ecoeff, ncoeff, M1);

    org_gaussian_charge_n(phi, rhon, drhon_dV);
    rhon.clear();
	
	jac = A1;
	for(int i=0; i<nnodes; i++){
		jac[i][i] -= M1[i][i]*drhon_dV[i];
	}	

    int indexT = _nTrees-1;
    std::tuple<int, int, func_quad> tupla1(0,2,[](tmesh::quadrant_iterator quad, tmesh::idx_t i){return 0.0;}),
                                    tupla2(indexT,3,[](tmesh::quadrant_iterator quad, tmesh::idx_t i){return 1.0;});
                                                
    dirichlet_bcs_quad  bcs;
    bcs.push_back(tupla1);
    bcs.push_back(tupla2);
	
	/// BCs Dirichlet type: dphi(-t_semic) = 0 ; dphi(t_ins) = 1;
	bim2a_dirichlet_bc (_msh, bcs, jac, f, default_ord);
	
	/// Assembling rhs term:
	dphi = f;
	
	mumps mumps_solver;
      
	std::vector<double> vals(nnodes,0);
	std::vector<int> 	irow(nnodes,0),
						jcol(nnodes,0);
	  
	jac.aij(vals, irow, jcol, mumps_solver.get_index_base ());
	  
	mumps_solver.set_lhs_structure (jac.rows(), irow, jcol);
		
	mumps_solver.analyze ();
	mumps_solver.set_lhs_data (vals);
      
	mumps_solver.factorize ();

	mumps_solver.set_rhs (dphi);
      
	mumps_solver.solve ();
	mumps_solver.cleanup ();
	
	for(int i=0; i<nnodes; i++){
		dphi[i] *= (-1);
	}	
	
	/// Riassemblo le matrici
	bim2a_advection_diffusion (_msh, epsilon, psi, A2);
	bim2a_reaction (_msh, ecoeff, ncoeff, M2);
	
	jac = A2;
	for(int i=0; i<nnodes; i++){
		jac[i][i] -= M2[i][i]*drhon_dV[i];
	}
	
	for(unsigned i=0; i<_dnodes[1].size(); i++){
		sum = 0;
		for(J = jac[i].begin(); J!= jac[i].end(); ++J){
			sum += jac[i][jac.col_idx(J)]*dphi[jac.col_idx(J)];
		}
		S[i] = sum;
	}
	
	dQ = 0;
	for(unsigned i=0; i<S.size(); i++){
		dQ += S[i];
	}
	
	C = _L*dQ;
	
	return C;
};