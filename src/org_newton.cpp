/*! \file org_newton.cpp
  \brief Implementation of class Newton's methods
*/

#include "newton.h"

/**
 *	Compute residual vector, which is the output of the method
 */
std::vector<double>							
Newton::org_secs2d_newton_residual(	Probl& P, std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I,		
									std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0,	
									double deltat, BCS_CIRC& bcs, std::vector<int>& indexingV, std::vector<int>& indexingn, std::vector<int>& indexingF, 
									std::vector<int>& indexingI)
{
	int	nnodes = P.get_msh_nodes(),
		nelements = P.get_msh_elem(),
        ndofs,
		j = 0;
	double	eps_semic = P._eps_semic,
			eps_ins = P._eps_ins,
			q = P._q,
			Vth = P._Vth,
			Vshift = P._Vshift,
			L = P._L,
			tins = P._t_ins,
			tsemic = P._t_semic,
			section = P._section,
			PhiB = P._PhiB,
			s1, s2, s = 0;
	unsigned int	numcontacts;
	bool	ins = P._ins;
	std::array<int,2>	pins = P._pins;
	std::vector< std::vector<int> >	dnodes = P._dnodes;
	std::vector<int> 	insulator = P._insulator,
						scnodes = P._scnodes;
	std::vector<double>	epsilon(insulator.size(),eps_semic),
						rowscaling = P._rowscaling,
						res;
	
	numcontacts = pins.size();
	ndofs = 2 * nnodes + F.size() + numcontacts;		
	
	/// COMPUTING COEFFICIENTS
	if(ins){
		for(unsigned i=0; i<insulator.size(); i++){
			if(insulator[i] == 1){
				epsilon[i] = eps_ins;
			}
		}
	}

	///	COMPUTING FIRST ROW
	sparse_matrix 	M,
					A11,
					A12;
					
	M.resize(nnodes);
	A11.resize(nnodes);
	A12.resize(nnodes);
	
	std::vector<double> eones(insulator.size(),1.0),
						not_ins(insulator.size(),0.0),
						ones(nnodes,1.0),
						psi(nnodes,0.0);
	
	for(unsigned i=0; i<insulator.size(); i++){
		if(insulator[i] == 0){
			not_ins[i] = 1;
		}
	}
	
	bim2a_reaction (P._msh, not_ins, ones, A12);
	for(int i=0; i<nnodes; i++){
		A12[i][i] *= q ;			
	}

	bim2a_reaction (P._msh, eones, ones, M);		

	bim2a_advection_diffusion (P._msh, epsilon, psi, A11);		
	
	///	Physical models.
	std::vector<double>	mobn(nelements,0.0),
						alphan(n.size(),0.0),
						out(nelements,0.0);

	org_physical_models2d(n, P, mobn, alphan, out);
	out.clear();				/// fare l'overload del method!!!

	
	///	ASSEMBLING FIRST ROW
	
	res.resize(ndofs);
	std::vector<double>	sum1(V.size(),0.0),
						sum2(n.size(),0.0),
						resV(V.size(),0.0);

	sum1 = A11*V;	
	sum2 = A12*n;
	
	for(unsigned i=0; i<resV.size(); i++){
		if(sum1.size() != resV.size() || sum2.size() != resV.size()){
			std::cout<<"error: resV; wrong dimensions!"<<std::endl;
			exit(EXIT_FAILURE);
		}
		else{
			resV[i] = sum1[i] + sum2[i];
		}
	}

	///	ENFORCING BCs ON FIRST ROW
	
	std::vector<double>	BCbulk(V.size(),0.0),
						BCgate(V.size(),0.0);
						
	for(unsigned i=0; i<V.size(); i++){
		BCbulk[i] = (V[i] - F[pins[0]]) - PhiB;
	}
	BCbulk = M*BCbulk;
	
	for(unsigned i=0; i<V.size(); i++){
		BCgate[i] = (V[i] - F[pins[1]]) - (PhiB + Vshift);
	}
	BCgate = M*BCgate;

	int indexT = _nTrees-1;	
	std::tuple<int, int, func_quad>	tupla1(0,2,[&resV,&BCbulk](tmesh::quadrant_iterator quad, tmesh::idx_t i)
																{return (resV[quad->gt(i)]+BCbulk[quad->gt(i)]);}),
									tupla2(indexT,3,[&resV,&BCgate](tmesh::quadrant_iterator quad, tmesh::idx_t i)
																{return (resV[quad->gt(i)]+BCgate[quad->gt(i)];});
	dirichlet_bcs_quad	bcsV;
	bcsV.push_back(tupla1);
	bcsV.push_back(tupla2);
	
	bim2a_dirichlet_bc (P._msh,bcsV,M,resV);	/// che matrice metto qui ??
	
	for(unsigned i=0; i<indexingV.size(); i++){
		res[indexingV[i]] = resV[i];
	}
	resV.clear();
	
	///	COMPUTING SECOND ROW
	sparse_matrix	A22,
					Aa;
	std::vector<double>	resn(nnodes,0.0),
						rhs(nnodes, 0.0),
						alpha(insulator.size(),0.0),
						a_el(insulator.size(),0.0),
						eta(scnodes.size(),0.0),
						beta(nnodes,0.0),
						gamma(scnodes.size(),1.0),
						mob_nodes(scnodes.size(),0.0);
  
	A22.resize(nnodes);
	Aa.resize(nnodes);
		
	for(unsigned i=0; i<insulator.size(); i++){
		if(insulator[i] == 0){
			alpha[i] = mobn[i]*Vth;
		}
	}

	int index;
	for (	auto quadrant = P._msh.begin_quadrant_sweep (); 
			quadrant != P._msh.end_quadrant_sweep (); 
			++quadrant)
	{
		for( int i=0; i<4; i++ ){
			index = quadrant->gt(i);
			if(scnodes[index] == 1){
				beta[index] = 	V[index]/Vth ;
			}
		}
	}

	bim2a_advection_diffusion (	P._msh, alpha, beta, A22);
	
	for(unsigned i=0; i<ones.size(); i++){
		ones[i] = ones[i]/deltat;
	}
	
	bim2a_reaction (P._msh, not_ins, ones, Aa);

	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 1){
				A22[i][i] += Aa[i][i];
		}
	}
	
	resn = A22 * n;

	bim2a_rhs (P._msh, not_ins, ones, rhs);
	
	// Avoid cancellation errors.
	
	for(unsigned i=0; i<rhs.size(); i++){
		resn[i] -= rhs[i]*n0[i];	
	}
	rhs.clear();
	
	// ///	ENFORCING BCs ON SECOND ROW
	
	// std::vector<double>	rho,	// forse posso non specificare le dimensioni...
						// nimposed,
						// vec(n.size(),0.0);
	
	// // out.resize(nelements);		forse non è necessario
	// org_gaussian_charge_n(V, P, rho, out);	/// fare overload del metodo!!!
	// out.clear();
	
	// nimposed = rho;
	// for(unsigned i=0; i<nimposed.size(); i++){
		// nimposed[i] *= (-1)/q;
	// }
	
	// BCbulk.clear();
	// BCbulk.resize(n.size());
	
	// for(unsigned i=0; i<n.size(); i++){
		// BCbulk[i] = (n[i] - nimposed[i]);
	// }
	// BCbulk = M*BCbulk;
	
	// std::tuple<int, int, func_quad>	tuplan(0,2,[&resn,&BCbulk](tmesh::quadrant_iterator quad, tmesh::idx_t i)
																// {return (resn[quad->gt(i)]+BCbulk[quad->gt(i)]);});
	// dirichlet_bcs_quad	bcsn;
	// bcsn.push_back(tupla1);
	
	// bim2a_dirichlet_bc (P._msh,bcsn,M,resn);	/// che matrice metto qui ??
	
	/// ADJUST FOR ZERO INSULATOR CHARGE
	vec = M*n;
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 0){
			resn[i] = vec[i];
		}
	}
	
	for(unsigned i=0; i<indexingn.size(); i++){
		res[indexingn[i]] = resn[i];
	}
	resn.clear();
	
	
	///	ASSEMBLING THIRD ROW
	sparse_matrix	A, r;
	
	std::vector<double>	diff(F.size(),0),
						resF(F.size(),0),
						C;
	double	s;
	A.resize(4);
	r.resize(4);
	bcs.get_A(A);
	bcs.get_r(r); 
	
	C = bcs.get_C();
	
	for(unsigned i=0; i<F.size(); i++){
		diff[i] = (F[i]-F0[i])/deltat;
	}
	resF = A*diff;
	diff.clear();
	
	for(unsigned i=0; i<C.size(); i++){
		resF[i] += C[i];
	}

	for(unsigned i=0; i<r.size(); i++){
		s = 0;
		for(unsigned j=0; j<I.size(); j++){
			s += r[i][j]*I[j];
		}
		resF[i] += s;
	}
	
	for(unsigned i=0; i<indexingF.size(); i++){
		res[indexingF[i]] = resF[i];
	}
	resF.clear();
	C.clear();

	///	COMPUTING FOURTH ROW
	std::vector<int>	rr;
	
	sum1.clear();
	sum1.resize(A11.rows());
	
	sum2.clear();
	sum2.resize(A12.rows());
	
	for(unsigned i=0; i<indexingI.size(); i++){
		res[indexingI[i]] = I[i];
	}
	
	for(unsigned i=0; i<numcontacts; i++){
		rr.resize(dnodes[pins[i]].size());
		rr = dnodes[pins[i]];

		// Displacement current. 
		diff.clear();
		diff.resize(V.size());
		for(unsigned j=0; j<V.size(); j++){
			diff[j] = (V[j]-V0[j])/deltat;
		}
		for(unsigned j=0; j<rr.size(); j++){
			for(unsigned k=0; k<A11[rr[j]].size(); k++){
				sum1[k] += A11[rr[j]][k];
			}
		}
		s=0;
		for(unsigned j=0; j<V.size(); j++){
			s += section*sum1[j]*diff[j];
		}
		res[indexingI[i]] -= s;

		diff.clear();
		diff.resize(n.size());
		for(unsigned j=0; j<n.size(); j++){
			diff[j] = (n[j]-n0[j])/deltat;
		}
		for(unsigned j=0; j<rr.size(); j++){			
			for(unsigned k=0; k<A12[rr[j]].size(); k++){
				sum2[k] += A12[rr[j]][k];
			}
		}
		s=0;
		for(unsigned j=0; j<n.size(); j++){
			s += section*sum2[j]*diff[j];
		}
		res[indexingI[i]] -= s;

		// Electron current.
		s=0;
		for(unsigned j=0; j<rr.size(); j++){
			s += section * q * r22[rr[j]];
		}
		res[indexingI[i]] += s;
	}	
	rr.clear();
	sum1.clear();
	sum2.clear();
	diff.clear();
	r22.clear();
	
	for(unsigned i=0; i<indexingV.size(); i++){
		res[indexingV[i]] /= rowscaling[0];
	}
	for(unsigned i=0; i<indexingn.size(); i++){
		res[indexingn[i]] /= rowscaling[1];
	}
	for(unsigned i=0; i<indexingF.size(); i++){
		res[indexingF[i]] /= rowscaling[2];
	}
	for(unsigned i=0; i<indexingI.size(); i++){
		res[indexingI[i]] /= rowscaling[3];
	}
	
	return res;	
 };

 
 /**
 *	Assemble jacobian matrix: the matrix is passed to the method as reference,
 *	then it is updated during the computation
 */				 
void 
Newton::org_secs2d_newton_jacobian(	Probl& P, std::vector<double>& V, std::vector<double>& n, std::vector<double>& F,			
									double deltat, BCS_CIRC& bcs, std::vector<int>& indexingV, std::vector<int>& indexingn,
									std::vector<int>& indexingF, std::vector<int>& indexingI,sparse_matrix& jacobian)
{
	int	nnodes = P.get_msh_nodes(),
		nelements = P.get_msh_elem(),
        numcontacts, ndofs, numscnodes = 0, j = 0;
		
	double	eps_semic = P._eps_semic,
			eps_ins = P._eps_ins,
			q = P._q,
			Vth = P._Vth,
			section = P._section,
			L = P._L,
			tins = P._t_ins,
			tsemic = P._t_semic;
			
	bool	ins = P._ins;
			
	std::vector<int>::iterator it;
	std::vector< std::vector<int> >	dnodes = P._dnodes;
	//std::vector< std::vector<double> >	r;
	std::vector<double>	epsilon(nelements,eps_semic),
						rowscaling = P._rowscaling,
						colscaling = P._colscaling;
	std::array<int,2>	pins = P._pins;
	std::vector<int>	//ppins(pins.size(),0),
						scnodes = P._scnodes,
						//alldnodes = P._alldnodes,
						insulator = P._insulator;
						//intnodes(2*nnodes,0),
						//intdofsV,
						//intdofsn;
	sparse_matrix 	A, B, r;
	
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i]==1)
			numscnodes++;
	}
	
	// for(unsigned i=0; i<pins.size(); i++){
		// ppins[i] = pins[i];
	// }
	// it = std::unique (ppins.begin(), ppins.end());	
	// ppins.resize( std::distance(ppins.begin(),it) );
	// assert (ppins.size() == 2);
	// ppins.clear();
	
	numcontacts = pins.size();
	ndofs = 2 * nnodes + F.size() + numcontacts;
	
	// j=0;
	// for(int i=0; i<nnodes; i++){
		// it = std::find (alldnodes.begin(), alldnodes.end(), i);
		// if (it == alldnodes.end()){
			// intnodes[j] = i;			 // Assembling vector of internal nodes.
			// j++;
		// }
	// }
	// intnodes.resize(j);
	// std::sort(intnodes.begin(),intnodes.end());
	
	// intdofsV.resize(intnodes.size());
	// intdofsn.resize(intnodes.size());
	
	// for(unsigned i=0; i<intnodes.size(); i++){
		// intdofsV[i] = indexingV[intnodes[i]];
		// intdofsn[i] = indexingn[intnodes[i]];
	// }
	

	///	COMPUTING COEFFICIENTS	
	if(ins){
		for(unsigned i=0; i<insulator.size(); i++){	
			if(insulator[i] == 1){
				epsilon[i] = eps_ins;
			}
		}
	}
	

	/// COMPUTING FIRST ROW
	sparse_matrix	M,
					A11,
					A12;

	std::vector<double> eones(insulator.size(),1.0),
						not_ins(insulator.size(),0.0),
						ones(nnodes,1.0),
						psi(nnodes,0.0);
	
	M.resize(nnodes);
	A11.resize(nnodes);
	A12.resize(nnodes);
	
	for(unsigned i=0; i<insulator.size(); i++){
		if(insulator[i] == 0){
			not_ins[i] = 1;
		}
	}

	bim2a_reaction (P._msh, not_ins, ones, A12);
	for(int i=0; i<nnodes; i++){
		A12[i][i] *= q ;			
	}
	
	bim2a_reaction (P._msh, eones, ones, M);

	bim2a_advection_diffusion (P._msh, epsilon, psi, A11);

	///	Physical models.
	std::vector<double>	mobn(nelements,0.0),
						alphan(n.size(),0.0),
						der_dalpha_n(nelements,0.0);

	org_physical_models2d(n, P, mobn, alphan, der_dalpha_n);
  
	jacobian.resize(ndofs);

	///	ASSEMBLING FIRST ROW
	// sparse_matrix::col_iterator J2;
	
	// for(unsigned i=0; i<intdofsV.size(); i++){
		
		// J = A11[intnodes[i]].begin ();
		// J2 = A12[intnodes[i]].begin ();
		
		// for(unsigned j=0; j<indexingV.size(); j++){	
			// if( A11.col_idx(J) == j ){
				// jacobian[intdofsV[i]][indexingV[j]] += A11[intnodes[i]][j];	
				// J++;
			// }
			// if( A12.col_idx(J2) == j ){
				// jacobian[intdofsV[i]][indexingn[j]] += A12[intnodes[i]][j];
				// J2++;
			// }
		// }
	// }
	

	// ///	ENFORCING BCs ON FIRST ROW
	// std::vector<int>	ddofsV;
	// std::vector<double>	vec;
	// int	pindof = 0;

	// for(int i=0; i<numcontacts; i++){
		// ddofsV.resize(dnodes[i].size());
		// vec.resize(dnodes[i].size());
		// for(unsigned j=0; j<dnodes[i].size(); j++){
			// ddofsV[j] = indexingV[dnodes[i][j]] ;
			// vec[j] = M[ dnodes[i][j] ][ dnodes[i][j] ] ;
		// }	
		
		// pindof = indexingF[pins[i]];
		// for(unsigned j=0; j<ddofsV.size(); j++){
			// jacobian[ ddofsV[j] ][ ddofsV[j] ] += vec[j];
			// jacobian[ ddofsV[j] ][pindof] -= vec[j];
		// }
		// ddofsV.clear();
		// vec.clear();
	// }
	

	///	COMPUTING SECOND ROW

	sparse_matrix	A21,
					A22,
					R;
					
	std::vector<double>	fluxn(n.size(),0.0),
						alfa(insulator.size(),0.0),
						mobn_n(insulator.size(),0.0);
			
	A21.resize(nnodes);
	A22.resize(nnodes);
	R.resize(nnodes);	

	j=0;
	for (auto quadrant = P._msh.begin_quadrant_sweep ();
        quadrant != P._msh.end_quadrant_sweep ();
        ++quadrant)
	{		
		if(insulator[j]==1){
			// Adjust for zero insulator charge.
			mobn_n[j] = 0;
		}
		else{	

			for(int i=0; i<4; i++){
				mobn_n[j] += n[ quadrant->gt(i) ];
			}		
			mobn_n[j] /= 4;
		}
		j++;
	}
	
	for(unsigned i=0; i<insulator.size(); i++){
		mobn_n[i] *= -mobn[i];
		if(insulator[i] == 0){
			alfa[i] = mobn[i]*Vth;
		}
	}
	
	bim2a_advection_diffusion (P._msh, mobn_n, psi, A21);
		
	if(alphan.size() == n.size() && alphan.size() == V.size()){
		fluxn = alphan;
		for(unsigned i=0; i<fluxn.size(); i++){
			fluxn[i] *= (-1);
			fluxn[i] += V[i]/Vth;
		}
	
	}
	else{
		std::cout<<"error: org_secs2d_newton_jacobian, dimensions mismatch"<<std::endl;
	}
	
	bim2a_advection_diffusion( P._msh, alfa, fluxn, A22); 
	
	for(unsigned i=0; i<ones.size(); i++){
		ones[i] = ones[i]/deltat;
	}
	
	bim2a_reaction (P._msh, not_ins, ones, R);

	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i]==1){
					A22[i][i] += R[i][i];
		}
	}
	
	///	ENFORCING BCs ON SECOND ROW: both JAC and RES
	
	std::vector<double>	resn(nnodes,0.0);
	
	for(unsigned i=0; i<indexingn.size(); i++){
		resn[i] = _res[indexingn[i]];
	}
	
	std::vector<double>	rho,	// forse posso non specificare le dimensioni...
						nimposed,
						vec(n.size(),0.0);
	
	// out.resize(nelements);		forse non è necessario
	org_gaussian_charge_n(V, P, rho, out);	/// fare overload del metodo!!!
	out.clear();
	
	nimposed = rho;
	for(unsigned i=0; i<nimposed.size(); i++){
		nimposed[i] *= (-1)/q;
	}
	
	BCbulk.clear();
	BCbulk.resize(n.size());
	
	for(unsigned i=0; i<n.size(); i++){
		BCbulk[i] = (n[i] - nimposed[i]);
	}
	BCbulk = M*BCbulk;
	
	std::tuple<int, int, func_quad>	tuplan(0,2,[&resn,&BCbulk](tmesh::quadrant_iterator quad, tmesh::idx_t i)
																{return (resn[quad->gt(i)]+BCbulk[quad->gt(i)]);});
	dirichlet_bcs_quad	bcsn;
	bcsn.push_back(tupla1);
	
	bim2a_dirichlet_bc (P._msh, bcsn, A22, resn);

	///	ADJUST FOR ZERO INSULATOR CHARGE: -----> anche sia per jac che per res?
	std::vector<int>	insn(scnodes.size(),0);
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 0){
			insn[i] = 1;
		}
	}
	
	vec.clear();
	vec.resize(insn.size());
	j = 0;
	for(unsigned i=0; i<insn.size(); i++){
		if(insn[i] == 1){
			vec[j] = M[i][i];
			j++;
		}
	}
	vec.resize(j);
	
	sparse_matrix diag;
	diag.resize(vec.size());
	for(unsigned j=0; j<vec.size(); j++){
		diag[j][j] = vec[j];
	}
	
	int h = 0, i = 0;
	for(unsigned j=0; j<insn.size(); j++){
		if(insn[j] == 1){
			for (J = jacobian[indexingn[j]].begin (); J != jacobian[indexingn[j]].end (); ++J){
				jacobian[ indexingn[j] ][ jacobian.col_idx (J) ] = 0;
			}
			for(unsigned k=0; k<insn.size(); k++){
				if(insn[k] == 1 && i==h){
					jacobian[indexingn[j]][indexingn[k]] += diag[i][h];
					h++;
				}
			}
			i++;
		}
	}
	
	///	ASSEMBLING SECOND ROW
	
	// res
	for(unsigned i=0; i<indexingn.size(); i++){
		_res[indexingn[i]] = resn[i];
	}
	
	// jacobian
	for(unsigned i=0; i<indexingn.size(); i++){
		
		J = A21[i].begin ();
		J2 = A22[i].begin ();
		
		for(unsigned j=0; j<indexingV.size(); j++){
			if( A21.col_idx (J)==j ){
				jacobian[indexingn[i]][indexingV[j]] += A21[i][j];	
				J++;
			}
			if( A22.col_idx (J2)==j ){
				jacobian[indexingn[i]][indexingn[j]] += A22[i][j];
				J2++;
			}
		}
	}

	
	///	ASSEMBLING THIRD ROW
	A.resize(4);
	bcs.get_A(A);
	
	B.resize(4);
	bcs.get_B(B);
	
	r.resize(4);
	bcs.get_r(r);
	
	for(unsigned i=0; i<indexingF.size(); i++){
	
		J = B[i].begin();
	
		for(unsigned j=0; j<indexingF.size(); j++){
			jacobian[indexingF[i]][indexingF[j]] = (A[i][j]/deltat);
			if( B.col_idx(J) == j ){
				jacobian[indexingF[i]][indexingF[j]] += B[i][j];
				J++;
			}
		}
		for(unsigned j=0; j<indexingI.size(); j++){
			jacobian[indexingF[i]][indexingI[j]] = r[i][j];
		}
	}
	

	///	COMPUTING FOURTH ROW
	std::vector<int>	rr;
	std::vector<double>	s1(indexingV.size(),0),
						s2(indexingV.size(),0),
						zeros(indexingV.size(),0);
	sparse_matrix	eye;
	eye.resize(indexingI.size());
	for(unsigned i=0; i<indexingI.size(); i++){
		eye[i][i] = 1;
		for(unsigned j=0; j<indexingI.size(); j++){
			jacobian[indexingI[i]][indexingI[j]] = eye[i][j];
		}
	}
	
	for (int i=0; i<numcontacts; i++){
		rr.resize(dnodes[pins[i]].size());
		for(unsigned j=0; j<dnodes[pins[i]].size(); j++){
			rr[j] = dnodes[pins[i]][j] ;
		}
		s1 = zeros;
		s2 = zeros;
		for(unsigned k=0; k<rr.size(); k++){
			for(unsigned j=0; j<indexingV.size(); j++){
				s1[j] += (section*(A11[rr[k]][j])/deltat);
				s2[j] += (section*(A12[rr[k]][j])/deltat);
			}
		}
	
		// Displacement current.
		for(unsigned j=0; j<indexingV.size(); j++){
			jacobian[indexingI[i]][indexingV[j]] -= s1[j];
		}
		for(unsigned j=0; j<indexingn.size(); j++){
			jacobian[indexingI[i]][indexingn[j]] -= s2[j];
		}
		
		// Electron current.
		s1 = zeros;
		s2 = zeros;
		for(unsigned k=0; k<rr.size(); k++){
			for(unsigned j=0; j<indexingV.size(); j++){
				s1[j] += (-section * q * (A21[rr[k]][j]));
				s2[j] += (-section * q * (A22[rr[k]][j]));
			}
		}
		
		for(unsigned j=0; j<indexingV.size(); j++){
			jacobian[indexingI[i]][indexingV[j]] -= s1[j];
		}
		for(unsigned j=0; j<indexingn.size(); j++){
			jacobian[indexingI[i]][indexingn[j]] -= s2[j];
		}
		rr.clear();
	}
	zeros.clear();
	s1.clear();
	s2.clear();

		
	for(unsigned i=0; i<indexingV.size(); i++){
		for (J = jacobian[indexingV[i]].begin (); J != jacobian[indexingV[i]].end (); ++J){
			jacobian[indexingV[i]][jacobian.col_idx (J)] /= rowscaling[0];					
		}
	}
	for(unsigned i=0; i<indexingn.size(); i++){
		for (J = jacobian[indexingn[i]].begin (); J != jacobian[indexingn[i]].end (); ++J){
			jacobian[indexingn[i]][jacobian.col_idx (J)] /= rowscaling[1];					
		}
	}
	for(unsigned i=0; i<indexingF.size(); i++){
		for (J = jacobian[indexingF[i]].begin (); J != jacobian[indexingF[i]].end (); ++J){
			jacobian[indexingF[i]][jacobian.col_idx (J)] /= rowscaling[2];					
		}
	}
	for(unsigned i=0; i<indexingI.size(); i++){
		for (J = jacobian[indexingI[i]].begin (); J != jacobian[indexingI[i]].end (); ++J){
			jacobian[indexingI[i]][jacobian.col_idx (J)] /= rowscaling[3];					
		}
	}
	
 
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		for(unsigned j=0; j<indexingV.size(); j++){
			if( jacobian.col_idx(J) == indexingV[j] ){
				jacobian[i][indexingV[j]] *= colscaling[0];
				J++;
			}
		}
	}
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		for(unsigned j=0; j<indexingn.size(); j++){
			if( jacobian.col_idx(J) == indexingn[j] ){
				jacobian[i][indexingn[j]] *= colscaling[1];
				J++;
			}
		}
	}
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		for(unsigned j=0; j<indexingF.size(); j++){
			if( jacobian.col_idx(J) == indexingF[j] ){
				jacobian[i][indexingF[j]] *= colscaling[2];
				J++;
			}
		}
	}
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		for(unsigned j=0; j<indexingI.size(); j++){
			if( jacobian.col_idx(J) == indexingI[j] ){
				jacobian[i][indexingI[j]] *= colscaling[3];
				J++;
			}
		}
	}

};
 
 
void
Newton::org_physical_models2d (	std::vector<double>& n, Probl& P,
								// output
								std::vector<double>& mobilityn, std::vector<double>& alphan, std::vector<double>& der_dalpha_n)
 {
	std::vector<double>	nm,
						out;
	int vpe = 4;	//3;
	int nelements = P.get_msh_elem();
	double	s = 0, j = 0, i;

	nm.resize(nelements);
    for (auto quadrant = P._msh.begin_quadrant_sweep (); quadrant != P._msh.end_quadrant_sweep (); ++quadrant){
		s = 0;
		for (int ii = 0; ii < 4; ++ii){
			i= n[	quadrant->gt(ii) ];
			s += 1/i;
        }
		nm[j] = (vpe/s);
		j++;
    }	

	// Mobility model.
	mobilityn = org_secs_mobility_EGDM (P._mu0n, nm, P._N0, P._sigman_kT, P);
	
	// Generalized Einstein relation.
	out.resize(n.size());
	org_einstein_n(n, P, alphan, out);	// output: alphan, Piece-wise linear.
	out.clear();	
	
	// Generalized Einstein relation.
    org_einstein_n (nm, P, out, der_dalpha_n); // output: der_dalpha_n, Element-wise constant.
	out.clear();

};

std::vector<double> 
Newton::org_secs_mobility_EGDM (double mu0, std::vector<double>& c, double C0, double sigma, Probl& P)
{
	for(unsigned i=0; i<c.size(); i++){
		c[i] = std::min(c[i], 0.1 * C0);
	}
	int j = 0;
	std::vector<double> data_n = P._data_n,
						udata_n, udata_phi_lumo, out, mu;
	std::vector<int>	ind;

	std::vector<double>	phi_lumo(c.size(),0),
						g1;

	double	q, Vth, Kb, T0, Efield;
	q = P._q;
	Vth = P._Vth;
	Kb = P._Kb;
	T0 = P._T0;
	Efield = P._Efield;

	for(unsigned i=0; i<data_n.size(); i++){
		data_n[i] *= (-1)/q;
	}
	if (udata_n.size()==0 || ind.size()==0 || udata_phi_lumo.size()==0){

		udata_n = data_n;
		std::vector<double>::iterator it;
		it = std::unique (udata_n.begin(), udata_n.end());			///lento
		udata_n.resize( std::distance(udata_n.begin(),it) );		///lento
		std::sort(udata_n.begin(), udata_n.end());					///lento
		j = 0;
		ind.resize(udata_n.size());
		for(unsigned i=udata_n.size()-1; i>0; i--){
			if(i==1){
				if(udata_n[i]!=udata_n[i-1]){
					ind[j] = i;
					ind[j+1] = 0;
					j += 2;
				}
				else{
					ind[j] = 0;
					j++;
				}
			}
			else{
				if(udata_n[i]!=udata_n[i-1]){
					ind[j] = i;
					j++;
				}
			}
		}
		ind.resize(j);
		std::sort(ind.begin(), ind.end());
		
		udata_phi_lumo.resize(ind.size());
		for(unsigned i=0; i<ind.size(); i++){
			udata_phi_lumo[i] = P.get_data_phi_lumo()[ ind[i] ] ;
		}
	}
	for(unsigned i=0; i<c.size(); i++){
		phi_lumo[i] = interp1( udata_n, udata_phi_lumo, c[i], true );
	}
	
	int phi_lumo_ref = -2;
	double n_ref;

	n_ref = org_gaussian_charge_n( phi_lumo_ref, P.mat(), P.cnst(), P.quad());
	n_ref = -n_ref/q;
  
	double g1_ref = exp (phi_lumo_ref / Vth + 0.5 * pow((sigma * Kb * T0),2)) / (n_ref / C0);

	// Density enhancement factor.
	g1.resize(phi_lumo.size());
	for(unsigned i=0; i<phi_lumo.size(); i++){
		if(phi_lumo.size() != c.size()){
			std::cout<<"error: org_secs_mobility_EGDM, dimensions mismatch"<<std::endl;
			break;
		}
		if(c[i] == 0){
			g1[i] = 0;
		}
		else{
			g1[i] = ( exp (phi_lumo[i] / Vth + 0.5 * pow((sigma * Kb * T0),2)) / (c[i] / C0) ) / g1_ref;
		}
	}

	// Field enhancement factor.
	double	g2;
	mu.resize(g1.size());
			
	g2 = exp((0.44*(pow(sigma,1.5)-2.2))*(sqrt(1+0.8*pow(std::fmin(q*1/std::cbrt(C0)*Efield/(sigma*Kb*T0),2),2))-1));

	mu = g1;
	for(unsigned i=0; i<g1.size(); i++){
		mu[i] *= mu0*g2;
	}
	g1.clear();
	
	return mu;
};

void
Newton::org_gaussian_charge_n(	std::vector<double>& V, Probl& P,
								// output
								std::vector<double>& rhon, std::vector<double>& drhon_dV)
{
	double	q = P._q;
    std::vector<double> n = n_approx(V,P);
	
	rhon = n;
    for(unsigned i=0; i<n.size(); i++){
		rhon[i] *= -q;
    }

    std::vector<double> dn_dV = dn_dV_approx(V,P);
	
	drhon_dV = dn_dV;
    for(unsigned i=0; i<dn_dV.size(); i++){
		drhon_dV[i] *= -q;
    }
};

double
Newton::org_gaussian_charge_n( double V, Probl& P)
{
	double	q = P._q,
			n, rhon;
    n = n_approx(V,P);
    rhon = -q*n;
	return rhon;
};

std::vector<double> 
Newton::n_approx(std::vector<double>& V, Probl& P)
{
    std::vector<double> coeff(V.size(),0),
						n(V.size(),0);
    double	q = P._q,
			kT = P._Kb*P._T0,
			sigman = P._sigman,
			N0 = P._N0,
			denom;

	std::vector<double>	gx = P._gx;
	std::vector<double>	gw = P._gw;
	
    for(unsigned i=0; i<gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){		
            coeff[j] = (sqrt(2) * sigman * gx[i] - q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            n[j] += N0 / sqrt(M_PI) * gw[i] / denom;
        }
    }
	gx.clear();
	gw.clear();
	coeff.clear();	
	
	return n;
};

double
Newton::n_approx( double V, Probl& P)
{
    double	coeff, n = 0, denom,
			q = P._q,
			kT = P._Kb*P._T0,
			sigman = P._sigman,
			N0 = P._N0;

	std::vector<double>	gx = P._gx,
						gw = P._gw;	

    for(unsigned i=0; i<gx.size(); i++){
            coeff = (sqrt(2) * sigman * gx[i] - q * V) / kT ;
            denom = 1+exp(coeff);
            n += N0 / sqrt(M_PI) * gw[i] / denom;
    }
	gx.clear();
	gw.clear();
	return n;
};

std::vector<double>
Newton::dn_dV_approx(std::vector<double>& V, Probl& P)
{
    std::vector<double> coeff(V.size(),0),
						dn_dV(V.size(),0);
    double	kT = P._Kb*P._T0,
					sigman = P._sigman,
					N0 = P._N0,
					q = P._q;
	double	denom;
	std::vector<double>	gx = P._gx;
	std::vector<double>	gw = P._gw;

    for(unsigned i=0; i<gx.size(); i++){
        for(unsigned j=0; j<V.size(); j++){
            coeff[j] = (sqrt(2) * sigman * gx[i] - q * V[j]) / kT ;
            denom = 1+exp(coeff[j]);
            dn_dV[j] += - q * N0 / sigman * sqrt(2/M_PI) * gw[i]*gx[i] / denom;
		}
    }
	gx.clear();
	gw.clear();
	coeff.clear();
	return dn_dV;
};

void
Newton::org_einstein_n (std::vector<double>& n, Probl& P,
						// output
						std::vector<double>& alphan, std::vector<double>& dalphan_dn)
{
	std::vector<double> phi,
						data_n,
						udata_n,
						ind,
						udata_alphan,
						out,
						udata_dalphan_dn,
						data_alphan,
						data_dalphan_dn;

	std::vector<int> cutoff;
	double	q = P._q,
			N0 = P._N0;
	
	int j = 0;
	double	a, b, N, step;
	a = -10;
	b = 10;
	N = 1001;
	step = (b-a)/(N-1);
	phi.resize(N);
	j = 0;
	for(auto i=a; i<=b; i+=step){
		phi[j] = i;
		j++;
	}
	// phi = 0 corresponds to n = N0 / 2.

	data_n.resize(phi.size());
	out.resize(phi.size());
	org_gaussian_charge_n( phi,P,data_n, out);
	out.clear();
	for(unsigned i=0; i<data_n.size(); i++){
		data_n[i] *= (-1)/ q;
	}
	
	data_alphan = alphan_fun(phi,P);
	
	data_dalphan_dn = dalphan_dn_fun (phi,data_alphan,P);	

	// Extract unique data to improve robustness.
  if (udata_n.size() == 0 || ind.size() == 0 || udata_alphan.size() == 0 || udata_dalphan_dn.size() == 0){

	udata_n = data_n;
	std::vector<double>::iterator it;
	it = std::unique (udata_n.begin(), udata_n.end());
	udata_n.resize( std::distance(udata_n.begin(),it) );
	std::sort(udata_n.begin(), udata_n.end());

	ind.resize(udata_n.size());
	j = 0;
	for(unsigned i=udata_n.size()-1; i>0; --i){
		if(i==1){
			if(udata_n[i]!=udata_n[i-1]){
				ind[j] = i;
				ind[j+1] = 0;
				j += 2;
			}
			else{
				ind[j] = 0;
				j++;
			}
		}
		else{
			if(udata_n[i]!=udata_n[i-1]){
				ind[j] = i;
				j++;
			}
		}
	}
	ind.resize(j);
	std::sort(ind.begin(), ind.end());

	udata_alphan.resize(ind.size());
	udata_dalphan_dn.resize(ind.size());
    for (unsigned i=0; i<ind.size(); i++){
        udata_alphan[i] = data_alphan[ ind[i] ];
        udata_dalphan_dn[i] = data_dalphan_dn[ ind[i] ];
    }

    // Fix for asymptotic values.
	j = 0;
	cutoff.resize(udata_alphan.size());
    for (unsigned i=0; i<udata_alphan.size(); i++){
        if(udata_alphan[i] <= 1){
            cutoff[j] = i;
			j++;
        }
    }
	cutoff.resize(j);
    for(unsigned i=0; i<cutoff.size(); i++){
        udata_alphan[cutoff[i]] = udata_alphan[cutoff.back()];
        udata_dalphan_dn[cutoff[i]] = 0;
    }

    cutoff.clear();
	cutoff.resize(udata_n.size());
	j = 0;
    for (unsigned i=0; i<udata_n.size(); i++){
        if(udata_n[i] <= 1e5){
            cutoff[j] = i;
			j++;
        }
    }
	cutoff.resize(j);
    for(unsigned i=0; i<cutoff.size(); i++){
        udata_alphan[cutoff[i]] = udata_alphan[cutoff.back()];
        udata_dalphan_dn[cutoff[i]] = 0;
    }

    cutoff.clear();
	cutoff.resize(udata_n.size());
	j = 0;
    for (unsigned i=0; i<udata_n.size(); i++){
        if(udata_n[i] >= 0.99 * N0){
            cutoff[j] = i;
			j++;
        }
    }
	cutoff.resize(j);

    for(unsigned i=0; i<cutoff.size(); i++){
        udata_alphan[cutoff[i]] = udata_alphan[cutoff[0]];
        udata_dalphan_dn[cutoff[i]] = 0;
    }
	cutoff.clear();
  }

	// Compute alphan by interpolation.
	for(unsigned i=0; i<n.size(); i++){
		alphan[i] = interp1( udata_n, udata_alphan, n[i], true );
		dalphan_dn[i] = interp1( udata_n, udata_dalphan_dn, n[i], true );
	}
};

std::vector<double>
Newton::alphan_fun (std::vector<double>& phi, Probl& P)
{
	std::vector<double> n(phi.size(),0),
						dn_dphi,
						alphan(phi.size(),0),
						out(phi.size(),0);
	double	q = P._q,
			Vth = P._Vth;
  
	org_gaussian_charge_n(phi,P,n,out);
	for(unsigned i=0; i<n.size(); i++){
		n[i] *= (-1)/q;
	}

	dn_dphi = dn_dphi_approx (phi,P);
	
	for(unsigned i=0; i<n.size(); i++){
		alphan[i] = (1 / Vth) * n[i]/dn_dphi[i];
	}
	n.clear();
	dn_dphi.clear();
	out.clear();

	return alphan;
};

std::vector<double>
Newton::dn_dphi_approx (std::vector<double>& phi, Probl& P)
{						
	const double  	N0 		= P._N0,
					sigman	= P._sigman,
					kT		= P._Kb * P._T0,
					q		= P._q,
					Vth		= P._Vth;

	std::vector<double>	coeff(phi.size(),0),
						dn_dphi(phi.size(),0),
						gx = P._gx,
						gw = P._gw;

	for (unsigned i=0; i<gx.size(); i++){
		for(unsigned j=0; j<phi.size(); j++){
			coeff[j] = (sqrt(2) * sigman * gx[i] - q * phi[j]) / kT;

			dn_dphi[j] += N0 / (Vth * sqrt(M_PI)) * gw[i] * exp(coeff[j]) / pow((1 + exp( coeff[j] )),2);
		}
	}
	coeff.clear();
	gx.clear();
	gw.clear();
	return dn_dphi;
};

std::vector<double>
Newton::d2n_dphi2_approx (std::vector<double>& phi, Probl& P)
{
	const double	N0		= P._N0,
					sigman	= P._sigman,
					kT		= P._Kb * P._T0,
					q		= P._q,
					Vth		= P._Vth;
  
	std::vector<double>	coeff(phi.size(),0.0),
						d2n_dphi2(phi.size(),0.0),
						gx = P._gx,
						gw = P._gw;
	
	for(unsigned i=0; i<gx.size(); i++){
		for(unsigned j=0; j<phi.size(); j++){
			coeff[j] = (sqrt(2)*sigman*gx[i] - q*phi[j])/(kT);
			
			d2n_dphi2[j] += N0 / (pow(Vth,2)*sqrt(M_PI))*gw[i]*exp(coeff[j])/pow((1 + exp(coeff[j])),2) * //
							(2*exp(coeff[j])/(1+exp(coeff[j])) - 1);
		}
	}
	coeff.clear();
	gx.clear();
	gw.clear();
	return d2n_dphi2;
};

std::vector<double>
Newton::dalphan_dn_fun (std::vector<double>& phi, std::vector<double>& alphan, Probl& P)
{
	double	Vth = P._Vth;
	std::vector<double>	dn_dphi, d2n_dphi2,
						dalphan_dn(phi.size(),0);

	dn_dphi = dn_dphi_approx (phi,P);

	d2n_dphi2 = d2n_dphi2_approx (phi,P);
		
	for(unsigned i=0; i<phi.size(); i++){
		dalphan_dn[i] = (1 / Vth) * (1 / dn_dphi[i]) - alphan[i] * d2n_dphi2[i] / pow(dn_dphi[i],2);
	}
	return dalphan_dn;
};

void 
Newton::diff (sparse_matrix& in, double a, double b, std::vector<double>& out)
{	
	out.resize( in.size() );
    for(unsigned i=0; i<in.size(); i++){
        out[i] = (in[i][b] - in[i][a]);
    }
};

void
Newton::org_secs_state_predict (	Probl& P, sparse_matrix& V, sparse_matrix& n, sparse_matrix& F,
									sparse_matrix& I, int& tstep, std::vector<double>& tout,
									// output
									std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0)
{
  if (tstep > 2){

    double  dt, difft;
    std::vector<double> it(2,0), t(2,0),
						vdiff1, dndt, dlndt, dVdt, dIdt, dFdt;
	sparse_matrix  n_2, ln;

    it[0] = (tstep - 2);
    it[1] = (tstep - 1);

    t[0] = (tout[it[0]] - tout[tstep - 2]);
    t[1] = (tout[it[1]] - tout[tstep - 2]);
	
    dt = tout[tstep] - tout[tstep - 1];

	///diff(n,it[0],it[1], vdiff1);	 ???
	diff(n,0,1, vdiff1);
    difft = t[1]-t[0];
	dndt = vdiff1;
    for(unsigned i=0; i<dndt.size(); i++){
		dndt[i] /= difft;
    }
	vdiff1.clear();

    // n
	n_2.resize(P._scnodes.size());
    for(unsigned i=0; i<P._scnodes.size(); i++){
        if(P._scnodes[i] == 1){
			n_2[i][0] = n[it[0]][i];
			n_2[i][1] = n[it[1]][i];
        }
    }
	
    //ln = log (n(device.scnodes, it(1:2)));
	ln.resize(n_2.size());
    for(int j=0; j<2; j++){
        for(unsigned i=0; i<n_2.size(); i++){
			ln[i][j] = log(n_2[i][j]);
        }
    }

    //dlndt  = diff (ln, 1, 2) ./ diff (t);
    vdiff1.clear();
	diff(ln,0.0,1.0, vdiff1);
	
	dlndt = vdiff1;
	for(unsigned i=0; i<dlndt.size(); i++){
		dlndt[i] /= difft;
    }
	vdiff1.clear();
	
	//n0 = n[tstep-1];
	for(unsigned i=0; i<n.size(); i++){
		n0[i] = n[i][1];	// colonna tstep-1, indice 1 o 2 ????
	}

    // V, F, I
    //dVdt  = diff (V(:, it(1:2)), 1, 2) ./ diff (t);
    vdiff1.clear();
	///diff(V,it[0],it[1], vdiff1);
	diff(V,0,1, vdiff1);
	dVdt = vdiff1;
	for(unsigned i=0; i<dVdt.size(); i++){
		dVdt[i] /= difft; 
    }

    //dFdt  = diff (F(:, it(1:2)), 1, 2) ./ diff (t);
    vdiff1.clear();
	///diff(F,it[0],it[1], vdiff1);
	diff(F,0,1, vdiff1);
	dFdt = vdiff1;
	for(unsigned i=0; i<dFdt.size(); i++){
		dFdt[i] /= difft;
    }

    //dIdt  = diff (I(:, it(1:2)), 1, 2) ./ diff (t);
    vdiff1.clear();
	///diff(I,it[0],it[1], vdiff1);
	diff(I,0,1, vdiff1);
	dIdt = vdiff1;
	for(unsigned i=0; i<dIdt.size(); i++){
		dIdt[i] /= difft;
    }
	vdiff1.clear();

    //V0 = V(:, it(2)) + dVdt * dt;
	V0.resize(dVdt.size());
	for(unsigned i=0; i<dVdt.size(); i++){
		V0[i] = ( dVdt[i]*dt );
	}
	if(V0.size() != V.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<V.size(); i++){
			V0[i] += V[i][it[1]];
		}
	}

    //n0(device.scnodes) = exp (ln(:, 2) + dlndt * dt);
	for(unsigned i=0; i<dlndt.size(); i++){
		dlndt[i] *= dt ;
	}
	
	if(dlndt.size() != ln.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<dlndt.size(); i++){
			dlndt[i] += ln[i][1];
		}	
		for(unsigned i=0; i<P._scnodes.size(); i++){
			if(P._scnodes[i] == 1){
				n0[i] = std::exp(dlndt[i]);
			}
		}
	}

    //F0 = F(:, it(2)) + dFdt * dt;
	F0.resize(dFdt.size());
	for(unsigned i=0; i<dFdt.size(); i++){
		F0[i] = ( dFdt[i]*dt );	
	}
	if(F0.size() != F.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<F.size(); i++){
			F0[i] += F[i][it[1]];
		}
	}

    //I0 = I(:, it(2)) + dIdt * dt;
	I0.resize(dIdt.size());
	for(unsigned i=0; i<dIdt.size(); i++){
		I0[i] = ( dIdt[i]*dt );	
	}
	if(I0.size() != I.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<I.size(); i++){
			I0[i] += I[i][it[1]];
		}
	}
  }
  else{
  
	V0.resize(V.size());
	n0.resize(n.size());
	if(V0.size() != n0.size() ){
		std::cout<<"Error: dimensions mismatch, org_secs_state_predict, tstep<1"<<std::endl;
	}
	else{
	for(unsigned i=0; i<V.size(); i++){
		V0[i] = V[i][1];
		n0[i] = n[i][1];
	}
	F0.resize(F.size());
	for(unsigned i=0; i<F.size(); i++){
		F0[i] = F[i][1];
	}
	I0.resize(I.size());
	for(unsigned i=0; i<I.size(); i++){
		I0[i] = I[i][1];
	}
	}
  }
};

void 
Newton::compute_residual_norm (	double& resnrm, int& whichone, std::vector<double>& resall, std::vector<double>& res,
								std::vector<int>& idxV, std::vector<int>& idxn, std::vector<int>& idxF, std::vector<int>& idxI)
{
	std::vector<double>	aux(idxV.size(),0);
	
	for(unsigned i=0; i<idxV.size(); i++){
		aux[i] = std::abs( res[ idxV[i] ] );
	}
	resall[0] = *std::max_element( aux.begin(),aux.end() );
  	
	aux.resize(idxn.size());
	for(unsigned i=0; i<idxn.size(); i++){
		aux[i] = std::abs( res[idxn[i]] );
	} 
	resall[1] = *std::max_element( aux.begin(),aux.end() );  
  
	aux.resize(idxF.size());
	for(unsigned i=0; i<idxF.size(); i++){
		aux[i] = std::abs( res[idxF[i]] );
	}
	resall[2] = *std::max_element( aux.begin(),aux.end() );

	aux.resize(idxI.size());
	for(unsigned i=0; i<idxI.size(); i++){
		aux[i] = std::abs( res[idxI[i]] );
	}
	resall[3] = *std::max_element( aux.begin(),aux.end() );
  
	aux.clear();
	resnrm = *std::max_element(resall.begin(), resall.end());
		for(unsigned i=0; i<resall.size(); i++){
			if(resall[i] == resnrm){
				whichone = i;
				break;
			}
		}
};

void 
Newton::compute_variation (	std::vector<double>& Va, std::vector<double>& na, std::vector<double>& Fa, std::vector<double>& Ia, 
							std::vector<double>& Vb, std::vector<double>& nb, std::vector<double>& Fb, std::vector<double>& Ib, 
							Probl& P, double clamping,
							// output
							double& incrV, double& incrn, double& incrF, double& incrI)						
{
  std::vector<int>  cl = P._clamping,
                    sc = P._scnodes;
  std::vector<double>   diff(Vb.size(),0),
                        q1,q2;
  double    ni = P._ni,
			Vth = P._Vth;

  int j=0;

  if(Vb.size() != Va.size()){
    std::cout<<"error: compute_variation, different sizes"<<std::endl;
  }
  else{
    for(unsigned i=0; i<Vb.size(); i++){
        diff[i] = Vb[i]-Va[i];
    }
  }
  
	incrV = norm(diff,0) / (norm(Va,0) * clamping + cl[0]);

  // incrn = constants.Vth * norm (log (nb(sc) ./ na(sc)), inf) / ...
  //         (norm (Va(sc) - constants.Vth * log (na(sc) ./ ni), inf) * clamping + c[1]);
  diff.clear();
  q1.resize(sc.size());
  q2.resize(sc.size());
  j=0;
  for(unsigned i=0; i<sc.size(); i++){
    if(sc[i] == 1){
        q1[j] = log(nb[i]/na[i]);
        q2[j] = log(na[i]/ni);
		j++;
    }
  }
  q1.resize(j);
  q2.resize(j);
  
  j=0;
  diff.resize(sc.size());
  for(unsigned i=0; i<sc.size(); i++){
    if(sc[i] == 1){
        diff[j] = Va[i] - (Vth * q2[i]);
		j++;
    }
  }
  diff.resize(j);
  
	incrn = Vth * norm(q1,0) / (norm(diff,0) * clamping + cl[1]);	

	diff.clear();
	if(Fb.size() != Fa.size()){
		std::cout<<"error: compute_variation, different sizes"<<std::endl;
	}
	else{
		diff.resize(Fb.size());
		for(unsigned i=0; i<Fb.size(); i++){
			diff[i] = Fb[i]-Fa[i];
		}
	}
	incrF = norm(diff, 0) / (norm(Fa, 0) * clamping + cl[2]);
	
	diff.clear();
	if(Ib.size() != Ia.size()){
		std::cout<<"error: compute_variation, different size"<<std::endl;
	}
	else{
		diff.resize(Ib.size());
		for(unsigned i=0; i<Ib.size(); i++){
			diff[i] = Ib[i]-Ia[i];
		}
	}
	incrI = norm(diff, 0) / (norm(Ia, 0) * clamping + cl[3]);
	
	diff.clear();
};

void 
Newton::org_secs_safe_increment (	std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0,
									std::vector<double>& dV, std::vector<double>& dn, std::vector<double>& dF, std::vector<double>& dI,
									Probl& P,
									std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I,
									int& clamp, double& tauk)
{
	double tk = 1;
	std::vector<double>	Min, n0_2;
	std::vector<int>	scnodes = P._scnodes,
						clamping = P._clamping;
	
  
	std::vector<int>  where(n0.size(),0);

	if (n0.size() != dn.size()){
		std::cout<<"error: org_secs_safe_increment, vectors with different size"<<std::endl;
	}
	for(unsigned i=0; i<n0.size(); i++){
		if(n0[i]+dn[i] <= 0)
			where[i] = 1;
	}
	int s = 0 , length = 0;
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 1){
			s++;
			if(where[i] == 1){
				length++;
			}
		}
	}

	if ( length == s){
		Min.resize(length);
		for(unsigned i=0; i<scnodes.size(); i++){
			if(scnodes[i]==1){
				if(where[i]==1){
					Min[i] = n0[i]/std::abs(dn[i]);
				}
			}
		}
		tk = 0.9 * (*std::min_element(Min.begin(),Min.end()));
		Min.clear();
	}

  clamp = 1;

  if (P._clampOnOff){
    // V

    if (any(dV)){
		clamp = std::fmin (1, clamping[0] / norm(dV,0));
    }
	

    // n
    std::vector<double> dn2(dn.size(),0),
						dn3(dn.size(),0);
    for(unsigned i=0; i<dn.size(); i++){
        if(dn[i] > 0)
            dn2[i] = 1;
        else if(dn[i] < 0)
            dn3[i] = 1;
    }

    if (any (dn2)){
      for(unsigned i=0; i<dn.size(); i++){
        dn2[i] = dn2[i]*(dn[i]/n0[i]);
      }
      clamp = std::fmin (clamp, (exp (clamping[1] / P._Vth) - 1) / norm (dn2,0) );
    }

    if (any (dn3)){
      for(unsigned i=0; i<dn.size(); i++){
        dn3[i] = dn3[i]*(dn[i]/n0[i]);
      }
      clamp = std::fmin (clamp, (1 - exp (-clamping[1]/P._Vth)) / norm (dn3, 0) );
    }

    // F
    if (any(dF)){
      clamp = std::fmin (clamp, clamping[2] / norm (dF, 0));
    }

    // I
    if (any(dI)){
      clamp = std::fmin (clamp, clamping[3] / norm (dI, 0));
    }
  }
  else{
    std::cout<<"No clamping applied"<<std::endl;
  }
  
  tauk = tk;
  tk = std::fmin (tauk, clamp);


  if(V0.size() != dV.size()){
    std::cout<<"error: org_secs_safe_increment, dimensions mismatch (V)"<<std::endl;
  }
  else{
	V.resize(dV.size());
    for(unsigned i=0; i<dV.size(); i++){
        V[i] = (V0[i] + tk * dV[i]);
    }
  }

  
  if(n0.size() != dn.size()){
    std::cout<<"error: org_secs_safe_increment, dimensions mismatch (n)"<<std::endl;
  }
  else{
	n.resize(dn.size());
    for(unsigned i=0; i<dn.size(); i++){
        n[i] = (n0[i] + tk * dn[i]);
    }
  }

  
  if(F0.size() != dF.size()){
    std::cout<<"error: org_secs_safe_increment, dimensions mismatch (F)"<<std::endl;
  }
  else{
	F.resize(dF.size());
    for(unsigned i=0; i<dF.size(); i++){
        F[i] = (F0[i] + tk * dF[i]);
    }
  }

  //I = I0 + tk * dI;
  if(I0.size() != dI.size()){
    std::cout<<"error: org_secs_safe_increment, dimensions mismatch (I)"<<std::endl;
  }
  else{
	I.resize(dI.size());
    for(unsigned i=0; i<dI.size(); i++){
        I[i] = (I0[i] + tk * dI[i]);
    }
  }

  if (tk <= 0){
    std::cout<<"SECS: nonpositive-relax"<<std::endl;
    std::cout<<"Relaxation parameter is too small: tau = "<<tk<<std::endl;
    std::cout<<" reducing time step!"<<std::endl;
  }

  n0_2.resize(s);
  int j=0;
  for(unsigned i=0; i<scnodes.size(); i++){
	if(scnodes[i]==1){
		if(n[i] <= 0){
            n0_2[j] = 1;
			j++;
        }
		else{
            n0_2[j] = 0;
			j++;
		}
	}
  }
  
  //if (any (n(_scnodes) <= 0))
	if(any(n0_2)){
		std::cout<<"error: org_secs_safe_increment, negative charge density"<<std::endl;
		exit(EXIT_FAILURE);
	}	
};

void
Newton::CONV_MSG (int tstep, int Nstep, int mNstep, double t, std::string reason, int field, double incr, std::vector<double>& res)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep-1<<", modified Newton iteration "<<mNstep-1<<", model time "<<t   //
            <<": convergence reached ("<<reason<<")"<<std::endl;
  std::cout<<"incr ("<<field<<") = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }

  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};

void 
Newton::MAXINCR_MSG (int tstep, double t, int Nstep, int field, double incr, std::vector<double>& res, Probl& P)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<" model time "<<t<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep-1;
  std::cout<<" the increment in field "<<field<<" has grown too large: "<<std::endl;
  std::cout<<" "<<incr<<" > "<<P._maxnpincr<<std::endl;
  std::cout<<"incr0 ( "<<field<<" ) = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }

  std::cout<<"]"<<res.back()<<std::endl;
  std::cout<<" "<<std::endl;
};

void
Newton::DIV_MSG (	int tstep, double t, int Nstep, int field, std::vector<double>& incrhist, //
					double incr, std::vector<double>& res, int nsteps_check)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<" model time "<<t<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep<<", the Newton algorithm is diverging: "<<std::endl;
  std::cout<<"the increment in field "<<field<<" is not decreasing : "<<std::endl;
  std::cout<<incrhist[Nstep-1]<<" > "<<incrhist[Nstep - nsteps_check-1]<<std::endl;
  std::cout<<"incr ("<<field<<") = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }
  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};

void
Newton::DIV_MN_MSG (	int tstep, int t, int Nstep, int mNstep, int field, std::vector<double>& incrhist,	//
						double incr, std::vector<double>& res, int nsteps_check)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<"model time "<<t<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep<<", modified Newton iteration "<<mNstep<<std::endl;
  std::cout<<"the Modified Newton algorithm is diverging: "<<std::endl;
  std::cout<<"the increment in field "<<field<<" is not decreasing : "<<std::endl;
  std::cout<<incrhist[mNstep-1]<<" > "<<incrhist[mNstep - nsteps_check-1]<<std::endl;
  std::cout<<"incr ("<<field<<") = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }
  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};