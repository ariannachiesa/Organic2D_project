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
									double deltat, BCS_CIRC& bcs, ordering& ordV, ordering& ordn, std::vector<int>& indexingF, std::vector<int>& indexingI)
{
	int	nnodes = P.get_msh_nodes(),
		nelements = P.get_msh_elem(),
        ndofs;
	double	eps_semic = P._eps_semic,
			eps_ins = P._eps_ins,
			q = P._q,
			Vth = P._Vth,
			section = P._section,
			s = 0.0;
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
						alphan(n.size(),0.0);

	org_physical_models2d(n, P, mobn, alphan);
	
	///	ASSEMBLING FIRST ROW
	
	res.resize(ndofs);
	std::vector<double>	sum1(V.size(),0.0),
						sum2(n.size(),0.0),
						resV(V.size(),0.0),
						f(V.size(),0.0);

	sum1 = A11*V;	
	sum2 = A12*n;
	
	for(unsigned i=0; i<resV.size(); i++){
		if(sum1.size() != resV.size() || sum2.size() != resV.size()){
			std::cout<<"error: resV; wrong dimensions!"<<std::endl;
			exit(EXIT_FAILURE);
		}
		else{
			resV[i] = -sum1[i] + sum2[i];
		}
	}

	for(int i=0; i<nnodes; i++){
		res[ordV(i)] = resV[i];
	}
	resV.clear();
	
	///	COMPUTING SECOND ROW
	sparse_matrix	A22,
					R;
	std::vector<double>	resn(nnodes,0.0),
						rhs(nnodes, 0.0),
						alpha(insulator.size(),0.0),
						a_el(insulator.size(),0.0),
						eta(scnodes.size(),0.0),
						beta(nnodes,0.0),
						gamma(scnodes.size(),1.0),
						mob_nodes(scnodes.size(),0.0);
  
	A22.resize(nnodes);
	R.resize(nnodes);
		
	for(unsigned i=0; i<insulator.size(); i++){
		if(insulator[i] == 0){
			alpha[i] = mobn[i]*Vth;
		}
	}

	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 1){
			beta[i] = V[i]/Vth;
		}
	}

	bim2a_advection_diffusion (	P._msh, alpha, beta, A22);

	for(unsigned i=0; i<ones.size(); i++){
		ones[i] = ones[i]/deltat;
	}
	
	bim2a_reaction (P._msh, not_ins, ones, R);

	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 1){
				A22[i][i] += R[i][i];
		}
	}
	
	resn = A22 * n;

	bim2a_rhs (P._msh, not_ins, ones, rhs);
	
	// Avoid cancellation errors.
	
	for(unsigned i=0; i<rhs.size(); i++){
		resn[i] -= rhs[i]*n0[i];
	}
	rhs.clear();
	
	/// ADJUST FOR ZERO INSULATOR CHARGE	
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 0){
			resn[i] = 0;
		}
	}
	
	for(int i=0; i<nnodes; i++){
		res[ordn(i)] = resn[i];
	}
	
	
	///	ASSEMBLING THIRD ROW
	sparse_matrix	A, r;
	
	std::vector<double>	diff(F.size(),0),
						resF(F.size(),0),
						C;
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
			diff[j] = (V[j]-V0[j]);
		}
		sum1.resize(rr.size());
		for(unsigned j=0; j<rr.size(); j++){
			s = 0;
			for(unsigned k=0; k<A11[rr[j]].size(); k++){
				s += A11[rr[j]][k]*diff[k];
			}
			sum1[j] = s;
		}
		s = 0;
		for(unsigned j=0; j<sum1.size(); j++){
			s += sum1[j];
		}
		sum1.clear();
		res[indexingI[i]] -= (section * s / deltat);

		diff.clear();
		diff.resize(n.size());
		for(unsigned j=0; j<n.size(); j++){
			diff[j] = (n[j]-n0[j]);
		}
		sum1.resize(rr.size());
		for(unsigned j=0; j<rr.size(); j++){
			s = 0;
			for(unsigned k=0; k<A12[rr[j]].size(); k++){
				s += A12[rr[j]][k]*diff[k];
			}
			sum1[j] = s;
		}
		s = 0;
		for(unsigned j=0; j<sum1.size(); j++){
			s += sum1[j];
		}
		sum1.clear();
		res[indexingI[i]] -= section * s / deltat;

		// Electron current.
		s=0;
		for(unsigned j=0; j<rr.size(); j++){
			s += resn[rr[j]];
		}
		res[indexingI[i]] += section * q * s;
	}	
	rr.clear();
	sum1.clear();
	diff.clear();
	resn.clear();
	
	
	for(int i=0; i<nnodes; i++){
		res[ordV(i)] /= rowscaling[0];
	}
	for(int i=0; i<nnodes; i++){
		res[ordn(i)] /= rowscaling[1];
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
Newton::org_secs2d_newton_jacobian(	Probl& P, std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, double deltat, BCS_CIRC& bcs, 
									ordering& ordV, ordering& ordn, std::vector<int>& indexingF, std::vector<int>& indexingI, sparse_matrix& jacobian)
{
	int	nnodes = P.get_msh_nodes(),
		nelements = P.get_msh_elem(),
        numcontacts, ndofs, j = 0;
		
	double	eps_semic = P._eps_semic,
			eps_ins = P._eps_ins,
			q = P._q,
			Vth = P._Vth,
			section = P._section;
			
	bool	ins = P._ins;
			
	std::vector<int>::iterator it;
	
	std::vector< std::vector<int> >	dnodes = P._dnodes;
	
	std::vector<double>	epsilon(nelements,eps_semic),
						rowscaling = P._rowscaling,
						colscaling = P._colscaling;
						
	std::array<int,2>	pins = P._pins;
	
	std::vector<int>	scnodes = P._scnodes,
						insulator = P._insulator;
						
	sparse_matrix 	A, B, r;
	
	//numscnodes = std::accumulate( scnodes.begin(), scnodes.end(), 0.0);
	
	numcontacts = pins.size();
	ndofs =  2 * nnodes + F.size() + numcontacts;	

	jacobian.resize(ndofs);
	
	///	COMPUTING COEFFICIENTS	
	if(ins){
		for(unsigned i=0; i<insulator.size(); i++){	
			if(insulator[i] == 1){
				epsilon[i] = eps_ins;
			}
		}
	}
	
	/// COMPUTING FIRST ROW
	sparse_matrix	M;

	std::vector<double> eones(insulator.size(),1.0),
						not_ins(insulator.size(),0.0),
						ones(nnodes,1.0),
						psi(nnodes,0.0);
	
	M.resize(nnodes);
	
	for(unsigned i=0; i<insulator.size(); i++){
		if(insulator[i] == 0){
			not_ins[i] = 1;
		}
	}

	bim2a_reaction (P._msh, not_ins, ones, jacobian, ordV, ordn);
	for(int i=0; i<nnodes; i++){		
		jacobian[ordV(i)][ordn(i)] *= q;
	}
	
	bim2a_reaction (P._msh, eones, ones, M);

	bim2a_advection_diffusion (P._msh, epsilon, psi, jacobian, ordV, ordV);
	
	///	Physical models.
	std::vector<double>	mobn(nelements,0.0),
						alphan(n.size(),0.0),
						der_dalpha_n(nelements,0.0);

	org_physical_models2d(n, P, mobn, alphan, der_dalpha_n);
	
	///	COMPUTING SECOND ROW
	sparse_matrix::col_iterator J;
					
	std::vector<double>	fluxn(n.size(),0.0),
						alfa(insulator.size(),0.0),
						mobn_n(insulator.size(),0.0);	

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
	
	bim2a_advection_diffusion (P._msh, mobn_n, psi, jacobian, ordn, ordV);
		
	if(alphan.size() == n.size() && alphan.size() == V.size()){
		fluxn = alphan;
		for(unsigned i=0; i<scnodes.size(); i++){
			if(scnodes[i] == 1){
				fluxn[i] *= (-1);
				fluxn[i] += V[i]/Vth;
			}
		}
	
	}
	else{
		std::cout<<"error: org_secs2d_newton_jacobian, dimensions mismatch"<<std::endl;
	}
	
	bim2a_advection_diffusion( P._msh, alfa, fluxn, jacobian, ordn, ordn);

	for(unsigned i=0; i<ones.size(); i++){
		ones[i] = ones[i]/deltat;
	}
	
	bim2a_reaction (P._msh, not_ins, ones, jacobian, ordn, ordn);

	///	ADJUST FOR ZERO INSULATOR CHARGE
	
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 0){
			for (J = jacobian[ordn(i)].begin (); J != jacobian[ordn(i)].end (); ++J){
				jacobian[ ordn(i) ][ jacobian.col_idx (J) ] = 0;
			}
		}
	}
	
	for(unsigned i=0; i<scnodes.size(); i++){
		if(scnodes[i] == 0){
			jacobian[ ordn(i) ][ ordn(i) ] += M[i][i];
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
			if( (unsigned)B.col_idx(J) == j ){
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
	std::vector<double>	s1(nnodes,0.0),
						s2(nnodes,0.0),
						zeros(nnodes,0.0);
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
			for(int j=0; j<nnodes; j++){
				s1[j] += jacobian[ordV(rr[k])][ordV(j)];
				s2[j] += jacobian[ordV(rr[k])][ordn(j)];
			}
		}
	
		// Displacement current.
		for(int j=0; j<nnodes; j++){
			jacobian[indexingI[i]][ordV(j)] -= section * s1[j] / deltat;
		}
		for(int j=0; j<nnodes; j++){
			jacobian[indexingI[i]][ordn(j)] -= section * s2[j] / deltat;
		}
		
		// Electron current.
		s1 = zeros;
		s2 = zeros;
		for(unsigned k=0; k<rr.size(); k++){
			for(int j=0; j<nnodes; j++){
				s1[j] += jacobian[ordn(rr[k])][ordV(j)];
				s2[j] += jacobian[ordn(rr[k])][ordn(j)];
			}
		}
		
		for(int j=0; j<nnodes; j++){
			jacobian[indexingI[i]][ordV(j)] -= (-section * q * s1[j]);
		}
		for(int j=0; j<nnodes; j++){
			jacobian[indexingI[i]][ordn(j)] -= (-section * q * s2[j]);
		}
		rr.clear();
	}
	zeros.clear();
	s1.clear();
	s2.clear();

	
	for(int i=0; i<nnodes; i++){
		for (J = jacobian[ordV(i)].begin (); J != jacobian[ordV(i)].end (); ++J){
			jacobian[ordV(i)][jacobian.col_idx (J)] /= rowscaling[0];					
		}
	}
	
	for(int i=0; i<nnodes; i++){
		for (J = jacobian[ordn(i)].begin (); J != jacobian[ordn(i)].end (); ++J){
			jacobian[ordn(i)][jacobian.col_idx (J)] /= rowscaling[1];					
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
		while( jacobian.col_idx(J)%2 != 0){
			J++;
			if( J == jacobian[i].end () ){
				break;
			}
		}
		if(J != jacobian[i].end ()){
		for(int j=0; j<nnodes; j++){
			if( (unsigned)jacobian.col_idx(J) == ordV(j) ){
				jacobian[i][ordV(j)] *= colscaling[0];
				J++;
				while( jacobian.col_idx(J)%2 != 0){
					J++;
					if( J == jacobian[i].end () ){
						j = nnodes;
						break;
					}
				}
			}
		}
		}
	}
	
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		while( jacobian.col_idx(J)%2 == 0 ){
			J++;
			if( J == jacobian[i].end () ){
				break;
			}
		}
		if(J != jacobian[i].end ()){
		for(int j=0; j<nnodes; j++){
			if( (unsigned)jacobian.col_idx(J) == ordn(j) ){
				jacobian[i][ordn(j)] *= colscaling[1];
				J++;
				while( jacobian.col_idx(J)%2 == 0 ){
					J++;
					if( J == jacobian[i].end () ){
						j = nnodes;
						break;
					}
				}
			}
		}
		}
	}

	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		while(jacobian.col_idx(J) < indexingF[0]){
			J++;
			if(J == jacobian[i].end ()){
				break;
			}
		}
		if(J != jacobian[i].end ()){
		for(unsigned j=0; j<indexingF.size(); j++){
			if( jacobian.col_idx(J) == indexingF[j] ){
				jacobian[i][indexingF[j]] *= colscaling[2];
				J++;
				if( J == jacobian[i].end () || (unsigned)jacobian.col_idx(J) > indexingF.size() ){
					break;
				}
			}
		}
		}
	}
	for(unsigned i=0; i<jacobian.rows(); i++){
		J = jacobian[i].begin();
		while(jacobian.col_idx(J) < indexingI[0]){
			J++;
			if(J == jacobian[i].end ()){
				break;
			}
		}
		if(J != jacobian[i].end ()){
		for(unsigned j=0; j<indexingI.size(); j++){
			if( jacobian.col_idx(J) == indexingI[j] ){
				jacobian[i][indexingI[j]] *= colscaling[3];
				J++;
				if( J == jacobian[i].end () ){
					break;
				}				
			}
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

void
Newton::org_physical_models2d (	std::vector<double>& n, Probl& P,
								// output
								std::vector<double>& mobilityn, std::vector<double>& alphan)
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

	n_ref = org_gaussian_charge_n( phi_lumo_ref, P);
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

void
Newton::org_gaussian_charge_n(	std::vector<double>& V, Probl& P,
								// output
								std::vector<double>& rhon)
{
	double	q = P._q;
    std::vector<double> n = n_approx(V,P);
	
	rhon = n;
    for(unsigned i=0; i<n.size(); i++){
		rhon[i] *= -q;
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
Newton::org_secs_state_predict (	Probl& P, std::vector<double>& Vold, std::vector<double>& nold, std::vector<double>& Fold, std::vector<double>& Iold,
									std::vector<double>& Voldold, std::vector<double>& noldold, std::vector<double>& Foldold, std::vector<double>& Ioldold,
									int& tstep, std::vector<double>& tout,
									// output
									std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0)
{
  if (tstep > 2){

	int	length, j;
    double  dt, difft;
    std::vector<double> it(2,0), t(2,0), dndt, dlndt, dVdt, dIdt, dFdt, lnold, lnoldold;

    it[0] = (tstep - 2);
    it[1] = (tstep - 1);

    t[0] = (tout[it[0]] - tout[tstep - 2]);
    t[1] = (tout[it[1]] - tout[tstep - 2]);
	
    dt = tout[tstep] - tout[tstep - 1];

	dndt.resize(nold.size());
	for(unsigned i=0; i<nold.size(); i++){
        dndt[i] = (nold[i] - noldold[i]);
    }
	
    difft = t[1]-t[0];
    for(unsigned i=0; i<dndt.size(); i++){
		dndt[i] /= difft;
    }

    // n
    //ln = log (n(device.scnodes, it(1:2)));
	length = std::accumulate( P._scnodes.begin(), P._scnodes.end(), 0.0);
	lnold.resize(length);
	lnoldold.resize(length);
    j = 0;
	for(unsigned i=0; i<P._scnodes.size(); i++){
			if(P._scnodes[i] == 1){
				lnold[j] = log(nold[i]);
				lnoldold[j] = log(noldold[i]);
				j++;
			}
    }

    //dlndt  = diff (ln, 1, 2) ./ diff (t);
	dlndt.resize(lnold.size());
	for(unsigned i=0; i<lnold.size(); i++){
        dlndt[i] = (lnold[i] - lnoldold[i]);
    }
	
	for(unsigned i=0; i<dlndt.size(); i++){
		dlndt[i] /= difft;
    }
	
	//n0 = n[tstep-1];
	n0 = nold;

    // V, F, I
    //dVdt  = diff (V(:, it(1:2)), 1, 2) ./ diff (t);

	dVdt.resize(Vold.size());
	for(unsigned i=0; i<Vold.size(); i++){
        dVdt[i] = (Vold[i] - Voldold[i]);
    }

	for(unsigned i=0; i<dVdt.size(); i++){
		dVdt[i] /= difft; 
    }

    //dFdt  = diff (F(:, it(1:2)), 1, 2) ./ diff (t);
	dFdt.resize(Fold.size());
	for(unsigned i=0; i<Fold.size(); i++){
        dFdt[i] = (Fold[i] - Foldold[i]);
    }
	
	for(unsigned i=0; i<dFdt.size(); i++){
		dFdt[i] /= difft;
    }

    //dIdt  = diff (I(:, it(1:2)), 1, 2) ./ diff (t);
	dIdt.resize(Iold.size());
	for(unsigned i=0; i<Iold.size(); i++){
        dIdt[i] = (Iold[i] - Ioldold[i]);
    }
	
	for(unsigned i=0; i<dIdt.size(); i++){
		dIdt[i] /= difft;
    }

    //V0 = V(:, it(2)) + dVdt * dt;
	V0 = dVdt;
	for(unsigned i=0; i<V0.size(); i++){
		V0[i] *= dt ;
	}
	if(V0.size() != Vold.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<V0.size(); i++){
			V0[i] += Vold[i];
		}
	}

    //n0(device.scnodes) = exp (ln(:, 2) + dlndt * dt);
	for(unsigned i=0; i<dlndt.size(); i++){
		dlndt[i] *= dt ;
	}
	
	if(dlndt.size() != lnold.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<dlndt.size(); i++){
			dlndt[i] += lnold[i];
		}	
		for(unsigned i=0; i<P._scnodes.size(); i++){
			if(P._scnodes[i] == 1){
				n0[i] = std::exp(dlndt[i]);
			}
		}
	}

    //F0 = F(:, it(2)) + dFdt * dt;
	F0 = dFdt;
	for(unsigned i=0; i<F0.size(); i++){
		F0[i] *= dt ;	
	}
	if(F0.size() != Fold.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<F0.size(); i++){
			F0[i] += Fold[i];
		}
	}

    //I0 = I(:, it(2)) + dIdt * dt;
	I0 = dIdt;
	for(unsigned i=0; i<I0.size(); i++){
		I0[i] *= dt ;	
	}
	if(I0.size() != Iold.size()){
		std::cout<<"error: dimensions mismatch, org_state_predict"<<std::endl;
	}
	else{
		for(unsigned i=0; i<I0.size(); i++){
			I0[i] += Iold[i];
		}
	}
  }
  else{
	// passo ai vettori _0 i valori al passo temporale precedente (t_step-1)
	V0 = Vold;
	n0 = nold;
	F0 = Fold;
	I0 = Iold;
  }
};

void 
Newton::compute_residual_norm (	double& resnrm, int& whichone, std::vector<double>& resall, std::vector<double>& res,
								int nnodes, ordering& idxV, ordering& idxn, std::vector<int>& idxF, std::vector<int>& idxI)
{
	std::vector<double>	aux(nnodes,0);
	
	for(int i=0; i<nnodes; i++){
		aux[i] = std::abs( res[ idxV(i) ] );
	}
	resall[0] = *std::max_element( aux.begin(),aux.end() );
  	
	aux.resize(nnodes);
	for(int i=0; i<nnodes; i++){
		aux[i] = std::abs( res[idxn(i)] );
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
  
	incrV = infnorm(diff) / (infnorm(Va) * clamping + cl[0]);

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
  
	incrn = Vth * infnorm(q1) / (infnorm(diff) * clamping + cl[1]);	

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
	incrF = infnorm(diff) / (infnorm(Fa) * clamping + cl[2]);
	
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
	incrI = infnorm(diff) / (infnorm(Ia) * clamping + cl[3]);
	
	diff.clear();
};

void 
Newton::org_secs_safe_increment (	std::vector<double>& V0, std::vector<double>& n0, std::vector<double>& F0, std::vector<double>& I0,
									std::vector<double>& dV, std::vector<double>& dn, std::vector<double>& dF, std::vector<double>& dI,
									Probl& P,
									std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I,
									double& clamp, double& tauk	)
{
	double	tk = 1.0,
			min;
	std::vector<double>	n0_2;
	std::vector<int>	scnodes = P._scnodes,
						clamping = P._clamping,
						where(n0.size(),0);

	if (n0.size() != dn.size()){
		std::cout<<"error: org_secs_safe_increment, vectors with different size"<<std::endl;
	}
	for(unsigned i=0; i<n0.size(); i++){
		if(scnodes[i] == 1){
			if(n0[i]+dn[i] <= 0)
				where[i] = 1;
		}
	}
	
	int s = std::accumulate( scnodes.begin(), scnodes.end(), 0.0);
	where.resize(s);

	if ( any(where) ){
		for(unsigned i=0; i<where.size(); i++){
				if(where[i]==1){
					min = n0[i]/std::abs(dn[i]);
					if(tk > min){
						tk = min;
					}
				}
		}
		tk *= 0.9;
	}

  clamp = 1;

  if (P._clampOnOff){
    // V
    if (any(dV)){
		clamp = std::fmin (1, clamping[0] / infnorm(dV));
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
      clamp = std::fmin (clamp, (exp (clamping[1] / P._Vth) - 1) / infnorm (dn2) );
    }

    if (any (dn3)){
      for(unsigned i=0; i<dn.size(); i++){
        dn3[i] = dn3[i]*(dn[i]/n0[i]);
      }
      clamp = std::fmin (clamp, (1 - exp (-clamping[1]/P._Vth)) / infnorm (dn3) );
    }

    // F
    if (any(dF)){
      clamp = std::fmin (clamp, clamping[2] / infnorm (dF));
    }

    // I
    if (any(dI)){
      clamp = std::fmin (clamp, clamping[3] / infnorm (dI));
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
  else{
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
	}
  }	
};

void
Newton::CONV_MSG (int tstep, int Nstep, int mNstep, double t, std::string reason, int field, double incr, std::vector<double>& res)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep-1<<", modified Newton iteration "<<mNstep-1<<", model time "<<t   //
            <<": convergence reached ("<<reason<<")"<<std::endl;
  std::cout<<"incr ("<<field+1<<") = "<<incr<<" residual = [ ";

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
  std::cout<<" the increment in field "<<(field+1)<<" has grown too large: "<<std::endl;
  std::cout<<" "<<incr<<" > "<<P._maxnpincr<<std::endl;
  std::cout<<"incr0 ( "<<field+1<<" ) = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }

  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};

void
Newton::DIV_MSG (	int tstep, double t, int Nstep, int field, std::vector<double>& incrhist, //
					double incr, std::vector<double>& res, int nsteps_check)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<" model time "<<t<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep<<", the Newton algorithm is diverging: "<<std::endl;
  std::cout<<"the increment in field "<<field+1<<" is not decreasing : "<<std::endl;
  std::cout<<incrhist[Nstep-1]<<" > "<<incrhist[Nstep - nsteps_check-1]<<std::endl;
  std::cout<<"incr ("<<field+1<<") = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }
  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};

void
Newton::DIV_MN_MSG (	int tstep, double t, int Nstep, int mNstep, int field, std::vector<double>& incrhist,	//
						double incr, std::vector<double>& res, int nsteps_check)
{
  std::cout<<" "<<std::endl;
  std::cout<<"at time step "<<tstep<<"model time "<<t<<std::endl;
  std::cout<<"fixed point iteration "<<Nstep<<", modified Newton iteration "<<mNstep<<std::endl;
  std::cout<<"the Modified Newton algorithm is diverging: "<<std::endl;
  std::cout<<"the increment in field "<<field+1<<" is not decreasing : "<<std::endl;
  std::cout<<incrhist[mNstep-1]<<" > "<<incrhist[mNstep - nsteps_check-1]<<std::endl;
  std::cout<<"incr ("<<field+1<<") = "<<incr<<" residual = [ ";

  for(unsigned i=0; i<res.size()-1; i++){
    std::cout<<res[i]<<" ";
  }
  std::cout<<res.back()<<" ]"<<std::endl;
  std::cout<<" "<<std::endl;
};


void
Newton::saveNEWT (	std::vector<double>& Vold, std::vector<double>& nold, std::vector<double>& Fold, std::vector<double>& Iold, double told, 
					std::vector<double>& V, std::vector<double>& n, std::vector<double>& F, std::vector<double>& I, std::vector<double>& res,
					double t, double dt, int nsaves, int newton_solves, int modified_newton_solves, double freq)
{

  ColumnVector oct_V (V.size (), 0.0);
  ColumnVector oct_n (n.size (), 0.0);
  ColumnVector oct_F (F.size (), 0.0);
  ColumnVector oct_I (I.size (), 0.0);
  
  ColumnVector oct_Vold (Vold.size (), 0.0);
  ColumnVector oct_nold (nold.size (), 0.0);
  ColumnVector oct_Fold (Fold.size (), 0.0);
  ColumnVector oct_Iold (Iold.size (), 0.0);
  
  ColumnVector oct_res (res.size (), 0.0);

  std::copy_n (V.begin (), V.size (), oct_V.fortran_vec ());
  std::copy_n (n.begin (), n.size (), oct_n.fortran_vec ());
  std::copy_n (F.begin (), F.size (), oct_F.fortran_vec ());
  std::copy_n (I.begin (), I.size (), oct_I.fortran_vec ());
  
  std::copy_n (Vold.begin (), Vold.size (), oct_Vold.fortran_vec ());
  std::copy_n (nold.begin (), nold.size (), oct_nold.fortran_vec ());
  std::copy_n (Fold.begin (), Fold.size (), oct_Fold.fortran_vec ());
  std::copy_n (Iold.begin (), Iold.size (), oct_Iold.fortran_vec ());
  
  std::copy_n (res.begin (), res.size (), oct_res.fortran_vec ());
  
  octave_io_mode m = gz_write_mode;
  
  // Define filename.
  char FileName[255] = "";
  sprintf(FileName,"NEWT_freq_%f_output_%d.gz",freq,nsaves);
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);

  assert (octave_save ("told", told) == 0);
  assert (octave_save ("Vold", oct_Vold) == 0);
  assert (octave_save ("nold", oct_nold) == 0);
  assert (octave_save ("Fold", oct_Fold) == 0);
  assert (octave_save ("Iold", oct_Iold) == 0);

  assert (octave_save ("t", t) == 0);
  assert (octave_save ("V", oct_V) == 0);
  assert (octave_save ("n", oct_n) == 0);
  assert (octave_save ("F", oct_F) == 0);
  assert (octave_save ("I", oct_I) == 0);
  
  assert (octave_save ("res", oct_res) == 0);
  assert (octave_save ("dt", dt) == 0);
  assert (octave_save ("newton_solves", newton_solves) == 0);
  assert (octave_save ("modified_newton_solves", modified_newton_solves) == 0);

  assert (octave_io_close () == 0);

};

double
Newton::infnorm(std::vector<double>& in){
    double out = 0;
	std::vector<double>	v(in.size(),0);

	for(unsigned i=0; i<in.size(); i++){
		v[i] = std::abs(in[i]);
	}
	out = *std::max_element( v.begin(),v.end() );

	v.clear();
    return out;
};

bool
Newton::any(std::vector<int>& v)
{
	// int l=0;
	// for(unsigned i=0; i<v.size(); i++){
		// if(v[i]!=0)
			// l++;
	// }
	int l = std::accumulate( v.begin(), v.end(), 0.0);
	
	// if( (unsigned)l == v.size())
	if( l > 0 )
		return true;
	else
		return false;
};

bool
Newton::any(std::vector<double>& v)
{
	// int l=0;
	// for(unsigned i=0; i<v.size(); i++){
		// if(v[i]!=0)
			// l++;
	// }
	int l = std::accumulate( v.begin(), v.end(), 0.0);
	
	// if( (unsigned)l == v.size())
	if( l > 0 )
		return true;
	else
		return false;
};


void
Newton::saveJAC (int nrows, int ncols, std::vector<double>& vals)
{
  Matrix oct_jac (nrows, ncols, 0.0);
  
  std::copy_n (vals.begin (), vals.size (), oct_jac.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("jac", oct_jac);
  
  octave_io_mode m = gz_write_mode;
  
  // Define filename.
  char FileName[255] = "";
  sprintf(FileName,"NEWT_jac.gz");
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);
  assert (octave_save ("jac", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};


void
Newton::saveVn(std::vector<double>& V, std::vector<double>& n, const char* FileName)
{
  ColumnVector oct_V (V.size (), 0.0);
  ColumnVector oct_n (n.size (), 0.0);

  std::copy_n (V.begin (), V.size (), oct_V.fortran_vec ());
  std::copy_n (n.begin (), n.size (), oct_n.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("V", oct_V);
  the_map.assign ("n", oct_n);
  
  octave_io_mode m = gz_write_mode;
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);
  assert (octave_save ("Newt", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};

void
Newton::saveRES(std::vector<double>& res, const char* FileName)
{
  ColumnVector oct_res (res.size (), 0.0);

  std::copy_n (res.begin (), res.size (), oct_res.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("res", oct_res);
  
  octave_io_mode m = gz_write_mode;
  
  // Save to filename.
  assert (octave_io_open (FileName, m, &m) == 0);
  assert (octave_save ("Res", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};