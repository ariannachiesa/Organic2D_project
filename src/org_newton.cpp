/*! \file org_newton.cpp
  \brief Implementation of class Newton's methods
*/

#include "newton.h"

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