/*! \file newton12.cpp
  \brief Class Newton: 1st + 2nd eq. only
*/

#include "newton12.h"

Newton::Newton(	Probl& P, std::vector<double>& Vin, std::vector<double>& nin, 
				std::vector<double>& tspan, double freq)
{
	int	nnodes = P.get_msh_nodes(),
		newton_solves = 0,
        modified_newton_solves = 0,
		rejected = 0,
		j = 0,
		lastsaved = 0,
		nsaves = 0,
		nsteps_check = P._nsteps_check,
		maxit = P._maxit,
		firstfixtstep, tstep, in,
		whichone, totmn, iimn, told;

	double	t, dt, y0, yend, dtfact, rejmn,
			dt0 = P._dt0,
			dtmax = P._dtmax,
			dtmin = P._dtmin,
			clamp = 0,
			tauk = 0,
			incr0, incrlast, incr1, incrk,
			incr0v = 0, incr0n = 0,
			incr1v = 0, incr1n = 0,
			incrkv = 0, incrkn = 0;
			
	bool	reject, rejectmnewton, convergedmnewton;
	
	auto quad = P._msh.begin_quadrant_sweep();
	
	std::vector<double>
						delta, resnrm, resnrmk,
						inc_clamp, inck_clamp, incnrm, incnrmk,
						V0, n0,
						V1, n1,
						V2, n2,
						Vk, nk,
						Vold, nold,
						Voldold, noldold,
						dV, dn,
						resall(2,0),
						vecincr(2,0);
						
	P._rowscaling.resize(2);
	P._colscaling.resize(2);
			
	dt = dt0;
	dt = std::fmin(dt, dtmax);
	
	/// INIT STATE FROM SCRATCH
	tstep = 0;
	_tout.resize(3);
	_tout[0] = tspan[0];
	_tout[1] = tspan[0];
	_tout[2] = tspan[0];
    t = tspan[0];
	told = t;
	
	_V = Vin;
	Vold = Vin;
	Voldold = Vin;
	_n = nin;
	nold = nin;
	noldold = nin;

    firstfixtstep = 2;
	
	y0 = quad->p(1,0);
	std::cout<<"Connecting contact at y = "<<y0<<" to circuit pin "<<P._pins[0]<<std::endl;

	for (auto quadrant = P._msh.begin_quadrant_sweep (); quadrant != P._msh.end_quadrant_sweep (); ++quadrant){
		yend = quadrant->p(1,3);
	}
	std::cout<<"Connecting contact at y = "<<yend<<" to circuit pin "<<P._pins[1]<<std::endl;
	
	/// Node ordering
	ordering	ordV = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 0> (gt); },
				ordn = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 1> (gt); };
	
	//for (unsigned n_fix_tstep = firstfixtstep; n_fix_tstep<tspan.size()+1; n_fix_tstep++){   /// TIME FIXED SPAN
unsigned n_fix_tstep = firstfixtstep;
		t = tspan[n_fix_tstep - 2];
		
		while (t < tspan[n_fix_tstep-1]){ /// TIME STEP
		//while (tstep < 14){ /// TIME STEP
			
			std::cout<<"--"<<std::endl;
		
			tstep++;

			_tout[2] = std::fmin (t + dt, tspan[n_fix_tstep - 1]);

			t = _tout[2];

			dt = t - _tout[1];
			
			incr0 = 4 * P._maxnpincr;
			
			org_secs_state_predict (P, Vold, nold, Voldold, noldold, tstep, _tout, V0, n0);
			
			V2 = V0;
			n2 = n0;
			
			in = 0;
			rejmn = 0;	///
			totmn = 0;	///

			reject = false;
			
			incnrm.resize(maxit);
			inc_clamp.resize(maxit);
			resnrm.resize(maxit);
			
			std::cout << "new Newton step = "<<std::endl;
			std::cout << "tstep = "<<tstep<<std::endl;
			
			while (!reject && (in < P._maxit)){ /// NEWTON STEP
				in += 1;
				std::cout << "in = "<<in<<std::endl;

				V1 = V2;
				n1 = n2;					
				
				_res = org_secs2d_newton_residual(P, V2, n2, Vold, nold, dt, ordV, ordn);
				
				org_secs2d_newton_jacobian(	P, V2, n2, dt, ordV, ordn, _jac);
				
				/// Dirichlet BCs on V:
				int indexT = P._nTrees-1;
				double	PhiB = P._PhiB,
						Vshift = P._Vshift,
						Vgate = P._VG;
				std::tuple<int, int, func_quad>	tupla1(0,2,[&PhiB,&V2](tmesh::quadrant_iterator quad, tmesh::idx_t i)
														{return (PhiB-V2[quad->gt(i)]);}),				// impongo PhiB		
												tupla2(indexT,3,[&Vgate,&Vshift,&V2](tmesh::quadrant_iterator quad, tmesh::idx_t i)
														{return (Vgate+Vshift-V2[quad->gt(i)]);});		// impongo Vshift + Vgate
																			
				dirichlet_bcs_quad	bcsV;
				bcsV.push_back(tupla1);
				bcsV.push_back(tupla2);

				bim2a_dirichlet_bc (P._msh, bcsV, _jac, _res, ordV);
				
				/// Dirichlet BCs on n:
				double rho, nimposed;

				rho = org_gaussian_charge_n( V2[0], P);
				nimposed = rho * (-1)/P._q ;
	
				std::tuple<int, int, func_quad>	tuplan(0,2,[&nimposed,&n2](tmesh::quadrant_iterator quad, tmesh::idx_t i)
																{return (nimposed-n2[quad->gt(i)]);});
				dirichlet_bcs_quad	bcsn;
				bcsn.push_back(tuplan);
	
				bim2a_dirichlet_bc (P._msh, bcsn, _jac, _res, ordn);
				
				// for(int i=0; i<nnodes; i++){
					// std::cout<<"resV = "<<_res[ordn(i)]<<std::endl;
				// }				
				
				if (in == 1){
					whichone = 0;
					compute_residual_norm (resnrm[0],whichone,resall,_res,nnodes,ordV,ordn);
					
					// for(unsigned i=0; i<resall.size(); i++){
						// std::cout<<"resall = "<<resall[i]<<std::endl;
					// }
					// std::cout<<"resnrm = "<<resnrm[in-1]<<std::endl;
					
					if(P._rowscaling.size() == resall.size()){
						for(unsigned i=0; i<P._rowscaling.size(); i++){
							P._rowscaling[i] = P._rowscaling[i] * (resall[i]+1);
							// std::cout<<"P._rowscaling[i] = "<<P._rowscaling[i]<<std::endl;
						}
					}
					else{
						std::cout<<"error: dimensions mismatch"<<std::endl;
						break;
					}
					
					_res = org_secs2d_newton_residual(P, V2, n2, Vold, nold, dt, ordV, ordn);
					
					// re-impongo le BC sul nuovo vettore residuo
					bim2a_dirichlet_bc (P._msh, bcsV, _jac, _res, ordV,true);
					bim2a_dirichlet_bc (P._msh, bcsn, _jac, _res, ordn,true);
				}				
				
				// for(int i=0; i<nnodes; i++){
					// std::cout<<"resV = "<<_res[ordV(i)]<<std::endl;
				// }				
				
				resall.resize(2);
				compute_residual_norm (resnrm[in-1],whichone,resall,_res,nnodes,ordV,ordn);
					
				// for(unsigned i=0; i<resall.size(); i++){
					// std::cout<<"resall = "<<resall[i]<<std::endl;
				// }
				// std::cout<<"resnrm = "<<resnrm[in-1]<<std::endl;

				
				/// Solve non.linear system.
				std::cout << "Solving linear system."<<std::endl;
		
				delta.resize(_res.size());
				delta = _res;
		
				mumps mumps_solver;
      
				std::vector<double> vals(nnodes,0.0);
				std::vector<int> 	irow(nnodes,0),
									jcol(nnodes,0);
							
				mumps_solver.init();
	  
				_jac.aij(vals, irow, jcol, mumps_solver.get_index_base ());
		
				if(in == 1){	// only at the first iteration of the Newton method
					mumps_solver.set_lhs_structure (_jac.rows(), irow, jcol);
					mumps_solver.analyze ();
				}
		
				mumps_solver.set_lhs_data (vals);
				mumps_solver.set_rhs (delta);
				mumps_solver.factorize ();
		
				mumps_solver.solve ();
		
				for(unsigned i=0; i<delta.size(); i++){
					delta[i] *= (-1);
				}
				
				newton_solves +=1;
	
				dV.resize(nnodes);
				for(int i=0; i<nnodes; i++){
					dV[i] = delta[ordV(i)] * P._colscaling[0];
				}				

				dn.resize(nnodes);
				for(int i=0; i<nnodes; i++){
					dn[i] = delta[ordn(i)] * P._colscaling[1];
				}
				delta.clear();
				
				V2.clear(); n2.clear();
				org_secs_safe_increment (V1, n1, dV, dn, P, V2, n2, clamp, tauk);

				// std::cout<<"clamp = "<<clamp<<std::endl;
				// std::cout<<"tauk = "<<tauk<<std::endl;
												
				if ((clamp <= 0) || (tauk <= 0)){
					reject = true;
					break;
				}
		
				compute_variation ( V0, n0, V2, n2, P, std::fmin (clamp, tauk), incr0v, incr0n);
									
				vecincr[0] = incr0v;
				vecincr[1] = incr0n;

				incr0 = *std::max_element(vecincr.begin(),vecincr.end());
				if(incr0==incr0v) whichone = 0;
				if(incr0==incr0n) whichone = 1;
				
				// for(unsigned i=0; i<vecincr.size(); i++){
					// std::cout<<"vecincr0 = "<<vecincr[i]<<std::endl;
				// }
				
				if (incr0 > P._maxnpincr && dt > P._dtmin){
					MAXINCR_MSG (tstep, t, in, whichone, incr0, resall, P);
					reject = true;
					break;
				}

				compute_variation (	V1, n1, V2, n2, P, std::fmin (clamp, tauk), incr1v, incr1n);
										
				vecincr[0] = incr1v;
				vecincr[1] = incr1n;
				
				// for(unsigned i=0; i<vecincr.size(); i++){
					// std::cout<<"vecincr1 = "<<vecincr[i]<<std::endl;
				// }

				incr1 = *std::max_element(vecincr.begin(),vecincr.end());
				if(incr1==incr1v) whichone = 0;
				if(incr1==incr1n) whichone = 1;
				
				incnrm[in-1] = incr1;
				inc_clamp[in-1] = incr1 / clamp;
				incrlast = incr1;
	
				if (resnrm[in-1] < P._toll ){
					CONV_MSG (tstep, in, 1, t, "res", whichone, incr1, resall);
					break;
				}

				if (in > P._nsteps_check && incnrm[in-1] > incnrm[in - P._nsteps_check -1]	//
					&& resnrm[in-1] > resnrm[in - P._nsteps_check -1] && dt > dtmin){
					DIV_MSG (tstep, t, in, whichone, incnrm, incr1, resall, nsteps_check);
					reject = true;
					break;
				}

				std::cout<<"incr ("<<whichone<<") = "<<incr1;
				std::cout<<"   residual = [ "<<resall[0]<<" "<<resall[1]<<" ] "<<std::endl;
				
				if (incr1 < P._toll ){
					CONV_MSG(tstep, in, 1, t, "incr", whichone, incr1, resall);
					break;
				}
				
				/// MODIFIED NEWTON
				V1 = V2;
				n1 = n2;

				rejectmnewton = false;
				convergedmnewton = false;
				
				incnrmk.resize(P._maxit_mnewton+1);
				inck_clamp.resize(P._maxit_mnewton+1);
				resnrmk.resize(P._maxit_mnewton+1);
				
				for (int imn = 1; imn<P._maxit_mnewton+1; imn++) {	/// MODIFIED NEWTON STEP 

					iimn = imn;
		
					Vk = V2;
					nk = n2;

					// for(unsigned i=0; i<nnodes; i++){
						// std::cout<<"V2 = "<<V2[i]<<std::endl;
					// }					
					
					_res = org_secs2d_newton_residual(P, V2, n2, Vold, nold, dt, ordV, ordn);
					
					// for(unsigned i=0; i<nnodes; i++){
						// std::cout<<"resV = "<<_res[ordV(i)]<<std::endl;
					// }				
					
					/// Dirichlet BCs on V:
					bim2a_dirichlet_bc (P._msh, bcsV, _jac, _res, ordV, true);
					
					/// Dirichlet BCs on n:
					rho = org_gaussian_charge_n(V2[0], P);
					nimposed = rho * (-1)/P._q;
	
					bim2a_dirichlet_bc (P._msh, bcsn, _jac, _res, ordn, true);
														
					resall.resize(2);
					compute_residual_norm (resnrmk[imn-1], whichone, resall, _res, nnodes, ordV, ordn);
					
					// for(unsigned i=0; i<resall.size(); i++){
						// std::cout<<"resall = "<<resall[i]<<std::endl;
					// }
					// std::cout<<"resnrm = "<<resnrmk[imn-1]<<std::endl;					
				
					/// Solve non linear system.
					std::cout << "Solving linear system - Modified Newton."<<std::endl;
		
					delta = _res;
					
					mumps_solver.set_rhs (delta);
					
					mumps_solver.solve ();
	  
					modified_newton_solves +=1;
					
					for(int i=0; i<nnodes; i++){
						dV[i] = delta[ordV(i)] * P._colscaling[0];
					}
					for(int i=0; i<nnodes; i++){
						dn[i] = delta[ordn(i)] * P._colscaling[1];
					}

					V2.clear(); n2.clear();  
					org_secs_safe_increment (Vk, nk, dV, dn, P, V2, n2, clamp, tauk);
					
					if ((tauk <=0) || (clamp <= 0)){
						reject = true;
						break;
					}
		  
					compute_variation(Vk, nk, V2, n2, P, std::fmin (clamp, tauk), incrkv, incrkn);
		
					vecincr[0] = incrkv;
					vecincr[1] = incrkn;
		
					incrk = *std::max_element(vecincr.begin(),vecincr.end());

					if(incrk==incrkv) whichone = 0;
					if(incrk==incrkn) whichone = 1;
					
					incnrmk[imn - 1] = incrk;
					inck_clamp[imn - 1] = incrk / clamp;
					incrlast = incrk;

					if (resnrmk[imn - 1] < P._toll){
						CONV_MSG (tstep, in, imn, t, "res", whichone, incrk, resall);
						rejectmnewton = false;
						convergedmnewton = true;
						break;
					}

					if (imn > P._nsteps_check && (resnrmk[imn - 1] > resnrmk[imn - P._nsteps_check - 1])
						&& (incnrmk[imn - 1] > incnrmk[imn - P._nsteps_check - 1])){
						DIV_MN_MSG (tstep, t, in, imn, whichone, incnrmk, incrk, resall, nsteps_check);
						rejectmnewton = true;
						convergedmnewton = false;
						rejmn++;
						break;
					}		  
		  
					std::cout<<"incr ("<<whichone<<") = "<<incrk;
					std::cout<<"   residual = [ "<<resall[0]<<" "<<resall[1]<<" ] * "<<std::endl;

					if (incrk < P._toll){
						CONV_MSG (tstep, in, imn, t, "incr", whichone, incrk, resall);
						rejectmnewton = false;
						convergedmnewton = true;
						break;
					}
				} /// END MODIFIED NEWTON
				
				incnrmk.resize(iimn);
				inck_clamp.resize(iimn);
				resnrmk.resize(iimn);
				
				totmn += iimn;
        
				if (reject){
					break;
				}

				if (rejectmnewton == true){
					V2 = V1;
					n2 = n1;
				}
				else{
					if (convergedmnewton == true){
						break;
					}
				}
				
				if (in >= P._maxit){
					std::cout<<"maximum number of Newton iterations reached,"<<std::endl;
					std::cout<<"try reducing timestep..."<<std::endl;
					if (dt > dtmin){
						reject = true;
						break;
					}
				}
			} /// END NEWTON STEP
			
			incnrm.resize(in);
			inc_clamp.resize(in);
			resnrm.resize(in);
			
			if (reject){
				++rejected;
				tstep -= 1;
				t = _tout [1];	// riporto t al passo told
				dt = dt * P._dtcut;

				std::cout<<"reverting to time step "<<tstep<<std::endl;
				std::cout<<"reducing time step: ";
				std::cout<<"model time "<<t<<" s, ";
				std::cout<<"new dt "<<dt<<" s"<<std::endl;
			}
			else{
				// aggiorno la soluzione al passo corrente con quella appena calcolata
				Voldold = Vold;
				noldold = nold;
				
				Vold = V2;
				nold = n2;
				
				_V = V2;
				_n = n2;
				
				_tout[0] = _tout[1];
				_tout[1] = t;
		
				dtfact = std::fmin(0.8 * sqrt(P._maxnpincr / incr0), P._maxdtincr);
				dt = std::fmax( std::fmin(dtfact * dt, dtmax), dtmin);

				std::cout<<" "<<std::endl;
				std::cout<<"t = "<<t<<" ";
				std::cout<<"<= "<<*std::max_element(tspan.begin(),tspan.end())<<" ;";
				std::cout<<" dtfact = "<<dtfact<<", estimate for next time step size: ";
				std::cout<<"dt = "<<dt<<" ";
				std::cout<<"(incr0 = "<<incr0<<")"<<std::endl;
			}

			V0.clear(); n0.clear();
		} /// END TIME STEP 
		
		saveRES(V2, "V2");
		saveRES(n2, "n2");
		
		if (P._savedata)
		{
			lastsaved = tstep;

			told = _tout[1];
			nsaves++;

			_V = V2;
			_n = n2;
			_tout[2] = t;
			
			saveNEWT(Vold, nold, told, V2, n2, _res, t, dt, nsaves, newton_solves, modified_newton_solves, freq);

			Vold.clear();
			nold.clear();
		}
	// //} /// END TIME FIXED SPAN
	
	std::cout<<"total number of rejected time steps: "<<rejected<<std::endl;
};