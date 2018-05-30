/*! \file newton.cpp
  \brief Class Newton
*/

#include "newton.h"

Newton::Newton(	Probl& P, std::vector<double>& Vin, std::vector<double>& nin, 
				std::vector<double>& tspan, std::vector<double>& Fin, std::vector<double>& Iin, BCS_CIRC& bcs, double freq)
{
	int	nnodes = P.get_msh_nodes(),
		newton_solves = 0,
        modified_newton_solves = 0,
		rejected = 0,
		j = 0,
		clamp = 0,
		lastsaved = 0,
		nsaves = 0,
		firstfixtstep, Nextvars, NI, tstep, in,
		whichone, totmn, iimn, told;

	double	t, dt, y0, yend, dtfact, rejmn,
			dt0 = P._dt0,
			dtmax = P._dtmax,
			dtmin = P._dtmin,
			tauk = 0,
			incr0, incrlast, incr1, incrk,
			incr0v = 0, incr0n = 0, incr0F = 0, incr0I = 0,
			incr1v = 0, incr1n = 0, incr1F = 0, incr1I = 0,
			incrkv = 0, incrkn = 0, incrkF = 0, incrkI = 0;
			
	bool	reject, rejectmnewton, convergedmnewton;
	
	auto quad = P._msh.begin_quadrant_sweep();
	
	std::vector<double>	F(4,0.0),
						delta, resnrm, resnrmk,
						inc_clamp, inck_clamp, incnrm, incnrmk,
						V0, n0, F0, I0,	
						V1, n1, F1, I1,
						V2, n2, F2, I2,
						Vk, nk, Fk, Ik,
						Vold, nold, Iold, Fold,
						dV, dn, dF, dI,
						resall(4,0),
						vecincr(4,0);
			
	// assert (numel (unique (device.pins)) == 2);
	// penso possa essere evitato come controllo
			
	dt = dt0;
	dt = std::fmin(dt, dtmax);
	
	/// INIT STATE FROM SCRATCH
	tstep = 0;
	_tout[0] = tspan[0];
    t = tspan[0];
	
	// cambiare la struttura vector di vector!!!
	_V.resize(Vin.size());
	_n.resize(nin.size());
	for(unsigned i=0; i<Vin.size(); i++){
		_V[i][0] = Vin[i];
	}
	for(unsigned i=0; i<nin.size(); i++){
		_n[i][0] = nin[i];
	}
	
	_F.resize(4);
	bcs.assign(t, P._Csb, P._Vshift, F);
	
	F = bcs.get_F();
    Nextvars  = F.size();
	if( std::accumulate( Fin.begin(), Fin.end(), 0.0) != 0.0){
		for(unsigned i=0; i<F.size(); i++){
			_F[i][0] = Fin[i];
		}
	}
	else{
		for(unsigned i=0; i<F.size(); i++){
			_F[i][0] = F[i];
		}
	}
    
	NI  = P._pins.size();
	_I.resize(NI);
	
	if( std::accumulate( Iin.begin(), Iin.end(), 0.0) != 0.0){
		for(unsigned i=0; i<Iin.size(); i++){
			_I[i][0] = Iin[i];
		}
	}

    firstfixtstep = 2;
	
	y0 = quad->p(1,0);
	std::cout<<"Connecting contact at y = "<<y0<<" to circuit pin "<<P._pins[0]<<std::endl;

	for (auto quadrant = P._msh.begin_quadrant_sweep (); quadrant != P._msh.end_quadrant_sweep (); ++quadrant){
		yend = quadrant->p(1,3);
	}
	std::cout<<"Connecting contact at y = "<<yend<<" to circuit pin "<<P._pins[1]<<std::endl;
	
	/// node ordering
	std::vector<int>    indexingV(nnodes,0),
						indexingn(nnodes,0),
						indexingF(Nextvars,0),
						indexingI(NI,0);

	for (int i=0; i<2*nnodes; i=i+2){
		indexingV[j] = i;
		j++;
	}	
	j=0;
	for (int i=1; i<2*nnodes; i=i+2){
		indexingn[j] = i;
		j++;
	}
	for (int i=0; i<Nextvars; i++){
		indexingF[i] = i + 2 * nnodes;
	}
	for (int i=0; i<NI; i++){
		indexingI[i] = i + Nextvars + 2 * nnodes;
	}
	
	for (unsigned n_fix_tstep = firstfixtstep; n_fix_tstep<tspan.size()+1; n_fix_tstep++){   /// TIME FIXED SPAN
		t = tspan[n_fix_tstep - 2];
		
		while (t < tspan[n_fix_tstep-1]){ /// TIME STEP
			std::cout<<"--"<<std::endl;
		
			tstep++;
			// std::cout<<"tstep = "<<tstep<<std::endl;

			_tout[tstep] = std::fmin (t + dt, tspan[n_fix_tstep - 1]);
			//std::cout<<"tout[tstep] = "<<_tout[tstep]<<std::endl;
			t = _tout[tstep];
			//std::cout<<"t = "<<t<<std::endl;
			dt = t - _tout[tstep - 1];
			//std::cout<<"dt = "<<dt<<std::endl;
			
			/// tutti gli incr?
			
			incr0 = 4 * P._maxnpincr;
			
			for(unsigned i=0; i<F.size(); i++){
				F[i] = _F[i][1];
			}
			bcs.assign(t, P._Csb, P._Vshift, F);	// ma quando sono alla prima iterazione ? passo una roba non inizializzata
			
			org_secs_state_predict (P, _V, _n, _F, _I, tstep, _tout, V0, n0, F0, I0);
			V2 = V0;
			n2 = n0;
			F2 = F0;
			I2 = I0;
			
			in = 0;
			//rejmn = 0;
			//totmn = 0;

			reject = false;
			
			// incnrm.resize(maxit);
			// inc_clamp.resize(maxit);
			// resnrm.resize(maxit); ???
			
			while (!reject && (in < P._maxit)){ /// NEWTON STEP
				in += 1;

				V1 = V2;
				n1 = n2;
				F1 = F2;
				I1 = I2;

				bcs.assign(t, P._Csb, P._Vshift, F2);
				
				// _res = org_secs2d_newton_residual(	P, V2, n2, F2, I2, _V[tstep-1], _n[tstep-1], _F[tstep-1], _I[tstep-1],
														// dt, bcs, indexingV, indexingn, indexingF, indexingI);
				// /// ora sembra essere ok
				// //for(unsigned i=0; i<_res.size(); i++){
				// //	std::cout<<"res = "<<_res[i]<<std::endl;
				// //}
				
				if (in == 1){
				
					whichone = 0;
					// compute_residual_norm (resnrm[0],whichone,resall,_res, P.alg(),indexingV,indexingn,indexingF,indexingI);
					
					//for(unsigned i=0; i<resall.size(); i++){
					//	std::cout<<"resall = "<<resall[i]<<std::endl;
					//}
					
					if(P._rowscaling.size() == resall.size()){
						for(unsigned i=0; i<P._rowscaling.size(); i++){
							P._rowscaling[i] = P._rowscaling[i] * (resall[i]+1);
						}
					}
					else{
						std::cout<<"error: dimensions mismatch"<<std::endl;
						break;
					}
					
					// _res = org_secs2d_newton_residual(	P, V2, n2, F2, I2, _V[tstep-1], _n[tstep-1], _F[tstep-1], _I[tstep-1],
															// dt, bcs, indexingV, indexingn, indexingF, indexingI);
				}

				resall.resize(4);
				// compute_residual_norm (resnrm[in-1],whichone,resall,_res, P.alg() ,indexingV,indexingn,indexingF,indexingI);
				
				// org_secs2d_newton_jacobian(	P,	V2, n2, F2,	dt, bcs, indexingV, indexingn, indexingF, indexingI, _jac);
				
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
		
				// if(in == 1){	// factorization of the matrix only at the first iteration of the Newton method
					mumps_solver.set_lhs_structure (_jac.rows(), irow, jcol);
					mumps_solver.analyze ();
				// }
		
				mumps_solver.set_lhs_data (vals);
				mumps_solver.set_rhs (delta);
				mumps_solver.factorize ();
		
				mumps_solver.solve ();
		
				for(unsigned i=0; i<delta.size(); i++){
					delta[i] *= (-1);
					// std::cout<<"delta = "<<delta[i]<<std::endl;	// non viene. rimane uguale a res ... ???
				}
				
				newton_solves +=1;
	
				dV.resize(indexingV.size());
				for(unsigned i=0; i<indexingV.size(); i++){
					dV[i] = delta[indexingV[i]] * P._colscaling[0];			// decisamente sbagliato
					//std::cout<<"dV = "<<dV[i]<<std::endl;
				}

				dn.resize(indexingn.size());
				for(unsigned i=0; i<indexingn.size(); i++){
					dn[i] = delta[indexingn[i]] * P._colscaling[1];			// decisamente sbagliato
					//std::cout<<"dn = "<<dn[i]<<std::endl;
				}

				dF.resize(indexingF.size());
				for(unsigned i=0; i<indexingF.size(); i++){
					dF[i] = delta[indexingF[i]] * P._colscaling[2];			// decisamente sbagliato
					//std::cout<<"dF = "<<dF[i]<<std::endl;
				}

				dI.resize(indexingI.size());
				for(unsigned i=0; i<indexingI.size(); i++){
					dI[i] = delta[indexingI[i]] * P._colscaling[3];			// decisamente sbagliato
					//std::cout<<"dI = "<<dI[i]<<std::endl;
				}
				delta.clear();
				
				V2.clear(); n2.clear(); F2.clear(); I2.clear();
				// org_secs_safe_increment (   V1, n1, F1, I1, dV, dn, dF, dI, P.alg(), P.dev(), P.cnst(),
												// // output
												// V2, n2, F2, I2, clamp, tauk);
												
				if ((clamp <= 0) || (tauk <= 0)){
					reject = true;
					break;
				}
		
				// compute_variation ( V0, n0, F0, I0, V2, n2, F2, I2, P.alg(), P.dev(), P.mat(), P.cnst(), std::fmin (clamp, tauk),
									// // output
									// incr0v, incr0n, incr0F, incr0I);
									
				vecincr[0] = incr0v;
				vecincr[1] = incr0n;
				vecincr[2] = incr0F;
				vecincr[3] = incr0I;

				incr0 = *std::max_element(vecincr.begin(),vecincr.end());
				if(incr0==incr0v) whichone = 0;
				if(incr0==incr0n) whichone = 1;
				if(incr0==incr0F) whichone = 2;
				if(incr0==incr0I) whichone = 3;
				
				if (incr0 > P._maxnpincr && dt > P._dtmin){
					// MAXINCR_MSG (tstep, t, in, whichone, incr0, resall, P.alg());
					reject = true;
					// infowhyfinished[tstep-1] = -3;
					break;
				}
	
				// compute_variation (	V1, n1, F1, I1, V2, n2, F2, I2, P.alg(), P.dev(), P.mat(), P.cnst(), std::fmin (clamp, tauk),
										// // output
										// incr1v, incr1n, incr1F, incr1I);
										
				vecincr[0] = incr1v;
				vecincr[1] = incr1n;
				vecincr[2] = incr1F;
				vecincr[3] = incr1I;

				incr1 = *std::max_element(vecincr.begin(),vecincr.end());
				if(incr1==incr1v) whichone = 0;
				if(incr1==incr1n) whichone = 1;
				if(incr1==incr1F) whichone = 2;
				if(incr1==incr1I) whichone = 3;
				
				//incnrm[in-1] = incr1;
				//inc_clamp[in-1] = incr1 / clamp;
				incrlast = incr1;
	
				if (resnrm[in-1] < P._toll ){
					//CONV_MSG (tstep, in, 1, t, "res", whichone, incr1, resall);
					// infowhyfinished[tstep-1] = 1;
					break;
				}

				if (in > P._nsteps_check && incnrm[in-1] > incnrm[in - P._nsteps_check -1]	//
					&& resnrm[in-1] > resnrm[in - P._nsteps_check -1] && dt > dtmin){
					// DIV_MSG (tstep, t, in, whichone, incnrm, incr1, resall, nsteps_check);
					reject = true;
					// infowhyfinished[tstep-1] = -4;
					break;
				}

				std::cout<<"incr ("<<whichone<<") = "<<incr1<<std::endl;
				std::cout<<"residual = [ "<<resall[0]<<" "<<resall[1]<<" "<<resall[2]<<" "<<resall[3]<<" ] "<<std::endl;
				
				if (incr1 < P._toll ){								/// è sbagliato ma per lo meno entra nell'if giusto (per ora)
					// CONV_MSG(tstep, in, 1, t, "incr", whichone, incr1, resall);
					// infowhyfinished[tstep-1] = 2;
					break;
				}
				
				/// MODIFIED NEWTON
				V1 = V2;
				n1 = n2;
				F1 = F2;
				I1 = I2;

				rejectmnewton = false;
				convergedmnewton = false;
				
				// incnrmk.resize(P._maxit_mnewton+1);
				// inck_clamp.resize(P._maxit_mnewton+1);
				// resnrmk.resize(P._maxit_mnewton+1);
				for (int imn = 1; imn<P._maxit_mnewton+1; imn++) {	/// MODIFIED NEWTON STEP 
					iimn = imn;
		
					Vk = V2;
					nk = n2;
					Fk = F2;
					Ik = I2;
		
					bcs.assign(t, P._Csb, P._Vshift, F2);
	  
					// _res = org_secs2d_newton_residual(	P, V2, n2, F2, I2, _V[tstep-1], _n[tstep-1], _F[tstep-1], _I[tstep-1],	//
														// dt, bcs, indexingV, indexingn, indexingF, indexingI);
														
					resall.resize(4);
					// compute_residual_norm (resnrmk[imn-1], whichone, resall, _res, P.alg(), indexingV, indexingn, indexingF, indexingI);
				
					/// Solve non linear system.
					std::cout << "Solving linear system."<<std::endl;
		
					delta = _res;
					mumps_solver.set_rhs (delta);
					
					mumps_solver.solve ();
	  
					for(unsigned i=0; i<delta.size(); i++){
						delta[i] *= (-1);
					}
	  
					modified_newton_solves +=1;
					
					for(unsigned i=0; i<indexingV.size(); i++){
						dV[i] = delta[indexingV[i]] * P._colscaling[0];
					}
					for(unsigned i=0; i<indexingn.size(); i++){
						dn[i] = delta[indexingn[i]] * P._colscaling[1];
					}
					for(unsigned i=0; i<indexingF.size(); i++){
						dF[i] = delta[indexingF[i]] * P._colscaling[2];
					}
					for(unsigned i=0; i<indexingI.size(); i++){
						dI[i] = delta[indexingI[i]] * P._colscaling[3];
					}

					V2.clear(); n2.clear(); F2.clear(); I2.clear();	  
					// org_secs_safe_increment (	Vk, nk, Fk, Ik, dV, dn, dF, dI, P.alg(), P.dev(), P.cnst(),
												// // output
												// V2, n2, F2, I2, clamp, tauk);
					
					if ((tauk <=0) || (clamp <= 0)){
						reject = true;
						// infowhyfinished[tstep-1] = -5;
						break;
					}
		  
					// compute_variation(	Vk, nk, Fk, Ik, V2, n2, F2, I2, P.alg(), P.dev(), P.mat(), P.cnst(), 
										// std::fmin (clamp, tauk),
										// // output
										// incrkv, incrkn, incrkF, incrkI);
		
					vecincr[0] = incrkv;
					vecincr[1] = incrkn;
					vecincr[2] = incrkF;
					vecincr[3] = incrkI;
		
					incrk = *std::max_element(vecincr.begin(),vecincr.end());

					if(incrk==incrkv) whichone = 0;
					if(incrk==incrkn) whichone = 1;
					if(incrk==incrkF) whichone = 2;
					if(incrk==incrkI) whichone = 3;
					
					//incnrmk[imn - 1] = incrk;
					//inck_clamp[imn - 1] = incrk / clamp;
					incrlast = incrk;

					if (resnrmk[imn - 1] < P._toll){
						// CONV_MSG (tstep, in, imn, t, "res", whichone, incrk, resall);
						// infowhyfinished[tstep-1] = 3;
						rejectmnewton = false;
						convergedmnewton = true;
						break;
					}

					if (imn > P._nsteps_check && (resnrmk[imn - 1] > resnrmk[imn - P._nsteps_check - 1])
						&& (incnrmk[imn - 1] > incnrmk[imn - P._nsteps_check - 1])){
						// DIV_MN_MSG (tstep, t, in, imn, whichone, incnrmk, incrk, resall, nsteps_check);
						rejectmnewton = true;
						convergedmnewton = false;
						rejmn++;
						break;
					}		  
		  
					std::cout<<"incr ("<<whichone<<") = "<<incrk;
					std::cout<<"residual = [ "<<resall[0]<<" "<<resall[1]<<" "<<resall[2]<<" "<<resall[3]<<" ] * "<<std::endl;

					if (incrk < P._toll){
						// CONV_MSG (tstep, in, imn, t, "incr", whichone, incrk, resall);
						// infowhyfinished[tstep-1] = 4;
						rejectmnewton = false;
						convergedmnewton = true;
						break;
					}
				} /// END MODIFIED NEWTON
				
				//incnrmk.resize(iimn);
				//inck_clamp.resize(iimn);
				//resnrmk.resize(iimn);
				
				totmn += iimn;
        
				if (reject){
					break;
				}

				if (rejectmnewton == true){
					V2 = V1;
					n2 = n1;
					F2 = F1;
					I2 = I1;
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
						// infowhyfinished[tstep-1] = -6;
						break;
					}
				}
			} /// END NEWTON STEP
			
			// incnrm.resize(in);
			// inc_clamp.resize(in);
			// resnrm.resize(in);
			
			if (reject){
				++rejected;
				tstep -= 1;
				t = _tout [tstep];
				dt = dt * P._dtcut;

				std::cout<<"reverting to time step "<<tstep<<std::endl;
				std::cout<<"reducing time step: ";
				std::cout<<"model time "<<t<<" s, ";
				std::cout<<"new dt "<<dt<<" s"<<std::endl;
			}
			else{

				for(unsigned i=0; i<V2.size(); i++){
					_V[i][1] = V2[i];
				}
				for(unsigned i=0; i<n2.size(); i++){
					_n[i][1] = n2[i];
				}
				for(unsigned i=0; i<F2.size(); i++){
					_F[i][1] = F2[i];
				}
				for(unsigned i=0; i<I2.size(); i++){
					_I[i][1] = I2[i];
				}
		
				dtfact = std::fmin(0.8 * sqrt(P._maxnpincr / incr0), P._maxdtincr);
				dt = std::fmax( std::fmin(dtfact * dt, dtmax), dtmin);

				std::cout<<" "<<std::endl;
				std::cout<<"t = "<<t<<" ";
				std::cout<<"<= "<<*std::max_element(tspan.begin(),tspan.end())<<" ;";
				std::cout<<" dtfact = "<<dtfact<<", estimate for next time step size: ";
				std::cout<<"dt = "<<dt<<" ";
				std::cout<<"(incr0 = "<<incr0<<")"<<std::endl;
			}
			
			if (P._savedata && ((tstep - lastsaved) >= 100)){
				lastsaved = tstep;
		
				Vold.resize(_V.size());
				for(unsigned i=0; i<Vold.size(); i++){
					Vold[i] = _V[i][0];
				}
		
				nold.resize(_n.size());
				for(unsigned i=0; i<nold.size(); i++){
					nold[i] = _n[i][0];
				}
		
				Fold.resize(_F.size());
				for(unsigned i=0; i<Fold.size(); i++){
					Fold[i] = _F[i][0];
				}
		
				Iold.resize(_I.size());
				for(unsigned i=0; i<Iold.size(); i++){
					Iold[i] = _I[i][0];
				}

				told = _tout[tstep - 1];
	  
				// saveNEWT(	Vold, nold, Fold, Iold, told, V2, n2, F2, I2, _res, t, dt, 
							// nsaves, newton_solves, modified_newton_solves, freq);
				
				Vold.clear();
				nold.clear();
				Fold.clear();
				Iold.clear();
			}

			V0.clear(); n0.clear(); F0.clear(); I0.clear();
		} /// END TIME STEP 
		
		if (P._savedata)
		{
			lastsaved = tstep;

			Vold.resize(_V.size());
			for(unsigned i=0; i<Vold.size(); i++){
				Vold[i] = _V[i][0];
			}
		
			nold.resize(_n.size());
			for(unsigned i=0; i<nold.size(); i++){
				nold[i] = _n[i][0];
			}
		
			Fold.resize(_F.size());
			for(unsigned i=0; i<Fold.size(); i++){
				Fold[i] = _F[i][0];
			}
		
			Iold.resize(_I.size());
			for(unsigned i=0; i<Iold.size(); i++){
				Iold[i] = _I[i][0];
			}
	  
			told = _tout[tstep - 1];
			nsaves++;
	  
			// saveNEWT(	Vold, nold, Fold, Iold, told, V2, n2, F2, I2, _res, t, dt, 
						// nsaves, newton_solves, modified_newton_solves, freq);

			Vold.clear();
			nold.clear();
			Fold.clear();
			Iold.clear();
		}
	} /// END TIME FIXED SPAN
	
	std::cout<<"total number of rejected time steps: "<<rejected<<std::endl;
};