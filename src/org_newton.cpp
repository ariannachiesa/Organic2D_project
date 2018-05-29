/*! \file org_newton.cpp
  \brief Implementation of class Newton's methods
*/

#include "newton.h"

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