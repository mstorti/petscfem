//__INSERT_LICENSE__
// $Id: smoke.cpp,v 1.7 2003/12/22 02:11:37 mstorti Exp $

#include "./smoke.h"

smoke_ff::~smoke_ff() { tmp.clear(); }

void smoke_ff::start_chunk(int &ret_options) {
  int ierr;
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);     

  // Get element integer props
  elemset->elem_params(nel,ndof,nelprops);

  //o Frequency of oscillating local source. 
  EGETOPTDEF_ND(new_adv_dif_elemset,double,omega,0.);
  //o Coefficient scaling the reaction 
  EGETOPTDEF_ND(new_adv_dif_elemset,double,Cr,0.);
  //o Equilibrium value
  EGETOPTDEF_ND(new_adv_dif_elemset,double,phieq,1.);
  //o Dimension of problem
  EGETOPTDEF_ND(new_adv_dif_elemset,int,ndim,0);
  PETSCFEM_ASSERT0(ndim>0,"Dimension must be positive.");  

  assert(omega>0.);

  elemset->get_prop(u_prop,"u");
  elemset->get_prop(G_prop,"G");
  assert(u_prop.length==ndim);
  assert(G_prop.length==2);
  u.resize(1,ndim);
  // Tell `advdife' that we will use a scalar `tau'
  ret_options |= SCALAR_TAU;
  Cp.resize(2,ndof,ndof).eye();
  W_N.resize(2,nel,nel);
  A.resize(3,ndim,ndof,ndof);
}

void smoke_ff::element_hook(ElementIterator &element) {
  element_m = element;
  u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
  A.ir(2,1).ir(3,1).set(u).rs();
}

void smoke_ff::set_state(const FastMat2 &UU) { 
  U.set(UU);
  phi = U.get(1); 
}

void smoke_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0);			// fixme:= Not implemented yet, used for
				// absorbing boundary conditions
}

void smoke_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
}

inline double modulo(double k, double n, int *div=NULL) {
  double m = fmod(k,n);
  int d = int((k-m)/n);
  if (m < 0.) {
    m += n;
    d -= 1;
  }
  if (div) *div = d;
  return m;
}

void smoke_ff::compute_flux(COMPUTE_FLUX_ARGS) {
  u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
  double uu = sqrt(u.sum_square_all());
  const double *GG = new_adv_dif_elemset->prop_array(element_m,G_prop);
  double t = new_adv_dif_elemset->time();
  double T = 2.0*M_PI/omega;
  double G = GG[0] * sin(omega*t) + GG[1] * cos(omega*t);
  // Convective flux
  // flux(j,mu) = A(j,mu,nu) * U(nu)
  flux.prod(A,U,2,1,-1,-1);
  // Diffusive flux
  // fluxd(j,mu) = D(j,k,mu,nu) * grad_U(k,nu)
  fluxd.set(0.);
  // A_grad_U(mu) = A(j,mu,nu) * grad_U(j,nu)
  A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
  // Reaction term
  double phi = U.get(1);
  double r2 = -Cr*(phi*phi-phieq*phieq)*phi;
  double r1 = -Cr*phi;
  double alpha = (1+cos(omega*t))/2.;
  alpha = 1.;
  double r = alpha*r2+(1.0-alpha)*r1;
  // scalar jacobian of reaction term 
  double drdphi2 = -Cr*(2.*phi*phi-phieq*phieq);
  double drdphi1 = -Cr;
  drdphi = alpha*drdphi2+(1.0-alpha)*drdphi1;
  // G_source.set(G).axpy(U,Cjac);
  G_source.set(G+r);
  // Set to zero
  tau_supg.set(0.);
  // No shock capturing
  delta_sc = 0.;
  // maximum eigenvlue = absolute value of velocity
  lam_max = uu;
  // Intrinsic velocity
  Uintri.prod(iJaco,u,1,-1,-1);
  // This has scale of U/h
  double Uh = sqrt(Uintri.sum_square_all());

  // Consider pure advective flow
  double magic = 1., tau_fac=1.;

  // Set tau_(1,1) = scalar tau
  tau_supg.setel(tau_fac/Uh*magic,1,1);
}

void smoke_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
}

void smoke_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.ir(2,1).ir(4,1).prod(N,N,1,2).scale(-w*drdphi).rs();
}

void smoke_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			  FastMat2 &N,double w) {
  N_P_C.prod(P_supg,N,1,3,2).scale(-w*drdphi);
#if 0
  N_P_C.ir(3,1).ir(4,1).prod(N,P_supg,1,2).scale(-w*drdphi).rs();
#endif
}

void smoke_ff::enthalpy(FastMat2 &H) {  H.set(U); }

void smoke_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		 double w) {
  W_N.prod(W,N,1,2).scale(w);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

void smoke_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

