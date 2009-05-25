//__INSERT_LICENSE__
// $Id: smoke.cpp,v 1.7 2003/12/22 02:11:37 mstorti Exp $

#include "./smoke.h"

smoke_ff::~smoke_ff() { tmp.clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff::start_chunk(int &ret_options) {
  int ierr;
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);

  // Get element integer props
  elemset->elem_params(nel,ndof,nelprops);

  //o Frequency of oscillating local source. 
  EGETOPTDEF_ND(elemset,double,omega,0.);
  //o Coefficient scaling the reaction 
  EGETOPTDEF_ND(elemset,double,Cr,0.);
  //o Equilibrium value
  EGETOPTDEF_ND(elemset,double,phieq,1.);
  //o Dimension of problem
  EGETOPTDEF_ND(elemset,int,ndim,0);
  PETSCFEM_ASSERT0(ndim>0,"Dimension must be positive.");  

  //o Use nodal velocities
  EGETOPTDEF_ND(elemset,int,use_nodal_vel,0);

  //o Index of velocity in the H fields
  EGETOPTDEF_ND(elemset,int,nodal_vel_indx,1);

  if (!use_nodal_vel) {
    elemset->get_prop(u_prop,"u");
    PETSCFEM_ASSERT0(u_prop.length==ndim,
                     "u property must have length ndim");  
  }

  elemset->get_prop(G_prop,"G");
  PETSCFEM_ASSERT0(G_prop.length==0 || G_prop.length==2,
                   "G property, if given, must have length 2");  
  u.resize(1,ndim);

  //o Diffusivity
  EGETOPTDEF(elemset,double,diffusivity,0.0);
  PETSCFEM_ASSERT(diffusivity>=0.0,
                  "Diffusivity must be non-negative entereed %g",
                  diffusivity);  
  diff_max = diffusivity;

  //o Constant diffusivity
  EGETOPTDEF_ND(elemset,double,diffusivity0,0.0);
  PETSCFEM_ASSERT(diffusivity0>=0,
                  "diffusivity0 should be non-negative, entered %g",
                  diffusivity0);  

  //o Factor affecting stabilization term
  EGETOPTDEF_ND(elemset,double,tau_fac,1.0);

  // Tell `advdife' that we will use a scalar `tau'
  ret_options |= SCALAR_TAU;
  Cp.resize(2,ndof,ndof).eye();
  W_N.resize(2,nel,nel);
  A.resize(3,ndim,ndof,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::element_hook(ElementIterator &element) {
#if 0
  element_m = element;
  if (!use_nodal_vel)  
    u.set(elemset->prop_array(element_m,u_prop));
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff::set_state(const FastMat2 &UU) { 
  U.set(UU);
  phi = U.get(1); 
  diff = diff_max*mind(2,fabs(phi)/phieq,1.0) + diffusivity0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_A_jac_n(FastMat2 &A_jac_n, 
               FastMat2 &normal) {
  tmp3.prod(u,normal,-1,-1);
  double un = double(tmp3);
  A_jac_n.eye(un);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_A_grad_N(FastMat2 & A_grad_N,
                FastMat2 & grad_N) {
  A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff::compute_flux(COMPUTE_FLUX_ARGS) {
  if (use_nodal_vel) {
    int n1 = nodal_vel_indx;
    int n2 = n1 + ndim - 1;
    PETSCFEM_ASSERT0(H.dim(1)>=n2,"Not enough components in H field");  
    H.is(1,n1,n2);
    u.set(H);
    H.rs();
  } else {
    u.set(elemset->prop_array(element_m,u_prop));
  }
  A.ir(2,1).ir(3,1).set(u).rs();
  double vel = sqrt(u.sum_square_all());
  double G = 0.0;
  if (G_prop.length!=0) {
    assert(new_adv_dif_elemset);
    double t = new_adv_dif_elemset->time();
    const double *GG 
      = elemset->prop_array(element_m,G_prop);
    G = GG[0] * sin(omega*t) + GG[1] * cos(omega*t);
  }
  // Convective flux
  // flux(j,mu) = A(j,mu,nu) * U(nu)
  flux.prod(A,U,2,1,-1,-1);
  // Diffusive flux
  // fluxd(j,mu) = D(j,k,mu,nu) * grad_U(k,nu)
  fluxd.t().set(grad_U).scale(diff).rs();
  // A_grad_U(mu) = A(j,mu,nu) * grad_U(j,nu)
  A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
  // Reaction term
  double phi = U.get(1);
  double r = -Cr*(phi*phi-phieq*phieq)*phi;
  // scalar jacobian of reaction term 
  drdphi = -Cr*(3.0*phi*phi-phieq*phieq);
  // G_source.set(G).axpy(U,Cjac);
  G_source.set(G+r);
  // Set to zero
  tau_supg.set(0.);
  // No shock capturing
  delta_sc = 0.;
  // maximum eigenvlue = absolute value of velocity
  lam_max = vel;
  if (options & COMP_UPWIND) {
    // Intrinsic velocity
    Uintri.prod(iJaco,u,1,-1,-1);
    // This has scale of U/h, i.e. 1/T
    double tau, 
      Uh = sqrt(Uintri.sum_square_all()),
      h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());

    if (vel*vel > 20*Uh*diff) { 
      // remove singularity when D=0
      tau = 1.0/Uh;
    } else if (vel*vel > 1e-5*Uh*diff) {		
      double Pe  = vel*vel/(Uh*diff);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = 1.0/Uh*magic;
    } else {
      // remove singularity when v=0
      tau = h*h/(12.*diff);
    }

    // Set tau_(1,1) = scalar tau
    tau_supg.setel(tau_fac*tau,1,1);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
                       FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.ir(2,1).ir(4,1)
    .prod(dshapex,dshapex,-1,1,-1,2).scale(w*diff).rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_N_N_C(FastMat2 &N_N_C,
             FastMat2 &N,double w) {
  N_N_C.ir(2,1).ir(4,1).prod(N,N,1,2)
    .scale(-w*drdphi).rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
             FastMat2 &N,double w) {
  N_P_C.prod(P_supg,N,1,3,2).scale(-w*drdphi);
#if 0
  N_P_C.ir(3,1).ir(4,1).prod(N,P_supg,1,2).scale(-w*drdphi).rs();
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff::enthalpy(FastMat2 &H) {  H.set(U); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,
              const FastMat2 &N, double w) {
  W_N.prod(W,N,1,2).scale(w);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void smoke_ff
::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

