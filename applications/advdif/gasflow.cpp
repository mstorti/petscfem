//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.33.6.2 2005/03/27 22:06:39 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/generror.h>

#include "gasflow.h"

extern string fastmat2_stat_current_string;
#define FM2STAT(s) { fastmat2_stat_current_string = s; }

#define GF_GETOPTDEF_ND(type,var,def)				\
 { if (elemset) { EGETOPTDEF_ND(elemset,type,var,def); }	\
 else if (old_elemset)						\
      { TGETOPTDEF_ND(old_elemset->thash,type,var,def); }	\
 else assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::start_chunk(int &ret_options) {
  int ierr;
  FastMat2 tmp5;

  if (elemset) {
    new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
    elemset->elem_params(nel,ndof,nelprops);
  } else if (old_elemset) {
    nel = old_elemset->nel;
    ndof = old_elemset->ndof;
  } else assert(0);

  GF_GETOPTDEF_ND(int,ndim,0);
  GF_GETOPTDEF_ND(double,visco,0);
  GF_GETOPTDEF_ND(double,cond,0);

  // Turbulence parameters
  //o Activates LES turbulence model
  GF_GETOPTDEF_ND(int,LES,0);
  //o C_smag
  GF_GETOPTDEF_ND(double,C_smag,0.09);
  //o Turbulent Prandtl number
  GF_GETOPTDEF_ND(double,Pr_t,1.);
  //o Mask for stabilization term
  GF_GETOPTDEF_ND(double,tau_fac,1.);

  //o Threshold value for density. If density goes belows this value
  //  then a cutoff function is applied. This is a simple way to prevent
  //  crashing of the code due to taking square root of negative values.
  GF_GETOPTDEF_ND(double,rho_thrsh,0.);
  //o Threshold value for pressure. See doc for rho_thrsh.
  GF_GETOPTDEF_ND(double,p_thrsh,0.);
  //o If this flag is activated the code stops when a negative 
  //  value for density or pressure is found. Setting either
  //  #p_thrsh# or #rho_thrsh# deactivates #stop_on_neg_val# . 
  GF_GETOPTDEF_ND(int,stop_on_neg_val,1);
  PETSCFEM_ASSERT0(rho_thrsh>=0,"Density threshold should be non-negative");  
  PETSCFEM_ASSERT0(p_thrsh>=0,"Pressure threshold should be non-negative");  
  if (rho_thrsh>0 || p_thrsh>0) stop_on_neg_val=0;

  //o Adjust the stability parameters, taking into account
  // the time step. If the #steady# option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  // #temporal_stability_factor# is set to 0.
  GF_GETOPTDEF_ND(double,temporal_stability_factor,1.);

  // gamma coeficient
  GF_GETOPTDEF_ND(double,ga,1.4);
  // constant of a particular gas for ideal gas law (state equation for the gas)
  GF_GETOPTDEF_ND(double,Rgas,287.);

  // shock capturing scheme 
  // [0] standard [default]
  // [1] see Tezduyar & Senga WCCM VI (2004) 
  GF_GETOPTDEF_ND(int,shocap_scheme,0);
  // beta parameter for shock capturing, see Tezduyar & Senga WCCM VI (2004)
  GF_GETOPTDEF_ND(double,shocap_beta,1.0);

  // Sutherland law key
  GF_GETOPTDEF_ND(int,sutherland_law,0);
  // Sutherland law infinity Temperature
  GF_GETOPTDEF_ND(double,Tem_infty,0.0);
  // Sutherland law reference Temperature
  GF_GETOPTDEF_ND(double,Tem_ref,0.0);

  //o _T: double[ndim] _N: G_body _D: null vector
  // _DOC: Vector of gravity acceleration (must be constant). _END
#if 0
  if (elemset) {
    ierr = elemset->get_double("G_body",
			       *G_body.storage_begin(),1,ndim);
  }
#endif
  G_body.resize(1,ndim);
  G_body.set(0.);
  if (new_adv_dif_elemset) {
    new_adv_dif_elemset->get_prop(G_body_prop,"G_body");
    if (G_body_prop.length>0)
      assert(G_body_prop.length == ndim);
    new_adv_dif_elemset
      ->get_prop(G_body_scale_prop,"G_body_scale");
  }

  //o _T: double[ndof] _N: Uref _D: <none>
  // _DOC: reference state for absorbing b.c.'s. _END
  Uref.resize(1,ndof).set(0.);
  if (old_elemset) {
    const char *line;
    vector<double> Uref_v;
    old_elemset->thash->get_entry("G_body",line);
    if(line) {
      read_double_array(Uref_v,line);
    assert(Uref_v.size()==ndof);
    Uref.set(&Uref_v[0]);
    }
  } else {
    ierr = elemset->
      get_double("Uref",*Uref.storage_begin(),1,ndof);
  }
  dUabso.resize(1,ndof);

  //o Use linear approximation for the absorbing
  //  boundary condition. Gives fully absorbing b.c. in the
  //  linear range, but not truly Riemman Invariant in the
  //  nonlinear range. 
  GF_GETOPTDEF_ND(int,linear_abso,0);

  assert(ndim>0);
  assert(ndof==2+ndim);

  assert(ga>0.);
  assert(Rgas>0.);
  assert(visco>0.);
  assert(cond>0.);
  vl_indx = 2;
  vl_indxe = 2+ndim-1;
  vel.resize(1,ndim);
  vel_supg.resize(1,ndim);
  Cp.resize(2,ndof,ndof);
  Cpi.resize(2,ndof,ndof);
  Ajac.resize(3,ndim,ndof,ndof);
  Ajac_tmp.resize(3,ndim,ndof,ndof);
  Id_ndim.resize(2,ndim,ndim).eye();
  Amom.resize(2,ndim,ndim);
  Y.resize(2,ndim,ndim);

  Djac.resize(4,ndim,ndof,ndim,ndof);
  Djac_tmp.resize(4,ndim,ndof,ndim,ndof);
  tmp1.resize(4,nel,ndof,ndim,ndof);
  Cjac.resize(2,ndof,ndof);
  tmp2.resize(2,nel,nel);
  tmp3.resize(2,ndof,ndof);
  tmp20.resize(0);
  maktgsp.init(ndim);

  grad_vel.resize(2,ndim,ndim);
  strain_rate.resize(2,ndim,ndim);
  sigma.resize(2,ndim,ndim);
  grad_T.resize(1,ndim);
  grad_p.resize(1,ndim);

  IdId.resize(4,ndim,ndim,ndim,ndim);
  IdId.prod(Id_ndim,Id_ndim,1,3,2,4);
  //  IdId.prod(Id_ndim,Id_ndim,1,3,2,4);
  //  tmp5.prod(Id_ndim,Id_ndim,2,3,1,4);
  //  IdId.add(tmp5);

  uintri.resize(1,ndim);
  svec.resize(1,ndim);
  tmp9.resize(1,nel);
  W_N.resize(2,nel,nel);
  jvec.resize(1,ndim);

  vel_old.resize(1,ndim);

  tmp_vel.resize(1,ndim);

  if (sutherland_law>0) {
     assert(Tem_infty>0.0);
     assert(Tem_ref>0.0);
  }

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::element_hook(ElementIterator &element) {
  if (new_adv_dif_elemset) {
    const double *G_body_p
      = new_adv_dif_elemset
      ->prop_array(element,G_body_prop);
    G_body_scale = new_adv_dif_elemset
      ->prop_val(element,G_body_scale_prop,
		 new_adv_dif_elemset->time());
    if (G_body_prop.length>0)
      G_body.set(G_body_p).scale(G_body_scale);
    else G_body.set(0.);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
gasflow_ff::gasflow_ff(NewElemset *e)
  : AdvDifFFWEnth(e), old_elemset(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
gasflow_ff::gasflow_ff(Elemset *e)
  : old_elemset(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
gasflow_ff::~gasflow_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::set_state(const FastMat2 &UU) {
  U.set(UU);
  double dummy;
  // Scalar variables
  rho = U.get(1);
  p = U.get(ndof);

  if (stop_on_neg_val && (p<0 || rho<0)) 
    throw GenericError("negative pressure or density detected.");

  // Cutoff values
  rho = ctff(rho,dummy,rho_thrsh);
  p = ctff(p,dummy,p_thrsh);

  // Velocities
  U.is(1,vl_indx,vl_indxe);
  vel.set(U);
  U.rs();
  velmod = sqrt(vel.sum_square_all());
  g1=ga-1;
  rho_ene = p/g1+0.5*rho*square(velmod);
  ene = rho_ene/rho;
  entalpy=rho_ene+p;
  Cv=Rgas/g1;
  int_ene=p/g1/rho;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
set_state(const FastMat2 &U,
	  const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::enthalpy(FastMat2 &H) {

  H.set(0.);

  H.setel(rho,1);
  H.is(1,vl_indx,vl_indxe)
    .set(vel).scale(rho).rs();
  H.setel(rho_ene,ndof);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,
	      const FastMat2 &N,double weight) {
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
comp_P_Cp(FastMat2 &P_Cp,
	  const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
get_Cp(FastMat2 &Cp_a) {
  Cp_a.set(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
get_Ajac(FastMat2 &Ajac_a) {
  Ajac_a.set(Ajac);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
compute_tau(int ijob,double &delta_sc) {

  //    static FastMat2 svec,tmp9;

    const FastMat2 &grad_N = *advdf_e->grad_N();

    double tol=1.0e-16;
    h_supg=0;
    //    const FastMat2 &grad_N = *advdf_e->grad_N();
    velmod = sqrt(vel_supg.sum_square_all());
    FastMat2::branch();
    if(velmod>tol) {
      FastMat2::choose(0);
      // svec:= a streamline oriented unit vector
      svec.set(vel_supg).scale(1./velmod);
      h_supg = tmp9.prod(grad_N,svec,-1,1,-1).sum_abs_all();
      h_supg = (h_supg < tol ? tol : h_supg);
      h_supg = 2./h_supg;
    } else {
      h_supg = h_pspg;
    }
    FastMat2::leave();

    double Peclet = velmod * h_supg / (2. * visco_supg);
    double sonic_speed = sqrt(ga*p/rho);
    double velmax = velmod+sonic_speed;
    
    tau_supg_a = h_supg/2./velmax;

#if 0
    if (shocap_scheme==0) {
      double fz = (Peclet < 3. ? Peclet/3. : 1.);
      delta_sc = 0.5*h_supg*velmax*fz;
    } 
#endif
    // antes era shocap_scheme==1
    double tol_shoc = 1e-010;
    // compute j direction , along density gradient
    double h_shoc, grad_rho_mod = sqrt(grad_rho.sum_square_all());
    FastMat2::branch();
    if(grad_rho_mod>tol_shoc) {
      FastMat2::choose(0);
      //      grad_rho_mod = (grad_rho_mod <= 0 ? tol_shoc : grad_rho_mod);
      jvec.set(grad_rho).scale(1.0/grad_rho_mod);
      h_shoc = tmp9.prod(grad_N,jvec,-1,1,-1).sum_abs_all();
      h_shoc = (h_shoc < tol ? tol : h_shoc);
      h_shoc = 2./h_shoc;
    } else {
      FastMat2::choose(1);
      jvec.set(0.);
      h_shoc = h_supg;
    }
    FastMat2::leave();
    double fz = grad_rho_mod*h_shoc/rho;
    fz = pow(fz,shocap_beta);
    delta_sc_aniso = 0.5*h_shoc*velmax*fz;
    
    double fz2 = (Peclet < 3. ? Peclet/3. : 1.);
    delta_sc = 0.5*h_supg*velmax*fz2;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  double tmp00,tmp00_lam,tmp01,tmp02,tmp03,tmp04;

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau


  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Enthalpy Jacobian Cp

  // double velmod = sqrt(vel.sum_square_all());
  // double entalpy=rho_ene+p;

  Cp.set(0.);
  Cp.setel(1.0,1,1);
  Cp.is(1,vl_indx,vl_indxe).ir(2,1).set(vel).rs();
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(rho).rs();
  Cp.ir(1,ndof).is(2,vl_indx,vl_indxe).set(vel).scale(rho).rs();
  Cp.setel(1.0/g1,ndof,ndof);
  Cp.setel(0.5*square(velmod),ndof,1);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Advective fluxes
  flux.set(0.);
  flux.ir(1,1).set(vel).scale(rho).rs();
  Amom.prod(vel,vel,1,2);
  Y.set(0.).add(Amom).scale(rho).axpy(Id_ndim,p);
  flux.is(1,vl_indx,vl_indxe).add(Y).rs();
  flux.ir(1,ndof).set(vel).scale(entalpy).rs();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Adjective Jacobians
  Ajac.set(0.);

  // First row
  Ajac.ir(2,1).is(3,vl_indx,vl_indxe).set(Id_ndim).rs();
  // Last column
  Ajac.ir(3,ndof).is(2,vl_indx,vl_indxe).set(Id_ndim).scale(g1).rs();
  Ajac.ir(2,ndof).ir(3,ndof).set(vel).scale(ga).rs();

  // Momentum row and columns
  for (int j=1; j<=ndim; j++) {
    //  vel.ir(1,j);
  Y.set(0.).eye(vel.get(j));
  //  vel.rs();
  Y.ir(2,j).add(vel).rs();
  Y.ir(1,j).axpy(vel,-g1).rs();
  Ajac.ir(1,j).is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).add(Y).rs();
  }

  // Momentum rows and first column
  Y.eye(square(velmod)).scale(g1/2.);
  Ajac.ir(3,1).is(2,vl_indx,vl_indxe).add(Y).rs();
  Ajac.ir(3,1).is(2,vl_indx,vl_indxe).rest(Amom).rs();

  // Momentum columns and last row
  Ajac.ir(2,ndof).is(3,vl_indx,vl_indxe).rest(Amom).scale(g1);
  Y.eye(ene).scale(ga);
  Ajac.add(Y);
  Y.eye(square(velmod)).scale(g1/2.);
  Ajac.rest(Y).rs();

  // First column and last row
  tmp00 = g1*square(velmod)-ga*ene;
  Ajac.ir(2,ndof).ir(3,1).set(vel).scale(tmp00).rs();

  // Last column and last row
  Ajac.ir(2,ndof).ir(3,ndof).set(vel).scale(ga).rs();

  // A_i * Cp  transformed to the primitive state variables
  Ajac_tmp.set(Ajac);
  Ajac.prod(Ajac_tmp,Cp,1,2,-1,-1,3);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

  A_grad_U.prod(Ajac,grad_U,-1,1,-2,-1,-2);

  // Strain rate tensor
  grad_U.is(2,vl_indx,vl_indxe);
  grad_vel.set(grad_U);
  grad_U.rs();
  strain_rate.set(grad_vel);
  grad_vel.t();
  strain_rate.add(grad_vel).scale(0.5);
  grad_vel.rs();

  grad_vel.d(1,2);

  // REVISAR ESTO
  tmp10.sum(grad_vel,-1);
  grad_vel.rs();
  double divu = double(tmp10);
  // effective viscosity
  visco_t = 0.;  // aqui agregar la expresion algebraica para la misma
  cond_t = 0.;  // aqui agregar la expresion algebraica para la misma

  if (LES) {
    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    if (advdf_e) {
#define pi M_PI
      double Volume = advdf_e->volume();
      int axi = advdf_e->axi;
      double h_grid;
      
      if (ndim==2 | (ndim==3 && axi>0)) {
	h_grid = sqrt(4.*Volume/pi);
      } else if (ndim==3) {
	h_grid = cbrt(6*Volume/pi);
      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Only dimensions 2 and 3 allowed for this element.\n");
      }
      double strain_rate_abs = strain_rate.sum_square_all();
      
      strain_rate_abs = sqrt(strain_rate_abs);
      double nu_t = C_smag*strain_rate_abs*h_grid*h_grid;
      visco_t = rho*nu_t;
      double Cp = ga*Cv;
      cond_t = Cp*visco_t/Pr_t;
    }
  }

  // Sutherland law for thermal correction of laminar viscosity

  double sutherland_factor = 1.0;
  if (sutherland_law !=0) {
    sutherland_factor = pow(Tem/Tem_infty,1.5)
      *((Tem_infty+Tem_ref)/(Tem_ref+Tem));
  }
  double visco_l = sutherland_factor * visco;
  double visco_eff = (visco_l + visco_t);
  // effective thermal conductivity
  double cond_eff = sutherland_factor * cond + cond_t;
  // Aca ponemos visco (y no visco_eff). Ver paper de LESIEUR Y COMTE ...
  double visco_bar = visco_l - cond_eff/Cv;

  // Stress tensor
  sigma.set(0.);
  //  sigma.eye(divu).scale(-2./3.*visco_eff).add(strain_rate).scale(2.*visco_eff).rs();
  sigma.eye(divu).scale(-2./3.).add(strain_rate).scale(2.).rs();

  // Temperature gradient
  grad_U.ir(2,1);
  grad_rho.set(grad_U);
  grad_U.ir(2,ndof);
  grad_p.set(grad_U);
  grad_U.rs();
  grad_T.set(grad_p).axpy(grad_rho,-p/rho).scale(1./Rgas/rho).rs();

  //  LIMITAR POR VALORES NEGATIVOS DE RHO Y P , VISCO, CONDUCTIVITY , ETC

  // Diffusive fluxes
  fluxd.set(0.);
  fluxd.is(1,vl_indx,vl_indxe).set(sigma).scale(visco_eff).rs();
  // En el paper de LESIEUR Y COMTE usan viscosidad laminar en la ecuacion de energia
  viscous_work.prod(sigma,vel,1,-1,-1).scale(visco_l);
  heat_flux.set(grad_T).scale(-cond_eff);
  fluxd.ir(1,ndof).set(viscous_work).rest(heat_flux).rs();

  // Diffusive jacobians
  Djac.set(0.);
  Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,visco_eff).rs();
  tmp00 = 1./3.*visco_eff;
  tmp00_lam = 1./3.*visco_l;
  for (int j=1; j<=ndim; j++) {
     Djac.addel(tmp00,j,j+1,j,j+1);
  }
  // Momentum rows and first columns of K_ii
  tmp_vel.set(vel).scale(-tmp00);
  for (int j=1; j<=ndim; j++) {
    Djac.ir(1,j).ir(3,j).is(2,vl_indx,vl_indxe).ir(4,1).set(vel).scale(-visco_eff).rs();
    Djac.ir(1,j).ir(3,j).ir(4,1).ir(2,j+1).add(tmp_vel.get(j)).rs();
  }
  // Momentum columns and last row of K_ii
  tmp_vel.set(vel).scale(tmp00_lam);
  for (int j=1; j<=ndim; j++) {
    Djac.ir(1,j).ir(3,j).is(4,vl_indx,vl_indxe).ir(2,ndof).set(vel).scale(visco_bar).rs();
    Djac.ir(1,j).ir(3,j).ir(2,ndof).ir(4,j+1).add(tmp_vel.get(j)).rs();
  }

  // Last row and first and last column
  tmp02 = cond_eff/Cv;
  tmp01 = (0.5*square(velmod)-int_ene)*tmp02-visco_l*square(velmod);
  for (int j=1; j<=ndim; j++) {
    vel_j2 = square(double(vel.get(j)));
    //    vel.rs();
    tmp03 = tmp01 - tmp00_lam*vel_j2;
    Djac.setel(tmp02,j,ndof,j,ndof);
    Djac.setel(tmp03,j,ndof,j,1);
  }
  Djac.rs();

  // Diffusive Jacobians K_ij with i .ne. j
  tmp00 = 2.*visco_eff/3.;
  tmp00_lam = 2.*visco_l/3.;

  if(ndim==2) {
  int ip[] = {1,2,1};

  for (int j=0; j<ndim; j++) {

    Djac.setel(-tmp00,ip[j],ip[j]+1,ip[j+1],ip[j+1]+1);
    Djac.setel(visco_eff,ip[j],ip[j+1]+1,ip[j+1],ip[j]+1);

    tmp04 = -visco_l/3*double(vel.get(ip[j]))*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],1);

    tmp04 = tmp00*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ip[j]+1,ip[j+1],1);
    tmp04 = -visco_eff*double(vel.get(ip[j]));
    Djac.setel(tmp04,ip[j],ip[j+1]+1,ip[j+1],1);

    tmp04 = -tmp00_lam*double(vel.get(ip[j]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],ip[j+1]+1);
    tmp04 = visco_l*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],ip[j]+1);
  }

  } else {

  int ip[] = {1,2,3,1};

  for (int j=0; j<ndim; j++) {

    Djac.setel(-tmp00,ip[j],ip[j]+1,ip[j+1],ip[j+1]+1);
    Djac.setel(visco_eff,ip[j],ip[j+1]+1,ip[j+1],ip[j]+1);

    tmp04 = -visco_l/3*double(vel.get(ip[j]))*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],1);

    tmp04 = tmp00*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ip[j]+1,ip[j+1],1);
    tmp04 = -visco_eff*double(vel.get(ip[j]));
    Djac.setel(tmp04,ip[j],ip[j+1]+1,ip[j+1],1);

    tmp04 = -tmp00_lam*double(vel.get(ip[j]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],ip[j+1]+1);
    tmp04 = visco_l*double(vel.get(ip[j+1]));
    Djac.setel(tmp04,ip[j],ndof,ip[j+1],ip[j]+1);
  }

  }

  // D_ij * Cp  transformed to the primitive state variables
  Djac_tmp.set(Djac);
  Djac.prod(Djac_tmp,Cp,1,2,3,-1,-1,4);
  //  Djac.scale(1./rho);

  // Reactive terms
  //  there is no reactive terms in this formulation
  Cjac.set(0.);

  if (options & COMP_UPWIND) {
    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);
#define pi M_PI
    double Volume = advdf_e->volume();
    int axi = advdf_e->axi;

    if (ndim==2 | (ndim==3 && axi>0)) {
      h_pspg = sqrt(4.*Volume/pi);
    } else if (ndim==3) {
      h_pspg = cbrt(6*Volume/pi);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Only dimensions 2 and 3 allowed for this element.\n");
    }

    FastMat2 &Uo = (FastMat2 &) advdf_e->Uold();

    Uo.is(1,vl_indx,vl_indxe);
    vel_old.set(Uo);
    Uo.rs();

    const FastMat2 &grad_N = *advdf_e->grad_N();

    vel_supg.set(vel_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }

    visco_supg = visco_eff/rho;
    int ijob =0;
    compute_tau(ijob,delta_sc);

    if (tau_fac != 1.) {
      tau_supg_a *= tau_fac;
    }

    tau_supg.eye(tau_supg_a);

// postmultiply by inv(Cp)
    Cpi.set(0.);
    Cpi.setel(1.0,1,1);
    Cpi.is(1,vl_indx,vl_indxe).ir(2,1).set(vel).scale(-1.0/rho).rs();
    Cpi.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(1.0/rho).rs();
    Cpi.ir(1,ndof).is(2,vl_indx,vl_indxe).set(vel).scale(-g1).rs();
    Cpi.setel(g1,ndof,ndof);
    Cpi.setel(0.5*square(velmod)*g1,ndof,1);

//    tau_supg_c.set(tau_supg);
    tau_supg_c.prod(tau_supg,Cpi,1,-1,-1,2);

  }

  if (options & COMP_SOURCE) {
    G_source.set(0.);
    // Bouyancy forces
    G_source.is(1,vl_indx,vl_indxe).set(G_body).scale(rho).rs();
    tmp00 = double(tmp05.prod(G_body,vel,-1,-1).scale(rho));
    G_source.setel(tmp00,ndof);

  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(Ajac,normal,-1,1,2,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(Ajac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				     FastMat2 &dshapex,double w) {
  FM2STAT("ndn-0");
  tmp1.prod(Djac,dshapex,-1,2,3,4,-1,1).scale(w);
  FM2STAT("ndn-1");
  grad_N_D_grad_N.prod(tmp1,dshapex,1,2,-1,4,-1,3);
  FM2STAT("ndn-2");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,Cjac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,Cjac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}


#ifdef USE_COMP_P_SUPG
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::comp_P_supg(FastMat2 &P_supg) {

  double rho_m,tau;

    const FastMat2 &grad_N = *new_adv_dif_elemset->grad_N();

    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    P_supg.prod(Ao_grad_N,tau_supg_c,1,2,-1,-1,3);

}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void gasflow_ff::
Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
	    FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C) {
  maktgsp.make_tangent(normal);

  if (!linear_abso) {
    // Speed of sound
    double a = sqrt(ga*p/rho);
  
    tmp20.prod(vel,normal,-1,-1);
    double un = tmp20.get();
    // Riemman based b.c.'s
    double a2g = 2*a/g1;

    // Riemman Invariants
    Rie.setel(un-a2g,1);
    Rie.setel(un+a2g,2);
    double s = log(p)-ga*log(rho);
    Rie.setel(s,3);
    Rie.is(1,4,ndof)
      .prod(vel,maktgsp.tangent,-1,-1,1)
      .rs();

    // Jacobians
    drdU.set(0.);
    double agrho = a/(g1*rho);
    double agp = a/(g1*p);

    drdU.setel(+agrho,1,1);
    drdU.ir(1,1).is(2,2,ndim+1)
      .set(normal).rs();
    drdU.setel(-agp,1,ndof);

    drdU.setel(-agrho,2,1);
    drdU.ir(1,2).is(2,2,ndim+1)
      .set(normal).rs();
    drdU.setel(+agp,2,ndof);

    drdU.setel(-ga/rho,3,1);
    drdU.setel(1.0/p,3,ndof);

    drdU.is(1,4,ndof).is(2,2,ndim+1)
      .ctr(maktgsp.tangent,2,1).rs();

    // Characteristic speeds
    C.setel(un-a,1);
    C.setel(un+a,2);

    for (int k=3; k<=ndof; k++) 
      C.setel(un,k);
  } else {

    tmp20.prod(vel,normal,-1,-1);
    double un = tmp20.get();

    // Standard (linear) absorbing b.c.'s
    double rhoref = Uref.get(1);

    Uref.is(1,2,ndim+1);
    tmp20.prod(Uref,normal,-1,-1);
    Uref.rs();
    double uref = tmp20.get();
    double pref = Uref.get(ndof);
    double aref = sqrt(ga*pref/rhoref);
    double rhoaref = rhoref*aref;
    double aref2 = aref*aref;

    // These are the `equivalent' to
    // Riemman Invariants (not truly)
    Rie.setel(un-U.get(ndof)/rhoaref,1);
    Rie.setel(un+U.get(ndof)/rhoaref,2);
    Rie.setel(U.get(1)-U.get(ndof)/aref2,3);
    Rie.is(1,4,ndof)
      .prod(vel,maktgsp.tangent,-1,-1,1)
      .rs();

    // Jacobians
    drdU.set(0.);

    drdU.ir(1,1).is(2,2,ndim+1)
      .set(normal).rs();
    drdU.setel(-1.0/rhoaref,1,ndof);

    drdU.ir(1,2).is(2,2,ndim+1)
      .set(normal).rs();
    drdU.setel(+1.0/rhoaref,2,ndof);

    drdU.setel(1.0,3,1);
    drdU.setel(-1.0/aref2,3,ndof);

    drdU.is(1,4,ndof).is(2,2,ndim+1)
      .ctr(maktgsp.tangent,2,1).rs();

    // Characteristic speeds
  
    C.setel(uref-aref,1);
    C.setel(uref+aref,2);

    C.is(1,3,ndof).set(uref).rs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gasflow_ff::
compute_shock_cap_aniso(double &delta_aniso,
			FastMat2 &jvec_a) {
  delta_aniso = delta_sc_aniso;
  double vj = double(tmp_vj.prod(jvec,vel,-1,-1));
  jvec_a.set(jvec);
  FastMat2::branch();
  if (velmod>1e-10) {
    FastMat2::choose(0);
    jvec_a.axpy(vel,-vj/square(velmod));
  }
  FastMat2::leave();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void gasflow_ff::
get_C(FastMat2 &C) {
  C.set(0.);
}

extern void fastmat_prod_stat();

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void gasflow_ff::after_chunk() {
  // fastmat_prod_stat();
}
