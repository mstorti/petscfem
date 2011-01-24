// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elastld.h,v 1.10 2006/07/20 12:18:11 mstorti Exp $

#ifndef ELASTLD_H
#define ELASTLD_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  ld_elasticity : public adaptor { 
private:
  int use_displacement_formulation;
  double rho,E,nu,lambda,mu,cdamp;
  int ntens,nen, use_new_form;
  double newmark_gamma, newmark_beta;

  FastMat2 strain, stress, Jaco, iJaco, grad_u, F, 
    tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, 
    tmp7, tmp8, xnew, xmh, xph, vnew1, xold, G_body, 
    Id, ustar, vstar, a, res_pg, dv, dxdt,
    vold,aold, vnew, anew;

  dvector<int> elprpsindx; 
  int nprops;
  dvector<double> propel;
  int Young_modulus_indx;
  double Young_modulus_fac;
  ArgHandle state_mh_argh, state_ph_argh, res_h, mat_h;

public: 
  void init();
  void before_chunk(const char *jobinfo);
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);

  void element_connector_df(const FastMat2 &xloc,
                            const FastMat2 &state_old,
                            const FastMat2 &state_new,
                            FastMat2 &res,FastMat2 &mat);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class  ld_elasticity_load : public adaptor { 
private:
  dvector<int> elprpsindx; 
  int nprops, use_displacement_formulation;
  dvector<double> propel;
  int pressure_indx,use_new_form;
  FastMat2 nor, Jaco, tmp, xstar, state, force, tmp1;
  double defo_fac;
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
