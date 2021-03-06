// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast2.h,v 1.7 2006/04/04 15:32:43 mstorti Exp $

#ifndef ELASTICITY2_H
#define ELASTICITY2_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity2 : public adaptor { 
private:
  double rho,nu,cdamp;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,
    res_pg,mat_pg1,mat_pg2,mass_pg,dv,a,tmp,tmp2,
    xnew,xold,xstar,vnew,vold,vstar,tmp3,tmp4;
  FastMat2 Mass,Kstiff,xmh,xph,aold,anew,aux,aux2;

  dvector<int> elprpsindx; 
  int nprops;
  dvector<double> propel;
  int Young_modulus_indx;
  int use_displacement_formulation;
  double newmark_gamma,newmark_beta;
  ArgHandle state_mh_argh, state_ph_argh, res_h, mat_h;
  void element_connector_df(const FastMat2 &xloc,
                            const FastMat2 &state_old,
                            const FastMat2 &state_new,
                            FastMat2 &res,FastMat2 &mat);
public: 
  void init();
  void before_chunk(const char *jobinfo);
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the total, kinematic and potencial (volume) energy
    developed in the volume. */ 
class elast_energy_integrator : public gatherer {
private:
  double rho,E,nu;
  int ndim,ntens;
  FastMat2 C,e_total,e_kin,e_pot,strain,stress,tmp1;
public:
  /// perform several checks and initialization
  void init();
  /// add volume
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
