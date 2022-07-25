// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast.h,v 1.11 2006/04/04 15:53:02 mstorti Exp $

#ifndef ELASTICITY_H
#define ELASTICITY_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class elasticity : public adaptor { 
private:
  dvector<int> elprpsindx; 
  int nprops;
  dvector<double> propel;
  int Young_modulus_indx;

  double rho,nu,Young_modulus_fac;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,res_pg,mat_pg1,mat_pg2,
    x_new,G,fG,dshapex_scaled,G_body;
  FastMat2 x_def,Jaco_def,dJaco,tmp_elast,mat_pg3,epsilon_LC;
  double detJaco_min, Jaco_pow;
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  void clean();
};

#endif
