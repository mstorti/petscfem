// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast2.h,v 1.4 2006/02/04 14:25:07 mstorti Exp $

#ifndef ELASTICITY2_H
#define ELASTICITY2_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity2 : public adaptor { 
public: 
  double rho,E,nu;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,
    res_pg,mat_pg1,mat_pg2,mass_pg,dv,a,tmp,tmp2,
    xnew,xold,xstar,vnew,vold,vstar,tmp3,tmp4;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
