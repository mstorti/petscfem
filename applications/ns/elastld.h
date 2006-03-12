// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elastld.h,v 1.3 2006/03/12 03:32:06 mstorti Exp $

#ifndef ELASTLD_H
#define ELASTLD_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  ld_elasticity : public adaptor { 
public: 
  double rho,E,nu,lambda,mu;
  int ntens,nen;

  FastMat2 strain, stress, Jaco, iJaco, grad_u, F, 
    tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6,
    xnew, vnew, xold, vold, 
    Id, mass_pg, ustar, vstar, a, res_pg, dv;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
