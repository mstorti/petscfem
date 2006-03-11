// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elastld.h,v 1.1 2006/03/11 21:49:34 mstorti Exp $

#ifndef ELASTLD_H
#define ELASTLD_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  ld_elasticity : public adaptor { 
public: 
  double rho,E,nu,lambda,mu;
  int ntens,nen;

  FastMat2 strain, stress, grad_u, F, tmp1, tmp2, tmp3, 
    xnew, vnew, xold, vold, Id;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
