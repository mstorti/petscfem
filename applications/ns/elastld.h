// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elastld.h,v 1.5 2006/03/19 16:04:20 mstorti Exp $

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
    xnew, vnew, xold, vold, G_body, 
    Id, mass_pg, ustar, vstar, a, res_pg, dv;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class  ld_elasticity_load : public adaptor { 
public: 
  dvector<int> elprpsindx; 
  int nprops;
  dvector<double> propel;
  int pressure_indx;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
