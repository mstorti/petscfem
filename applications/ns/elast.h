// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast.h,v 1.3 2001/11/30 13:22:24 mstorti Exp $

#ifndef ELASTICITY_H
#define ELASTICITY_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity : public adaptor { 
public: 
  double rho,E,nu;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,res_pg,mat_pg1,mat_pg2;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
