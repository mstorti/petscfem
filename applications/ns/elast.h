// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast.h,v 1.5 2003/09/16 21:17:29 mstorti Exp $

#ifndef ELASTICITY_H
#define ELASTICITY_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity : public adaptor { 
public: 
  double rho,E,nu;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,res_pg,mat_pg1,mat_pg2,
    x_new,G,fG,dshapex_scaled;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  void element_init();
};

#endif
