// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast2.h,v 1.7 2006/04/04 15:32:43 mstorti Exp $

#ifndef PETSCFEM_RENORM3_H
#define PETSCFEM_RENORM3_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  renorm3: public adaptor { 
private:
  double creac, kond, mpenal, total_vol, t_area, vel;
  FastMat2 phi, phiold, tmp, tmp2, tmp3, tmp4, grad_phi, 
    phipgm, phipgmold, Jaco,iJaco, C,
    phirot, xrot, xlocc, Jaco_h,xarea;
  int use_lumped;
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
