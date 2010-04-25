// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast2.h,v 1.7 2006/04/04 15:32:43 mstorti Exp $

#ifndef PETSCFEM_RENORM2_H
#define PETSCFEM_RENORM2_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  renorm2 : public adaptor { 
private:
  FastMat2 phi, grad_phi, tmp, phipgm, tmp2, Jaco,iJaco;
  double c_reac, c_grad_phi, kond;
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
