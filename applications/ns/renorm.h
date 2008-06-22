// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast2.h,v 1.7 2006/04/04 15:32:43 mstorti Exp $

#ifndef PETSCFEM_RENORM_H
#define PETSCFEM_RENORM_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  renorm : public adaptor { 
private:
  double creac, kond, mpenal;
  FastMat2 phi, tmp, tmp2, grad_phi, 
    phipgm, Jaco,iJaco, resh;
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
