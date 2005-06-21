// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove2.h,v 1.3 2005/06/21 01:36:16 mstorti Exp $

#ifndef PETSCFEM_MMOVE2_H
#define PETSCFEM_MMOVE2_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move2 : public adaptor { 
private:
  FastMat2 u, A, x, xp, xedge, s, 
    tmp, dx, g, gp, gm;
  int nedge, nen;
  double vol_coef, qmax, qmin, chard, epsi;
  double disfun(FastMat2 &x);
  double hardfun(double q);
  void gdisfun(FastMat2 &x, FastMat2 &g);
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
