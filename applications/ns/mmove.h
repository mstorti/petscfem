// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove.h,v 1.3 2002/11/28 17:53:03 mstorti Exp $

#ifndef MMOVE_H
#define MMOVE_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move : public adaptor { 
public: 
  FastMat2 G, J, dNdxi, xlocp, xloc0, res_Dir;
#define USE_NEWMAT
#ifdef USE_NEWMAT
  SymmetricMatrix  GG;
  DiagonalMatrix D;
#else
  FastMat2 D;
#endif
  double c_volume, c_distor, distor_exp;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  double distor_fun(FastMat2 & xlocp);
};

#endif
