// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove.h,v 1.1.2.2 2001/12/20 18:21:55 mstorti Exp $

#ifndef MMOVE_H
#define MMOVE_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move : public adaptor { 
public: 
  FastMat2 G, J, dNdxi, xlocp, xloc0, res_Dir;
  Matrix GGG;
  SymmetricMatrix  GG;
  DiagonalMatrix D;
  double c_volume, c_distor, distor_exp;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  double distor_fun(FastMat2 & xlocp);
};

#endif
