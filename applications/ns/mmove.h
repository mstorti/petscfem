// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove.h,v 1.1.2.1 2001/12/20 02:32:24 mstorti Exp $

#ifndef MMOVE_H
#define MMOVE_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move : public adaptor { 
public: 
  FastMat2 G, J, dNdxi, xlocp;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
  double distor_fun(FastMat2 & xlocp);
};

#endif
