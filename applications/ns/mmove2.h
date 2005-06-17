// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: mmove2.h,v 1.1 2005/06/17 21:31:16 mstorti Exp $

#ifndef PETSCFEM_MMOVE2_H
#define PETSCFEM_MMOVE2_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  mesh_move2 : public adaptor { 
public: 
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
