// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: adaptor.h,v 1.1.2.1 2001/10/29 14:34:41 mstorti Exp $
#ifndef ADAPTOR_H
#define ADAPTOR_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class adaptor : public ns_volume_element { 
public: 
  double rec_Dt;
  int npg,ndim,elem;
  FastMat2 shape,dshapexi,wpg,dshapex;
  GlobParam *glob_param;

  ASSEMBLE_FUNCTION;
  virtual void init()=0;
  virtual void clean() {};
  virtual void element_connector(const FastMat2 &xloc,
				 const FastMat2 &state_old,
				 const FastMat2 &state_new,
				 FastMat2 &res,FastMat2 &mat)=0;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class adaptor_f : public adaptor { 
public: 
  // Double pointers for interface with Fortran
  double *shape_p,*dshapexi_p,*wpg_p,*dshapex_p;
  // As FORTRAN can't manipulate directly FastMat2 objects
  // we pass the corresponding storage area
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat) {
    // This `(FastMat2 &)' explicit cast is for removing the
    // `const' qualifiers
    element_connector_f(((FastMat2 &)xloc).storage_begin(),
			((FastMat2 &)state_old).storage_begin(),
			((FastMat2 &)state_new).storage_begin(),
			res.storage_begin(),
			mat.storage_begin());
  }
  virtual void element_connector_f(double *xloc,double *state_old,
				   double *state_new,double *res,
				   double *mat)=0;
};

#endif
