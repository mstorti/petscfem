// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: lagmul.h,v 1.1 2005/01/07 20:10:44 mstorti Exp $
#ifndef PETSCFEM_ADVDIF_LAGMUL_H
#define PETSCFEM_ADVDIF_LAGMUL_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class AdvdifLagrangeMult : public GLagrangeMult {
 public:
  ASK_FUNCTION;
  virtual 
  void get_comp_flags(const char *jobinfo,
		      int &comp_mat,int &comp_mat_res);

  virtual 
  void get_comp_flags(const char *jobinfo,
		      int &comp_mat,int &comp_mat_res);
  virtual
  get_data(double *&locst,double *&retval,
	   double *&retvalmat,double *&rec_Dt);

  virtual get_data(double *&retvalmat);
};

#endif
