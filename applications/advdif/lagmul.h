// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: lagmul.h,v 1.3 2005/01/26 20:02:22 mstorti Exp $
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

  void
  get_data(arg_data_list &arg_data_v,
	   arg_data *&stateo,
	   arg_data *&staten,
	   arg_data *&retval,
	   arg_data *&retvalmat);

  void 
  get_data(arg_data_list &arg_data_v,
	   arg_data *&retvalmat);
};

#endif
