// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nslagmul.h,v 1.1 2005/03/28 03:29:31 mstorti Exp $
#ifndef PETSCFEM_NSLAGMUL_H
#define PETSCFEM_NSLAGMUL_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class NSLagrangeMult : public GLagrangeMult {
 public:
  ASK_FUNCTION;
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
