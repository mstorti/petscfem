// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: advpenal.h,v 1.1 2005/11/06 00:35:55 mstorti Exp $
#ifndef PETSCFEM_ADVDIFPENAL_H
#define PETSCFEM_ADVDIFPENAL_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/penalize.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element 
    imposed via penalization for Advdif module. */ 
class AdvdifPenalize : public Penalize {
public:
  AdvdifPenalize(Restriction *r) : Penalize(r) { }
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class dl_penalize : public AdvdifPenalize {
public:
  // First two nodes are real nodes at both sides of the membrane. 
  // Other two nodes are lagrange multipliers. 
  dl_penalize() : AdvdifPenalize(new DLRestriction) { }
};

#endif
