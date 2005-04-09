// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nspenal.h,v 1.2 2005/04/09 17:42:03 mstorti Exp $
#ifndef PETSCFEM_NSPENAL_H
#define PETSCFEM_NSPENAL_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/penalize.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element 
    imposed via penalization for NS module. */ 
class NSPenalize : public Penalize {
public:
  NSPenalize(Restriction *r) : Penalize(r) { }
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
class dl_penalize : public NSPenalize {
public:
  // First two nodes are real nodes at both sides of the membrane. 
  // Other two nodes are lagrange multipliers. 
  dl_penalize() : NSPenalize(new DLRestriction) { }
};

#endif
