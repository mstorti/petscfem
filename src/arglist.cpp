//__INSERT_LICENSE__
//$Id: arglist.cpp,v 1.2 2001/04/01 01:35:06 mstorti Exp $
#include "arglist.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void arg_list::arg_add(void *,int,string)"
void arg_list::arg_add(void *arg,int options,string
		       arginfo="") {

  push_back(arg_entry(arg,options,arginfo));
  // vector<arg_entry>::push_back(arg_entry(arg,options,arginfo));

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "arg_data::arg_data(void)"
arg_data::arg_data(void) {
  x = NULL;
  ghost_vec=NULL;
  ghost_vals=NULL;
  locst=NULL;
  retval=NULL;
  sstate=NULL;
  vector_assoc=NULL;
}  

#undef __FUNC__
