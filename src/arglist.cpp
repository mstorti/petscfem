//__INSERT_LICENSE__
//$Id: arglist.cpp,v 1.9 2003/07/02 02:32:47 mstorti Exp $
#include "fem.h"
#include "dofmap.h"
#include "elemset.h"
#include "arglist.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
void arg_list::arg_add(void *arg,int options,string arginfo="") {
  push_back(arg_entry(arg,options,arginfo));
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
arg_data::arg_data(void) : time_data(NULL) {
  x = NULL;
  ghost_vec=NULL;
  ghost_vals=NULL;
  locst=NULL;
  retval=NULL;
  sstate=NULL;
  vector_assoc=NULL;
  profile=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
arg_data::~arg_data(void) {
//    delete time_data;
//    time_data=NULL;
}
