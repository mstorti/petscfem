//__INSERT_LICENSE__
//$Id: arglist.cpp,v 1.12 2004/09/24 20:44:27 mstorti Exp $
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
  must_flush = 0;
  x = NULL;
  ghost_vec=NULL;
  ghost_vals=NULL;
  locst=NULL;
  retval=NULL;
  sstate=NULL;
  vector_assoc=NULL;
  profile=NULL;
  was_set=0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
arg_data::~arg_data(void) {
//    delete time_data;
//    time_data=NULL;
}
