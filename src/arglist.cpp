//__INSERT_LICENSE__
//$Id: arglist.cpp,v 1.8 2002/03/13 02:04:59 mstorti Exp $
#include "fem.h"
#include "dofmap.h"
#include "elemset.h"
#include "arglist.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void arg_list::arg_add(void *arg,int options,string arginfo="") {
  push_back(arg_entry(arg,options,arginfo));
}

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
