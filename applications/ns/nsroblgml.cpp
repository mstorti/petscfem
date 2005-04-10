//__INSERT_LICENSE__
// $Id: nsroblgml.cpp,v 1.1 2005/04/10 08:48:17 mstorti Exp $

#include "./nslagmul.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NSROBLagrangeMult::ask"
int NSROBLagrangeMult::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_mat);
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_mat_res);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
NSROBLagrangeMult::
get_comp_flags(const char *jobinfo,
	       int &comp_mat_a,int &comp_mat_res_a) {
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat_res);
  comp_mat_a = comp_mat;
  comp_mat_res_a = comp_mat_res;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void
NSROBLagrangeMult::
get_data(arg_data_list &arg_data_v,
	   arg_data *&stateo,
	   arg_data *&staten,
	   arg_data *&retval,
	   arg_data *&retvalmat) {
  // The trapezoidal rule integration parameter
  // arg_data *staten,*stateo,*retval,*fdj_jac,*jac_prof,*Ajac;
  int j=-1;
  staten = &arg_data_v[++j];	//[0]
  stateo = &arg_data_v[++j];	//[1]
  retval  = &arg_data_v[++j];	//[2]
  retvalmat = &arg_data_v[++j];	//[3]
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void
NSROBLagrangeMult::
get_data(arg_data_list &arg_data_v,
	 arg_data *&retvalmat) {
  retvalmat = &arg_data_v[0];
}
