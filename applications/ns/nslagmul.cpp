//__INSERT_LICENSE__
/* $Id: nslagmul.cpp,v 1.2 2005/03/28 16:42:53 mstorti Exp $ */

#include "./nslagmul.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "NSLagrangeMult::ask"
int NSLagrangeMult::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_mat);
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_mat_res);
   return 0;
}

void 
NSLagrangeMult::
get_comp_flags(const char *jobinfo,
	       int &comp_mat_a,int &comp_mat_res_a) {
  GET_JOBINFO_FLAG(comp_mat);
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_mat_res);
  comp_mat_a = comp_mat;
  comp_mat_res_a = comp_mat_res;
}

void
NSLagrangeMult::
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

void
NSLagrangeMult::
get_data(arg_data_list &arg_data_v,
	 arg_data *&retvalmat) {
  retvalmat = &arg_data_v[0];
}
