//__INSERT_LICENSE__
/* $Id: advpenal.cpp,v 1.1 2005/11/06 00:35:55 mstorti Exp $ */

#include "./advpenal.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "AdvdifPenalize::ask"
int AdvdifPenalize::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
AdvdifPenalize::
get_comp_flags(const char *jobinfo,
	       int &comp_mat,int &comp_mat_res) {
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);
  comp_mat = comp_prof;
  comp_mat_res = comp_res;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void
AdvdifPenalize::
get_data(arg_data_list &arg_data_v,
	   arg_data *&stateo,
	   arg_data *&staten,
	   arg_data *&retval,
	   arg_data *&retvalmat) {
  // arg_data *staten,*stateo,*retval,*fdj_jac,*jac_prof,*Ajac;
  int j=-1;
  stateo = &arg_data_v[++j]; //[0]
  staten = &arg_data_v[++j]; //[1]
  retval  = &arg_data_v[++j];//[2]
  ++j;//[3] (Not used)
  retvalmat = &arg_data_v[++j];//[4]
#if 0
  glob_param = (GlobParam *)arg_data_v[++j].user_data;;
  rec_Dt = 1./glob_param->Dt;
  if (glob_param->steady) rec_Dt = 0.;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void
AdvdifPenalize::
get_data(arg_data_list &arg_data_v,
	 arg_data *&retvalmat) {
  retvalmat = &arg_data_v[0];
}
