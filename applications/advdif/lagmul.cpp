//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.4 2005/01/26 20:02:22 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#include "./lagmul.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int AdvdifLagrangeMult::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   DONT_SKIP_JOBINFO(comp_prof);
   return 0;
}

void 
AdvdifLagrangeMult::
get_comp_flags(const char *jobinfo,
	       int &comp_mat,int &comp_mat_res) {
  GET_JOBINFO_FLAG(comp_res);
  GET_JOBINFO_FLAG(comp_prof);
  comp_mat = comp_prof;
  comp_mat_res = comp_res;
}

void
AdvdifLagrangeMult::
get_data(arg_data_list &arg_data_v,
	   arg_data *&stateo,
	   arg_data *&staten,
	   arg_data *&retval,
	   arg_data *&retvalmat) {
  // The trapezoidal rule integration parameter
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

void
AdvdifLagrangeMult::
get_data(arg_data_list &arg_data_v,
	 arg_data *&retvalmat) {
  retvalmat = &arg_data_v[0];
}
