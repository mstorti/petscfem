//__INSERT_LICENSE__
/* $Id: lagmul.cpp,v 1.1 2005/01/07 19:54:29 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>

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
get_data(double *&locst,double *&retval,
	 double *&retvalmat,double *&rec_Dt) {
}

void
AdvdifLagrangeMult::
virtual get_data(double *&retvalmat) {
}
