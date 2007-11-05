//__INSERT_LICENSE__
//$Id: nsvol.cpp,v 1.7 2002/01/14 03:45:05 mstorti Exp $
  
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>

#include "nsi_tet.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int ns_volume_element::ask(char *,int &)"
int ns_volume_element::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(comp_prof);
  DONT_SKIP_JOBINFO(comp_res);
  DONT_SKIP_JOBINFO(comp_mat_res);
  DONT_SKIP_JOBINFO(get_nearest_wall_element);
  return 0;
}
