//__INSERT_LICENSE__
//$Id: blessf.cpp,v 1.1 2004/01/26 20:22:34 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./srfgath.h"

void bless_elemset0(char *type,Elemset *& elemset) {
  if (elemset) return;
  SET_ELEMSET_TYPE(field_surf_integrator)
    { }
}
