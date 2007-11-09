//__INSERT_LICENSE__
//$Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#include "./nsi_tet.h"
#include "./adaptor.h"
#include "./nodeload.h"

nodeload::nodeload() {
  printf("hihi en nodeload()\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void nodeload::init() {
  int ierr;

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define MAXPROPS 100
  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);
  
  //o Affect all node loads by this factor.
  TGETOPTDEF_ND(thash,double,load_fac,1.0);

  int iprop=0;
  load_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "loads",ndof);
  nprops = iprop;
  assert(nel==1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void nodeload::element_connector(const FastMat2 &xloc,
                              const FastMat2 &state_old,
                              const FastMat2 &state_new,
                              FastMat2 &res,FastMat2 &mat) {
  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));
  double *load_p = (propel.buff()+load_indx);
  
  res.set(load_p);
  if (load_fac!=1.0) res.scale(load_fac);
}
