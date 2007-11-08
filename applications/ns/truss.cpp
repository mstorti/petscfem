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
#include "./truss.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void truss::init() {
  int ierr;

  //o Perturbation scale length for increment in computing
  // the Jacobian with finite differences. 
  TGETOPTDEF(thash,int,ndim,-1);
  assert(ndim>=0);

  //o GLobal factor affecting the `krig' option
  TGETOPTDEF(thash,double,krig_fac,1.0);

  len.resize(1,ndim);
  int nel=2;
  xlocc.resize(2,nel,ndim);
  dx.resize(2,nel,ndim);

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define MAXPROPS 100
  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);
  
  // TGETOPTNDEF_ND(thash,int,ndim,none); //nd

  int iprop=0;
  krig_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "krig",1);
  nprops = iprop;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void truss::element_connector(const FastMat2 &xloc,
                              const FastMat2 &state_old,
                              const FastMat2 &state_new,
                              FastMat2 &res,FastMat2 &mat) {
  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));
  double krig = *(propel.buff()+krig_indx);
  // printf("elem %d, krig %f\n",elem,krig);
#if 0
  if (rand()%100==0) 
    printf("krig: %f, fac %f\n",
	   krig,krig_fac);
#endif

  xlocc.set(xloc);
  dx.set(state_new);

  xlocc.ir(1,1);
  len.set(xlocc);
  xlocc.ir(1,2);
  len.rest(xlocc);
  xlocc.rs();
  double len0 = len.norm_p_all();

  dx.ir(1,1);
  len.add(dx);
  dx.ir(1,2);
  len.rest(dx);
  dx.rs();
  double len1 = len.norm_p_all();
  
  double coef = krig*(1.0-len1/len0);
  res.ir(1,1).set(len).scale(coef);
  res.ir(1,2).set(len).scale(-coef);
  res.rs();
    
}
