//__INSERT_LICENSE__
//$Id mstorti-v6-branch-1.0.2-10-gac227da Tue Oct 23 16:55:19 2007 -0300$

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

  len.resize(1,ndim);
  int nel=2;
  xlocc.resize(2,nel,ndim);
  dx.resize(2,nel,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void truss::element_connector(const FastMat2 &xloc,
                              const FastMat2 &state_old,
                              const FastMat2 &state_new,
                              FastMat2 &res,FastMat2 &mat) {
  xlocc.set(xloc);
  dx.set(state_new);

  double krig = 1;
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
