//__INSERT_LICENSE__
//$Id: adappg.cpp,v 1.4 2001/12/03 02:59:49 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

void adaptor_pg::init() {
  xpg.resize(1,ndim);
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  grad_state_new_pg.resize(2,ndim,ndof);
  grad_state_old_pg.resize(2,ndim,ndof);
  state_new_pg.resize(1,ndof);
  state_old_pg.resize(1,ndof);
  res_pg.resize(2,nel,ndof);
  mat_pg.resize(4,nel,ndof,nel,ndof);

  elemset_init();
}

void adaptor_pg::clean() {
  xpg.clear();
  Jaco.clear();
  dshapex.clear();  
  grad_state_new_pg.clear();
  grad_state_old_pg.clear();
  state_new_pg.clear();
  state_old_pg.clear();
  res_pg.clear();
  mat_pg.clear();

  elemset_end();
}

void adaptor_pg::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  // loop over Gauss points
  elem_init();
  for (int ipg=0; ipg<npg; ipg++) {
    
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco = Jaco.det();
    if (detJaco <= 0.) {
      printf("Jacobian of element %d is negative or null\n"
	     " Jacobian: %f\n",elem,detJaco);
      PetscFinalize();
      exit(0);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);
    iJaco.inv(Jaco);
    dshapex.prod(iJaco,dshapexi,1,-1,-1,2);

    grad_state_new_pg.prod(dshapex,state_new,1,-1,-1,2);
    grad_state_old_pg.prod(dshapex,state_old,1,-1,-1,2);

    shape.ir(2,ipg+1);
    xpg.prod(shape,xloc,-1,-1,1);
    state_new_pg.prod(shape,state_new,-1,-1,1);
    state_old_pg.prod(shape,state_old,-1,-1,1);
    shape.rs();

    pg_connector(xpg,state_old_pg,grad_state_old_pg,
		 state_new_pg,grad_state_new_pg,res_pg,mat_pg);
    mat.axpy(mat_pg,wpgdet);
    res.axpy(res_pg,wpgdet);
    
  }
  elem_end();
}
