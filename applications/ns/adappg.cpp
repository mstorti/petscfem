//__INSERT_LICENSE__
//$Id: adappg.cpp,v 1.6 2003/02/23 16:49:38 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2 & adaptor_pg::normal() { return normal_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adaptor_pg::init() {
  // This covers may real cases, it only remains lines in 3D.
  // But maybe it should work also. 
  assert(ndimel==ndim || ndimel==ndim-1);
  xpg.resize(1,ndim);
  normal_m.resize(1,ndim);
  g.resize(2,ndimel,ndimel);
  ig.resize(2,ndimel,ndimel);
  Jaco.resize(2,ndimel,ndim);
  dshapex.resize(2,ndim,nel);  
  grad_state_new_pg.resize(2,ndim,ndof);
  grad_state_old_pg.resize(2,ndim,ndof);
  state_new_pg.resize(1,ndof);
  state_old_pg.resize(1,ndof);
  res_pg.resize(2,nel,ndof);
  mat_pg.resize(4,nel,ndof,nel,ndof);
  tmp.resize(2,ndimel,nel);

  elemset_init();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adaptor_pg::clean() {
  // Clean up local functions
  xpg.clear();
  Jaco.clear();
  dshapex.clear();  
  grad_state_new_pg.clear();
  grad_state_old_pg.clear();
  state_new_pg.clear();
  state_old_pg.clear();
  res_pg.clear();
  mat_pg.clear();
  normal_m.clear();
  g.clear();
  ig.clear();
  tmp.clear();

  elemset_end();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void adaptor_pg::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  // loop over Gauss points
  elem_init();
  for (int ipg=0; ipg<npg; ipg++) {
    
    // Select the column of dshapexi correponding to this GP
    dshapexi.ir(3,ipg+1);
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco;
    if (ndimel==ndim) {
      detJaco = Jaco.det();
    } else {
      detJaco = Jaco.detsur(&normal_m);
      normal_m.scale(-1.); // fixme:= This is to compensate a bug in mydetsur
    }
    if (detJaco <= 0.) {
      printf("Jacobian of element %d is negative or null\n"
	     " Jacobian: %f\n",elem,detJaco);
      PetscFinalize();
      exit(0);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);
    
    g.prod(Jaco,Jaco,1,-1,2,-1);
    ig.inv(g);
    tmp.prod(ig,dshapexi,1,-1,-1,2);
    dshapex.prod(Jaco,dshapexi,-1,1,-1,2);

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
