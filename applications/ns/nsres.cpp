//__INSERT_LICENSE__
/* $Id: nsres.cpp,v 1.2 2001/10/16 23:00:31 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "linear_restriction::init"
void linear_restriction::init() {
  int code,node,dof,j,k,ierr,ncoef;
  double c;
  FILE *fid;
  if (!was_loaded) {
    //o The file containing coefficients
    TGETOPTDEF_S(thash,string,coef_file,"coef_file.dat");
    fid = fopen(coef_file.c_str(),"r");
    assert(fid);

    code = fscanf(fid,"%d %d\n",&nres_m,&ncoef);
    assert(code==2);
    PetscPrintf(PETSC_COMM_WORLD,"reading restriction file %s\n"
		"number of restrictions: %d\n"
		"number of non null coefs. per restriction: %d\n",
		coef_file.c_str(),nres_m,ncoef);

    coef.resize(3,nres_m,nel,ndof).set(0.);
    w.resize(3,nres_m,nel,ndof).set(0.);
    b.resize(1,nres_m).set(0.);
    node_lm.resize(nres_m);
    dofs_lm.resize(nres_m);

    for (j=0; j<nres_m; j++) {
      code = fscanf(fid,"%d %d\n",&node,&dof);
      assert(code==2);
      node_lm[j] = node;
      dofs_lm[j] = dof;
      
      for (k=0; k<ncoef; k++) {
	code = fscanf(fid,"%d %d %lf\n",&node,&dof,&c);
	assert(code==3);
	coef.setel(c,j,node,dof);
      }
      code = fscanf(fid,"%lf\n",c);
      assert(code==1);

      b.setel(c,j);
    }

    for (j=0; j<nres_m; j++) {
      for (k=0; k<ncoef; k++) {
	code = fscanf(fid,"%d %d %lf\n",&node,&dof,&c);
	assert(code==3);
	w.setel(c,j,node,dof);
      }
    }
    
    was_loaded=1;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "linear_restriction::lag_mul_dof"
void linear_restriction::lag_mul_dof(int jr,int &node,int &dof) {
  node = node_lm[jr]; 
  dof = dofs_lm[jr]; 
}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)

void linear_restriction::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & w_arg,FastMat2 & jac) {
  w_arg.set(w);
  r.prod(coef,U,1,-1,-2,-1,-2).add(b);
  jac.set(coef).scale(-1.);
  
}

