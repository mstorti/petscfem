//__INSERT_LICENSE__
//$Id: genload.cpp,v 1.8.10.1 2001/12/21 00:13:32 mstorti Exp $
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
//#include "fastmat.h"
#include "newmat.h"
#include <src/getprop.h>
#include "genload.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "genload::assemble"
int genload::assemble(double *retval,Nodedata *nodedata,double *locst,
		    double *locst2,Dofmap *dofmap,int ijob,
		    char *jobinfo,int myrank,
		    int el_start,int el_last,int iter_mode) {

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  int double_layer=0;
  ierr = get_int(thash,"double_layer",&double_layer,1); CHKERRA(ierr);

  // Unpack Elemset
  int npg,ndim,ndimel;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  ierr = get_int(thash,"ndimel",&ndimel); CHKERRA(ierr);
  int nen = nel*ndof;

  int nu=nodedata->nu;

  int comp_mat  = !strcmp(jobinfo,"comp_mat");
  int comp_res  = !strcmp(jobinfo,"comp_res");

  // allocate local vecs
  int kdof, kloc, node, jdim, ipg, nel2;
  Matrix veccontr(nel,ndof),xloc(nel,ndim),
    locstate(nel,ndof),vecc2;

  if (ndof != 1) {
    PetscPrintf(PETSC_COMM_WORLD,"ndof != 1\n"); CHKERRA(1);
  }

  nen = nel*ndof;
  Matrix matloc(nen,nen),Jaco(ndimel,ndim);
  ColumnVector S(ndim);

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);

  GPdata gp_data(geometry.c_str(),ndimel,nel2,npg);
  
  RowVector flux(ndof),load(ndof);
  DiagonalMatrix hfilm(ndof);
  double detJ;
  Matrix locsta1,locsta2;
  ColumnVector u1,u2;
  

  if (double_layer) {
    nel2=nel/2;
    if (2*nel2!=nel)
      PFEMERRQ("Number of nodes per element has to be even for "
	      "double_layer mode");
    // state vector in both sides of the layer
    u1.ReSize(nel2),u2.ReSize(nel2);
    ierr = get_double(thash,"hfilm",hfilm.Store(),0,ndof);
    CHKERRA(ierr);
    load=0;
    ierr = get_double(thash,"load",load.Store(),1,ndof); CHKERRA(ierr);
    locsta1.ReSize(nel2,ndof),locsta2.ReSize(nel2,ndof);
    xloc.ReSize(nel2,ndim);
    vecc2.ReSize(nel2,ndof);
  } else {
    ierr = get_double(thash,"load",load.Store(),0,ndof); CHKERRA(ierr);
    nel2=nel;
  }
  
  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel2; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
    }

    matloc = 0;
    veccontr = 0;
    locstate << &(LOCST(ielh,0,0));

    if (double_layer) {
      locsta1 = locstate.Rows(1,nel2);
      locsta2 = locstate.Rows(nel2+1,nel);
    }

#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      Jaco = DSHAPEXI * xloc;

      if (ndimel==ndim) {
	detJ = mydet(Jaco);
      } else if (ndimel==ndim-1) {
	detJ = mydetsur(Jaco,S);
	S = -S; ; // fixme:= This is to compensate a bug in mydetsur
      } else {
	PFEMERRQ("Not considered yet ndimel<ndim-1\n");
      }

      if (!double_layer) {
	flux=load;
      } else {
	u1=SHAPE*locsta1;
	u2=SHAPE*locsta2;
	flux = (hfilm*(u2-u1)).t()+load;
      }
      
      double wpgdet = detJ*WPG;
      if (double_layer) {
	vecc2 = wpgdet * SHAPE.t() * flux;
	veccontr.Rows(1,nel2) += vecc2;
	veccontr.Rows(nel2+1,nel) -= vecc2;
      } else {
	veccontr  += wpgdet * SHAPE.t() * flux;
      }

    }

    if (ijob==COMP_VEC || ijob==COMP_FDJ || ijob==COMP_FDJ_PROF)
      veccontr >> &(RETVAL(ielh,0,0));
    if (ijob==COMP_MAT || ijob==COMP_MAT_PROF)
      matloc >> &(RETVALMAT(ielh,0,0,0,0));
    
  }
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
