//__INSERT_LICENSE__
//$Id: lapla.cpp,v 1.3 2001/04/01 01:34:47 mstorti Exp $
  
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
//#include "fastmat.h"
#include "newmat.h"
#include "../../src/getprop.h"
#include "lapla.h"

#define MAXPROP 100


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "lapla::assemble"
int lapla::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,Dofmap *dofmap,
		    const char *jobinfo,int myrank,
		    int el_start,int el_last,int iter_mode,
		    const TimeData *time_data) {

  int ierr=0;

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retvalmat,iele,j,nel,k,ndof,p,nel,q,ndof)
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
// #define IDENT(j,k) (ident[ndof*(j)+(k)]) 
// #define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  //  double *load;
  ident = dofmap->ident;
  neq = dofmap->neq;
  nnod = dofmap->nnod;
  //  load = dofmap->load;
  // fixme:= esto por ahora lo dejo asi. Que eventualmente pueda
  // ser diferente el ndof de 
  //  ndof = dofmap->ndof;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  int comp_res_mat  = !strcmp(jobinfo,"comp_res_mat");
  int comp_res  = !strcmp(jobinfo,"comp_res");
  int comp_prof  = !strcmp(jobinfo,"comp_prof");
  int set_to_pi  = !strcmp(jobinfo,"set_to_pi");
  double *locst,*retval,*retvalmat;

  if(set_to_pi) {
    locst=arg_data_v[0].locst;
  }    
  if(comp_res_mat) {
    locst=arg_data_v[0].locst;
    retval=arg_data_v[1].retval;
    retvalmat=arg_data_v[2].retval;
  } 
  if(comp_res) {
    locst=arg_data_v[0].locst;
    retval=arg_data_v[1].retval;
  } 
  if(comp_prof) {
    retvalmat=arg_data_v[0].retval;
  }

  // allocate local vecs
  int kdof, kloc, node, jdim, ipg;
  //      Matrix matloc(nen,nen), xloc(nel,ndim), veccontr(nel,ndof),
  //        locstate(nel,ndof);
  Matrix veccontr(nel,ndof),xloc(nel,ndim),
    locstate(nel,ndof);

//   if (ndof != 1) {
//     PetscPrintf(PETSC_COMM_WORLD,"ndof != 1\n"); CHKERRA(1);
//   }

  nen = nel*ndof;
  Matrix matloc(nen,nen),matlap(nel,nel),Jaco(ndim,ndim),iJaco,dshapex;
  Matrix seed;
  seed= Matrix(ndof,ndof);
  seed=0;
  for (int j=1; j<=ndof; j++) seed(j,j)=1;

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndim,nel,npg);

  int ielh=-1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;

    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
    }

    matlap=0;
    matloc = 0;
    veccontr = 0;
    if (comp_res_mat || comp_res) locstate << &(LOCST(ielh,0,0));

#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;
      Jaco = DSHAPEXI * xloc;
      //      cout << "Jaco: " << endl << Jaco << endl ;

      int elem = k;

      double detJaco = mydet(Jaco);
      if (detJaco <= 0.) {
	cout << "Jacobian of element " << k << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	exit(1);
      }
      // if (ipg==1) SHV(detJaco);
      double wpgdet = detJaco*WPG;
      iJaco = Jaco.i();
      dshapex = iJaco * DSHAPEXI;
      
      if (comp_res_mat || comp_prof || comp_res) {
	matlap += wpgdet * dshapex.t() * dshapex;
      } else if (set_to_pi) {
      } else {
	printf("Don't know how to compute jobinfo: %s\n",jobinfo);
	exit(1);
      }

    }

    matloc = kron(matlap,seed);
 
    if (comp_res_mat ) {
      // SHV(matloc);
      Matrix loc = locstate;
      reshape(loc,nen,1);
      veccontr = -matloc * loc;
      veccontr >> &(RETVAL(ielh,0,0));
      matloc >> &(RETVALMAT(ielh,0,0,0,0));
    } 

    if (comp_res) {
      Matrix loc = locstate;
      reshape(loc,nen,1);
      veccontr = -matloc * loc;
      veccontr >> &(RETVAL(ielh,0,0));
    } 


    if (comp_prof) {
      matloc >> &(RETVALMAT(ielh,0,0,0,0));
    }
    if (set_to_pi) {
      locstate << &(LOCST(ielh,0,0));
      locstate = 3.141516;
      locstate >> &(LOCST(ielh,0,0));
    } 
  }
  return 0;
}
#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
