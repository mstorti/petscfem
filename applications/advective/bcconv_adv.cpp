/*__INSERT_LICENSE__*/
//$Id: bcconv_adv.cpp,v 1.7 2001/05/30 03:58:38 mstorti Exp $

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
  
#include <vector>

#include <newmat.h>

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "advective.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "pvec"
/** Computes the vectorial product of two vectors in 3D*/
ColumnVector pvec(ColumnVector v1, ColumnVector v2) {
//  int m=Jaco.Nrows(), n=Jaco.Ncols();
  int n=v1.Nrows(),m=n-1;
  ColumnVector S(n);
  Matrix JJ(n,n);
  JJ.SubMatrix(1,1,1,n) = v1;
  JJ.SubMatrix(2,2,1,n) = v2;
  for (int k=1; k<=n; k++) {
    JJ.SubMatrix(n,n,1,n) = 0;
    JJ(n,k) = 1;
    S(k) = mydet(JJ);
  }
  return S;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int advective::ask(char *,int &)"
int BcconvAdv::ask(const char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);
   return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bcconv_adv::assemble"
int BcconvAdv::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			Dofmap *dofmap,const char *jobinfo,int myrank,
			int el_start,int el_last,int iter_mode,
			const TimeData *time_data) {

  GET_JOBINFO_FLAG(comp_res);
  //  GET_JOBINFO_FLAG(comp_mat_mass);
  //  GET_JOBINFO_FLAG(comp_diag_mat_mass);

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define RETVALMATT(iele,j,k,p,q) VEC5(retvalt,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  SGETOPTDEF(int,weak_form,1);

  int ndimel = ndim-1;
  int nen = nel*ndof;

  // Unpack Dofmap
  int neq,nnod;
  neq = dofmap->neq;
  nnod = dofmap->nnod;

  // Unpack nodedata
  int nu=nodedata->nu;
  int nH = nu-ndim;
  Matrix  Hloc(nel,nH),H(1,nH),grad_H;

  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  double *locst,*retval,*retvalt;
  //  if (comp_mat_mass || comp_diag_mat_mass) {
  //	printf("BC_CONV is not designed for matrices assembly \n");
  //	exit(1);
  //  }

  if (comp_res) {
    locst = arg_data_v[0].locst;
    retval = arg_data_v[1].retval;
    if (comp_mat_each_time_step_g) retvalt = arg_data_v[3].retval;
  }

  // allocate local vecs
  int kdof;
  Matrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),matloc(nen,nen); 

  nen = nel*ndof;

  // Gauss Point data
  //o Type of element geometry to define Gauss Point data
  TGETOPTDEF_S(thash,string,geometry,cartesian2d);
  GPdata gp_data(geometry.c_str(),ndimel,nel,npg);

  // Definiciones para descargar el lazo interno
  double detJaco, wpgdet, delta_sc;

  int elem, ipg,node, jdim, kloc,lloc,ldof,ret_options=0;

  Matrix Jaco(ndimel,ndim),flux(ndof,ndim),grad_U(ndim,ndof),
    A_grad_U(ndof,1),G_source(ndof,1),tau_supg(ndof,ndof),    
    fluxn(ndof,1),iJaco;
  ColumnVector normal(ndim);
  Matrix nor,lambda,Vr,Vr_inv;
  RowVector U(ndof);

  vector<Matrix *> A_jac;
  for (int jd=1; jd<=ndim; jd++) {
    A_jac.push_back(new Matrix(ndof,ndof));
  }

  //  DiagonalMatrix Id_ndof(ndof);
  //  Id_ndof = 1;

  int ielh=-1;
  int start_chunk=1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc(kloc+1,jdim+1) = NODEDATA(node-1,jdim);
      }
      for (int ih=1; ih<=nH; ih++) {
	Hloc(kloc+1,ih) = NODEDATA(node-1,ndim+ih-1);
      }

    }

    locstate << &(LOCST(ielh,0,0));

    matloc = 0;
    veccontr = 0;

    // DUDA: esto no se puede sacar fuera del lazo de los elementos o es lo mismo ???
#define DSHAPEXI (gp_data.dshapexi[ipg])
#define SHAPE    (gp_data.shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      Jaco = DSHAPEXI * xloc;
      detJaco = mydetsur(Jaco,normal);
      normal = -normal; // fixme:= This is to compensate a bug in mydetsur
      if (detJaco <= 0.) {
	cout << "bcconv: Jacobian of element " << k << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	assert(0);
      }
      // wpgdet = detJaco*WPG;

      H = SHAPE * Hloc;
      // grad_H = 0; // it is not used in this calling to flux_fun

      // state variables and gradient
      U = SHAPE * locstate;
      grad_U = 0;  // it is not used in this calling to flux_fun

      delta_sc=0;
      double lambda_max_pg;
      ierr =  (*flux_fun)(U,ndim,iJaco,H,grad_H,flux,A_jac,
			  A_grad_U,grad_U,G_source,tau_supg,delta_sc,
			  lambda_max_pg,thash,nor,lambda,Vr,Vr_inv,
			  &(ELEMPROPS(k,0)),NULL,DEFAULT,
			  start_chunk,ret_options);

      // normal = pvec(Jaco.SubMatrix(1,1,1,ndim).t(), Jaco.SubMatrix(2,2,1,ndim).t());
      fluxn = flux*normal;
      veccontr += WPG * SHAPE.t() * fluxn.t() ;

    }

    if (!weak_form) veccontr =0;

    veccontr >> &(RETVAL(ielh,0,0));
    if (comp_mat_each_time_step_g) 
      matloc >> &(RETVALMATT(ielh,0,0,0,0));
      
  }

  for (int jd=1; jd<=ndim; jd++) delete A_jac[jd-1];  
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
