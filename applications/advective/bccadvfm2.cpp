/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

extern int comp_mat_each_time_step_g,
  consistent_supg_matrix_g,
  local_time_step_g;
  
#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "advective.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int BcconvAdvFM2::ask(char *,int &)"
int BcconvAdvFM2::ask(char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_res);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "BcconvAdvFM2::assemble"
int BcconvAdvFM2::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			   Dofmap *dofmap,char *jobinfo,int myrank,
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
  FMatrix  Hloc(nel,nH),H(nH),grad_H;

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
  FMatrix veccontr(nel,ndof),xloc(nel,ndim),locstate(nel,ndof),matloc(nen,nen); 

  nen = nel*ndof;

  // Gauss Point data
  char *geom;
  thash->get_entry("geometry",geom);
  assert(geom!=NULL);
  
  GPdata gp_data(geom,ndimel,nel,npg,GP_FASTMAT2);

  // Definiciones para descargar el lazo interno
  double detJaco, wpgdet, delta_sc;

  int elem, ipg,node, jdim, kloc,lloc,ldof,ret_options=0;

  FMatrix Jaco(ndimel,ndim),flux(ndof,ndim),grad_U(ndim,ndof),
    A_grad_U(ndof),G_source(ndof),tau_supg(ndof,ndof),    
    fluxn(ndof),iJaco,normal(ndim),nor,lambda,Vr,Vr_inv,U(ndof);

  FastMat2 A_jac(3,ndim,ndof,ndof);
  FastMat2 tmp1;

#ifdef USE_FASTMAT2_CACHE
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);
#endif

  int ielh=-1;
  int start_chunk=1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc.setel(NODEDATA(node-1,jdim),kloc+1,jdim+1);
      }
      for (int ih=1; ih<=nH; ih++) {
	Hloc.setel(NODEDATA(node-1,ndim+ih-1),kloc+1,ih);
      }

    }

    locstate.set(&(LOCST(ielh,0,0)));

    matloc.set(0.);
    veccontr.set(0.);

    // DUDA: esto no se puede sacar fuera del lazo de los elementos o es lo mismo ???
#define DSHAPEXI (*gp_data.FM2_dshapexi[ipg])
#define SHAPE    (*gp_data.FM2_shape[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points

    for (ipg=0; ipg<npg; ipg++) {

      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);
      detJaco = mydetsur(Jaco,normal);
      if (detJaco <= 0.) {
	cout << "bcconv: Jacobian of element " << k << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	assert(0);
      }
      // wpgdet = detJaco*WPG;

      H.prod(SHAPE,Hloc,-1,-1,1);
      // grad_H = 0; // it is not used in this calling to flux_fun

      // state variables and gradient
      U.prod(SHAPE,locstate,-1,-1,1);
      grad_U.set(0.);  // it is not used in this calling to flux_fun

      delta_sc=0;
      double lambda_max_pg;
      ierr =  (*flux_fun)(U,ndim,iJaco,H,grad_H,flux,A_jac,
			  A_grad_U,grad_U,G_source,tau_supg,delta_sc,
			  lambda_max_pg,thash,nor,lambda,Vr,Vr_inv,
			  &(ELEMPROPS(k,0)),NULL,DEFAULT,
			  start_chunk,ret_options);

      // normal = pvec(Jaco.SubMatrix(1,1,1,ndim).t(), Jaco.SubMatrix(2,2,1,ndim).t());
      fluxn.prod(flux,normal,1,-1,-1);
      tmp1.prod(SHAPE,fluxn,1,2);
      veccontr.axpy(tmp1,WPG);

    }

    if (!weak_form) veccontr.set(0.);

    veccontr.export(&(RETVAL(ielh,0,0)));
    if (comp_mat_each_time_step_g) 
      matloc.export(&(RETVALMATT(ielh,0,0,0,0)));
      
  }
  FastMat2::void_cache();
  FastMat2::deactivate_cache();

}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
