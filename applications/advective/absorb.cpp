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

#include <cassert>  
#include <vector>

#include <newmat.h>

using namespace std;

#include "../../src/fem.h"
#include "../../src/utils.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"
#include "advective.h"

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCSTOLD(iele,j,k) VEC3(locstold,iele,j,nel,k,ndof)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int absorb::ask(char *,int &)"
int Absorb::ask(char *jobinfo,int &skip_elemset) {

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(absorb_bc_proj);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#undef __FUNC__
#define __FUNC__ "Absorb::assemble"
int Absorb::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
		     Dofmap *dofmap,char *jobinfo,int myrank,
		     int el_start,int el_last,int iter_mode,
		     const TimeData *time_data) {

  GET_JOBINFO_FLAG(absorb_bc_proj);

  int ierr=0;
  SGETOPTDEF(int,ndim,2);
  // ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);

  char *value;
  ColumnVector x1(ndim),x2(ndim),nor(ndim);
  RowVector U1(ndof),U2(ndof),U1old(ndof),
    U1av(ndof);
  Matrix locstate(nel,ndof), locstate_old(nel,ndof), flux(ndof,ndim),
    A_jac_nor(ndof,ndof),lambda_cmplx(ndof,2),lambda(ndof,2),
    Vr(ndof,ndof),Vr_inv(ndof,ndof),Vi(ndof,ndof),tau_supg, A_grad_U, grad_U;
  double delta_sc,lambda_max;

  // These enter only for the computation of the source term, so that
  // we set them to some arbtrary value (0 for most, identity for
  // ijaco) here. 
  int nu = nodedata->nu;
  int nH = nu - ndim;
  Matrix  H(1,nH),grad_H(ndim,nH),  G_source(ndof,1);;
  DiagonalMatrix iJaco(ndim);

  vector<Matrix *> A_jac;
  for (int jd=1; jd<=ndim; jd++) 
    A_jac.push_back(new Matrix(ndof,ndof));

  double *locst, *locstold;
  if (absorb_bc_proj) {
    locst = arg_data_v[0].locst;
    locstold = arg_data_v[1].locst;
  }

  int ielh=-1;
  assert(nel==1);

  int start_chunk=1;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    ielh++;

    int node1 = ICONE(k,0);
    for (int jd=1; jd<=ndim; jd++) {
      x1(jd) = NODEDATA(node1-1,jd-1);
    }
    for (int ih=1; ih<=nH; ih++) {
      H(1,ih) = NODEDATA(node1,ndim+ih-1);
    }

    nor << ELEMPROPS(k,0) << ELEMPROPS(k,1);
    double length = sqrt(nor.SumSquare());
    nor = nor/length;

    locstate << &(LOCST(ielh,0,0));
    locstate_old << &(LOCSTOLD(ielh,0,0));

    U1 = locstate.Row(1);
    U1old = locstate_old.Row(1);
    U1av = 0.5*(U1+U1old);
    
    // Si el application writer no calcula la descomposicion
    // en autovalores, Vr queda en cero. 
    Vr=0;
    int ret_options;
    ierr =  (*flux_fun)(U1old,ndim,iJaco,H,grad_H,flux,A_jac,
			A_grad_U, grad_U,
			G_source,tau_supg,delta_sc,
			lambda_max,
			thash,nor,lambda,Vr,Vr_inv,
			&(ELEMPROPS(k,0)),NULL,COMP_EIGENV,
			start_chunk,ret_options);

    if (Vr.Norm1() <= 0.) {
      A_jac_nor = 0;
      for (int jd=1; jd<=ndim; jd++) {
	A_jac_nor += nor(jd) * AJAC(jd);
      }

      non_symm_eigenvals(A_jac_nor,lambda,Vr,Vi);

    // Check that complex part is null
      if ((lambda.Column(2)).Norm1()>0 || Vi.Norm1()>0) {
//        SHV(lambda);
//        SHV(Vr);
//        SHV(Vi);
	for (int kdof=1; kdof<ndof; kdof++) {
	  double tol=1e-6;
	  if (lambda(kdof,2)!=0 
	      && fabs(lambda(kdof,1)-lambda(kdof+1,1))<tol 
	      && fabs(lambda(kdof,2)+lambda(kdof+1,2))<tol) {
	    Vr.Column(kdof+1)=Vi.Column(kdof);
	  }
	}
//        SHV(Vr);
      //        for (int jd=1; jd<=ndim; jd++) {
      //  	SHV(nor(jd));
      //  	SHV(AJAC(jd));
      //        }
      //        SHV(U1old);
      //        SHV(A_jac_nor);
	PetscPrintf(PETSC_COMM_WORLD,
		    "System is not hyperbolic. Imaginary part of eigenvalue is not null.\n");
      //        PFEMERRQ("System is not hyperbolic. Imaginary part of eigenvalue is not null.\n");
      }
      Vr_inv = Vr.i();
    }
      
    DiagonalMatrix K(ndof);
    Matrix Pi_ingoing(ndof,ndof);
    for (int k=1; k<=ndof; k++) 
      K(k,k) = (lambda(k,1)<0 ? 1 : 0);
    Pi_ingoing = Vr * K * Vr_inv;

    U1 = U1 - (U1 - U1old) * Pi_ingoing.t();

    if (absorb_bc_proj) {
      locstate.Row(1) = U1;
      locstate >> &(LOCST(ielh,0,0));
    }
    
  }

  for (int jd=1; jd<ndim; jd++) 
    A_jac[jd-1]->~Matrix();
  
}

#undef NODEDATA
