/* $Id: laplace.cpp,v 1.1 2000/12/28 12:54:43 mstorti Exp $ */

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
 
#include "../../src/fem.h"
#include "../../src/readmesh.h"
#include "../../src/arglist.h"
#include "../../src/utils.h"
#include "../../src/getprop.h"

#include "genload.h"
#include "lapla.h"
#include <time.h>

static char help[] = "Basic finite element program.\n\n";

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

int MyKSPMonitor(KSP ksp,int n,double rnorm,void *dummy)
{
  PetscPrintf(PETSC_COMM_WORLD,
	      "iteration %d KSP Residual norm %14.12e \n",n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
    SET_ELEMSET_TYPE(lapla)
    SET_ELEMSET_TYPE(genload)
	{
	printf("not known elemset \"type\": %s\n",type);
	exit(1);
	}
}

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {
  Vec     x, res;
  Mat     A;                          /* linear system matrix */
  SLES    sles;     /* linear solver context */
  PC      pc;           /* preconditioner context */
  KSP     ksp;        /* Krylov subspace method context */
  double  norm, *sol, scal; /* norm of solution error */
  int     flg, ierr,  size,  ndim, nel, nen, neq, myrank, its;
  double tol=2e-6;
  char fcase[FLEN+1];
  Dofmap *dofmap;
  Mesh *mesh;

  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  ierr = OptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);

  // Read the mesh
  read_mesh(mesh,fcase,dofmap,neq,size,myrank);
  dofmap->create_MPI_vector(x);

  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  string save_file = string("save_file.sal");
  get_string(mesh->global_options,"save_file",save_file,1,1);
  printf("retornado por get_string: \"%s\"\n",save_file.c_str());

  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);
  ierr = VecSet(&scal,res); CHKERRA(ierr);

  ierr = opt_read_vector(mesh,x,dofmap,myrank);

  arg_list argl;

  VOID_IT(argl);
  argl.arg_add(&A,PROFILE);
  ierr = assemble(mesh,argl,dofmap,"comp_prof"); CHKERRA(ierr);

  ierr = MatZeroEntries(A); CHKERRA(ierr);
      
  // SLES solver for the laplacian
  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
  ierr = SLESSetOperators(sles,A,
			  A,DIFFERENT_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRA(ierr);

  ierr = KSPSetType(ksp,KSPCG); CHKERRA(ierr);
  ierr = PCSetType(pc,PCJACOBI); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,
         PETSC_DEFAULT); CHKERRA(ierr);
  ierr = KSPSetMonitor(ksp,MyKSPMonitor,PETSC_NULL);

#if 0
  // Computes matrix by finite differences
  VOID_IT(argl);
  argl.arg_add(&x,IN_VECTOR);
  argl.arg_add(&res,OUT_VECTOR);
  ierr = assemble(mesh,argl,dofmap,"comp_res"); CHKERRA(ierr);

  VOID_IT(argl);
  argl.arg_add(&x,PERT_VECTOR);
  argl.arg_add(&A,OUT_MATRIX_FDJ);
  ierr = assemble(mesh,argl,dofmap,"comp_res"); CHKERRA(ierr);
#else
  // Computes analytic matrix 
  VOID_IT(argl);
  argl.arg_add(&x,IN_VECTOR);
  argl.arg_add(&res,OUT_VECTOR);
  argl.arg_add(&A,OUT_MATRIX);
  ierr = assemble(mesh,argl,dofmap,"comp_res_mat"); CHKERRA(ierr);
#endif

  ierr = SLESSolve(sles,res,x,&its); CHKERRA(ierr); 
  print_vector(save_file.c_str(),x,dofmap);

  PetscFinalize();
  exit(0);
}

