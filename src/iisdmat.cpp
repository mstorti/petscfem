//__INSERT_LICENSE__
//$Id: iisdmat.cpp,v 1.1.2.6 2001/12/28 21:13:17 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;
int SCHED_ALG=1;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

// for inclusion of asprintf
// #define _GNU_SOURCE
#include <typeinfo>
#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <mat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
// #pragma implementation "PFMat"
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>

PFMat::PFMat() {}

PFMat::~PFMat() {}

DofPartitioner::~DofPartitioner() {}

enum PETScFEMErrors {
  iisdmat_set_value_out_of_range
};

int PFPETScMat::solve(Vec &res,Vec &dx) {
  int retval;
  if (!factored) {
    build_sles();
    retval = factor_and_solve(res,dx);
    factored=1;
    return retval;
  } else {
    return solve_only(res,dx);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clear"
void PFPETScMat::clear() {
  if (sles_was_built) {
    int ierr = SLESDestroy(sles); 
    PETSCFEM_ASSERT0(ierr==0,"Error destroying SLES\n");
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PFPETScMat::PFPETScMat(MPI_Comm comm_) : sles_was_built(0), 
  A(NULL), P(NULL), factored(0), comm(comm_)  {}
#endif

PFPETScMat::~PFPETScMat() {clear();};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::duplicate"
int PFPETScMat::duplicate(MatDuplicateOption op,const PFMat &A) {
  PETSCFEM_ERROR("Not implemented yet!! duplicate operation \n"
		 "for \"%s\" PFMat derived type\n",
		 typeid(*this).name());
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::build_sles"
int PFPETScMat::build_sles() {

#define TGETOPTDEF_ND_PF(thash,type,name,default)		\
        name = default;						\
        get_option(#name,&name); 
  
#define TGETOPTDEF_S_ND_PF(thash,type,name,default)	\
        name = string(#default);			\
        get_option(#name,name); 
  
#define TGETOPTDEF_S_PF(thash,type,name,default)	\
        string name;					\
        TGETOPTDEF_S_ND_PF(thash,type,name,default)
  
  int ierr;
  //o Absolute tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,atol,1e-6);
  //o Relative tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,rtol,1e-3);
  //o Divergence tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,double,dtol,1e+3);
  //o Krylov space dimension in solving the monolithic linear
  // system (Newton linear subiteration) by GMRES.
  TGETOPTDEF_ND_PF(thash,int,Krylov_dim,50);
  //o Maximum iteration number in solving the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PF(thash,int,maxits,Krylov_dim);
  //o Prints convergence in the solution of the GMRES iteration. 
  TGETOPTDEF_ND_PF(thash,int,print_internal_loop_conv,0);
  //o Defines the KSP method
  TGETOPTDEF_S_ND_PF(thash,string,KSP_method,gmres);
  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S_PF(thash,string,preco_type,jacobi);

  ierr = SLESCreate(comm,&sles); CHKERRQ(ierr);
  ierr = SLESSetOperators(sles,A,
			  P,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRQ(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRQ(ierr);

  set_preco(preco_type);

  // warning:= avoiding `const' restriction!!
  ierr = KSPSetType(ksp,(char *)KSP_method.c_str()); CHKERRQ(ierr);

  ierr = KSPSetPreconditionerSide(ksp,PC_RIGHT);
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  ierr = KSPGMRESSetRestart(ksp,Krylov_dim); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits);

  ierr = KSPSetMonitor(ksp,PFPETScMat_default_monitor,this);
  sles_was_built = 1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_preco"
int PFPETScMat::set_preco(const string & preco_type) {
  // warning:= avoiding `const' restriction!!
  int ierr = PCSetType(pc,(char *)preco_type.c_str()); CHKERRQ(ierr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::clean_factor"
int PFPETScMat::clean_factor() {
  if (sles_was_built) {
    int ierr = SLESDestroy(sles); CHKERRQ(ierr);
    sles=NULL;
    sles_was_built = 0;
    factored = 0;
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clean_factor"
int IISDMat::clean_factor() {
  int ierr;
  if (factored && local_solver == PETSc) {
    ierr = SLESDestroy(sles_ll); CHKERRQ(ierr); 
  }
  if (sles_was_built) {
    PFPETScMat::clean_factor();
    if (local_solver == PETSc) {
      ierr = MatDestroy(A_LL); CHKERRQ(ierr); 
      A_LL = NULL;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PPETScFMat_default_monitor"
int PFPETScMat_default_monitor(KSP ksp,int n,double rnorm,void *A_) {
  PFPETScMat *A = dynamic_cast<PFPETScMat *>((PFMat *)A_);
  assert(A);
  return A->monitor(n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::monitor"
int PFPETScMat::monitor(int n,double rnorm) {
  int ierr;
  if (print_internal_loop_conv) {
    if (n==0) PetscPrintf(comm,
			  " Begin internal iterations "
			  "--------------------------------------\n");
    PetscPrintf(comm,
		"iteration %d KSP Residual_norm = %14.12e \n",n,rnorm);
  }
  return 0;
}

const int IISDMat::D=0;
const int IISDMat::O=1;
const int IISDMat::L=0;
const int IISDMat::I=1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "petscfem_null_monitor"
int petscfem_null_monitor(KSP ksp,int n,
			  double rnorm,void *A_) {return 0;}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
IISDMat::~IISDMat() {
  delete A_LL_other;
  A_LL_other = NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::local_solve_SLU"
int IISDMat::local_solve_SLU(Vec x_loc,Vec y_loc,int trans=0,double c=1.) {

  int ierr,j;
  double *a, *aa;

  assert(!trans);
  ierr = VecGetArray(y_loc,&aa); CHKERRQ(ierr); 
  // x_loc and y_loc may be aliased, but perhaps this is somewhat dangerous
  if (x_loc != y_loc) {
    ierr = VecGetArray(x_loc,&a); CHKERRQ(ierr); 
  } else {
    a = aa;
  }
  for (j = 0; j < n_loc; j++) a[j] = c*aa[j];
  A_LL_SLU.solve(a);
  if (x_loc != y_loc) {
    ierr = VecRestoreArray(x_loc,&a); CHKERRQ(ierr); 
  }
  ierr = VecRestoreArray(y_loc,&aa); CHKERRQ(ierr); 
  return 0;
}
  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::local_solve"
int IISDMat::local_solve(Vec x_loc,Vec y_loc,int trans=0,double c=1.) {
  int ierr,j;
  double *a,*aa;
  // ------------- IT SEEMS THIS COMMENTS ARE WRONG ????
  // For an MPI vector on the local part, solves the
  // `A_LL * x_loc = y_loc' problem.
  // `x_loc' and `y_loc' may be aliased. 
  // Localize on procesor and make y_loc_seq <- - A_LI * XI
  // -----------------------------------------
  // Solves for the local part, from MPI vectors that contain only the
  // local part. Localizes the vectors, copy on Petsc seq, vectors and
  // solve. 

  ierr = VecGetArray(y_loc,&a); CHKERRQ(ierr); 

  ierr = VecGetArray(y_loc_seq,&aa); CHKERRQ(ierr); 

  for (j = 0; j < n_loc; j++) aa[j] = c*a[j];

  ierr = VecRestoreArray(y_loc,&a);
  ierr = VecRestoreArray(y_loc_seq,&aa);

  // Solve local system: x_loc_seq <- XL
  if (trans) {
    ierr = SLESSolveTrans(sles_ll,y_loc_seq,x_loc_seq,&its_);
    CHKERRQ(ierr); 
  } else {
    ierr = SLESSolve(sles_ll,y_loc_seq,x_loc_seq,&its_);
    CHKERRQ(ierr); 
  }
  
  // Pass to global vector: x_loc <- XL
  ierr = VecGetArray(x_loc_seq,&aa); CHKERRQ(ierr); 
  ierr = VecGetArray(x_loc,&a); CHKERRQ(ierr); 

  for (j = 0; j < n_loc; j++) a[j] = aa[j];

  ierr = VecRestoreArray(x_loc_seq,&aa); CHKERRQ(ierr); 
  ierr = VecRestoreArray(x_loc,&a); CHKERRQ(ierr); 

  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::mult"
int IISDMat::mult(Vec x,Vec y) {

  int myrank;
  MPI_Comm_rank(comm, &myrank);
  // const int &neqp = dofmap->neqproc[myrank];
  // We are solving the kernel of the linear system
  //
  // ALL XL + ALI XI = 0
  // AIL XL + AII XI = RI
  //
  // XI comes in x (interface nodes) and we have to compute y <- RI 
  
  int j,ierr,its_;
  double *a,*aa;;

  // x_loc <- A_LI * XI
  ierr = MatMult(A_LI,x,x_loc); CHKERRQ(ierr); 
  if (local_solver == PETSc) {
    ierr = local_solve(x_loc,x_loc,0,-1.); CHKERRQ(ierr); 
  } else {
    ierr = local_solve_SLU(x_loc,x_loc,0,-1.); CHKERRQ(ierr); 
  }
    
  ierr = MatMult(A_II,x,y); CHKERRQ(ierr); 
  ierr = MatMultAdd(A_IL,x_loc,y,y); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::mult_trans"
int IISDMat::mult_trans(Vec x,Vec y) {

  int myrank;
  MPI_Comm_rank(comm, &myrank);
  // const int &neqp = dofmap->neqproc[myrank];
  // We are solving the kernel of the linear system
  //
  // ALL' XL + AIL' XI = 0
  // ALI' XL + AII' XI = RI
  //
  // XI comes in x (interface nodes) and we have to compute y = RI 
  
  int j,ierr,its_;
  double *a,*aa;;

  // x_loc <- A_IL' * XI
  ierr = MatMultTrans(A_IL,x,x_loc); CHKERRQ(ierr); 
  if (local_solver == PETSc) {
    ierr = local_solve(x_loc,x_loc,1,-1.); CHKERRQ(ierr); 
  } else {
    ierr = local_solve_SLU(x_loc,x_loc,1,-1.); CHKERRQ(ierr); 
  }
  ierr = MatMultTrans(A_II,x,y); CHKERRQ(ierr); 
  ierr = MatMultTransAdd(A_LI,x_loc,y,y); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::assembly_begin"
int IISDMat::assembly_begin(MatAssemblyType type) {
  int ierr;

  DistMat::const_iterator I,I1,I2;
  Row::const_iterator J,J1,J2;
  int row_indx,col_indx,row_t,col_t;
  double v,val;

  A_LL_other->scatter();

  I1 = A_LL_other->begin();
  I2 = A_LL_other->end();
  for (I = I1; I != I2; I++) {
    map_dof_fun(I->first,row_t,row_indx);
    row_indx -= n_locp;
    const Row &row = I->second;
    J1 = row.begin();
    J2 = row.end();
    for (J = J1; J != J2; J++) {
      map_dof_fun(J->first,col_t,col_indx);
      assert (row_t == L && col_t == L);
      col_indx -= n_locp;
      if (row_indx < 0 || row_indx >= n_loc
	  || col_indx < 0 || col_indx >= n_loc) {
	printf("[%d] LL element not in this proc, \n"
	       "global/local row: %d/%d, column %d/%d,  n_loc: %d, "
	       "n_locp: %d\n",MY_RANK,I->first,row_indx,
	       J->first,col_indx,n_loc,n_locp);
	// MPI_Abort(comm,iisdmat_set_value_out_of_range);
      } 
      // printf("[%d] debuffering (%d,%d) -> %f\n",
      // MY_RANK,I->first,J->first,J->second);
      v = J->second;
#ifndef DEBUG_IISD_DONT_SET_VALUES
      if (local_solver == PETSc ) {
	MatSetValues(A_LL,1,&row_indx,1,&col_indx,&v,insert_mode);
      } else {
	if (insert_mode==ADD_VALUES) {
	  val = A_LL_SLU.get(row_indx,col_indx) + v;
	} else {
	  val = v;
	}
	A_LL_SLU.set(row_indx,col_indx,val);
      }	
#endif
    }
  }
  A_LL_other->clear();

  ierr = MatAssemblyBegin(A_LI,type); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_II,type); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_IL,type); CHKERRQ(ierr);
  if (local_solver == PETSc) {
    ierr = MatAssemblyBegin(A_LL,type); CHKERRQ(ierr);
  }
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::assembly_end"
int IISDMat::assembly_end(MatAssemblyType type) {
  int ierr;
  ierr = MatAssemblyEnd(A_LI,type); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_II,type); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_IL,type); CHKERRQ(ierr);
  if (local_solver == PETSc) {
    ierr = MatAssemblyEnd(A_LL,type); CHKERRQ(ierr);
  }

  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::zero_entries"
int IISDMat::zero_entries() {
  int ierr;

  if (local_solver == PETSc) {
    if (A_LL) {
      ierr = MatDestroy(A_LL); CHKERRQ(ierr); 
      A_LL = NULL;
    }    
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_loc,n_loc,PETSC_NULL,
			   d_nnz_LL.begin(),&A_LL); CHKERRQ(ierr); 
    ierr=MatZeroEntries(A_LL); CHKERRQ(ierr);
  } else {
    A_LL_SLU.clear().resize(n_loc,n_loc);
  }
  ierr=MatZeroEntries(A_IL); CHKERRQ(ierr);
  ierr=MatZeroEntries(A_LI); CHKERRQ(ierr);
  ierr=MatZeroEntries(A_II); CHKERRQ(ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::view"
int IISDMat::view(Viewer viewer=VIEWER_STDOUT_WORLD) {
  int ierr;
  Viewer matlab;
  if (local_solver == PETSc) {
    for (int rank=0; rank<SIZE; rank++) {
      char f[10];
      sprintf(f,"a_ll_%03d",MY_RANK);
      ierr = ViewerASCIIOpen(PETSC_COMM_SELF,
			     f,&matlab); CHKERRQ(ierr);
      ierr = ViewerSetFormat(matlab,
			     VIEWER_FORMAT_ASCII_MATLAB,f); CHKERRQ(ierr);
      ierr = MatView(A_LL,matlab); CHKERRQ(ierr);
      ierr = ViewerDestroy(matlab);
    }
  } else {
    A_LL_SLU.print("L-L part");
  }
  ierr = MatView(A_LI,viewer); CHKERRQ(ierr);
  ierr = MatView(A_IL,viewer); CHKERRQ(ierr);
  ierr = MatView(A_II,viewer); CHKERRQ(ierr);

  return 0;
//    ViewerASCIIPrintf(viewer,"% IISD SLES\n");
//    ierr =  SLESView(sles,viewer);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clear"
void IISDMat::clear() {
  // P is not destroyed, since P points to A
  int ierr;

  PFPETScMat::clear();
  if (local_solver == PETSc) {
    if (A_LL) {
      int ierr = MatDestroy(A_LL); 
      PETSCFEM_ASSERT0(ierr==0,
		       "Error destroying PETSc matrix A_LL (loc-loc)\n");
      A_LL=NULL;
    }
  } else {
    A_LL_SLU.clear();
  }
  ierr = MatDestroy(A_LI); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_LI (loc-int)\n");
  ierr = MatDestroy(A_IL); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_IL (int-loc)\n");
  ierr = MatDestroy(A_II); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_II (int-int)\n");

  ierr = MatDestroy(A); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix shell A\n");
  delete A_LL_other;
  A_LL_other = NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::map_dof_fun"
void IISDMat::map_dof_fun(int gdof,int &block,int &ldof) {
  ldof = map_dof[gdof];
  if (ldof < n_loc_tot) {
    block = L;
  } else {
    block = I;
    ldof -= n_loc_tot;
  }
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::set_value"
void IISDMat::set_value(int row,int col,Scalar value,
			InsertMode mode=ADD_VALUES) {
  int row_indx,col_indx,row_t,col_t;
  double val;

  map_dof_fun(row,row_t,row_indx);
  map_dof_fun(col,col_t,col_indx);
  insert_mode = mode;  
  if (row_t == L && col_t == L) {
    row_indx -= n_locp;
    col_indx -= n_locp;
    if (row_indx < 0 || row_indx >= n_loc
	|| col_indx < 0 || col_indx >= n_loc) {
      // printf("[%d] buffering (%d,%d) -> %f\n",MY_RANK,row,col,value);
      A_LL_other->insert_val(row,col,value);
      return;
    } 
    if (local_solver == SuperLU) {
      if (mode==ADD_VALUES) {
	val = A_LL_SLU.get(row_indx,col_indx) + value;
      } else {
	val = value;
      }
      A_LL_SLU.set(row_indx,col_indx,val);
      return;
    }
  } 
#ifndef DEBUG_IISD_DONT_SET_VALUES
  MatSetValues(*(AA[row_t][col_t]),
	       1,&row_indx,1,&col_indx,&value,mode);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::solve"
int IISDMat::solve_only(Vec &res,Vec &dx) {
  return factor_and_solve(res,dx);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::factor_and_solve"
int IISDMat::factor_and_solve(Vec &res,Vec &dx) {

  int ierr,kloc,itss,j,jj;
  Viewer matlab;
  double *res_a,*res_i_a,*res_loc_a,*y_loc_seq_a,
    *x_loc_seq_a,*x_loc_a,*dx_a,scal,*x_a,*x_i_a;
  Vec res_i,x_i,res_loc,x_loc,res_loc_i;

  TGETOPTDEF_ND_PF(thash,double,pc_lu_fill,5.);

  if (n_int_tot > 0 ) {

#ifdef DEBUG_IISD
    if (MY_RANK==0) {
      FILE *fid = fopen("map.dat","w");
      fprintf(fid,
	      "# Created by PETSc-FEM \n# name: %s\n# type: matrix\n"
	      "# rows: %d\n# columns: %d\n","map",dofmap->neq,2);
      for (j=0; j<dofmap->neq; j++) 
	fprintf(fid,"%d %d\n",j,map_dof[j]+1);
      fprintf(fid,
	      "# Created by PETSc-FEM \n# name: %s\n# type: matrix\n"
	      "# rows: %d\n# columns: %d\n","n_int_v",SIZE+1,1);
      for (j=0; j<SIZE+1; j++) 
	fprintf(fid,"%d\n",n_int_v[j]);
      fprintf(fid,
	      "# Created by PETSc-FEM \n# name: %s\n# type: matrix\n"
	      "# rows: %d\n# columns: %d\n","n_loc_v",SIZE+1,1);
      for (j=0; j<SIZE+1; j++) 
	fprintf(fid,"%d\n",n_loc_v[j]);
      fclose(fid);
    }
    ierr = ViewerASCIIOpen(comm,
			   "debug_iisd.dat",&matlab); CHKERRQ(ierr);
    ierr = ViewerSetFormat(matlab,
			   VIEWER_FORMAT_ASCII_MATLAB,"resiisd"); CHKERRQ(ierr);
    ierr = VecView(res,matlab); CHKERRQ(ierr); 
#endif
    
    ierr = VecCreateMPI(comm,n_int,PETSC_DETERMINE,&res_i); CHKERRQ(ierr); 
    ierr = VecDuplicate(res_i,&x_i); CHKERRQ(ierr); 

    // Get diagonal part of A_II matrix for preconditioning
    ierr = VecDuplicate(res_i,&A_II_diag); CHKERRQ(ierr); 
    ierr = MatGetDiagonal(A_II,A_II_diag); CHKERRQ(ierr);

    if (!factored && local_solver == PETSc) {
    
      ierr = SLESCreate(PETSC_COMM_SELF,&sles_ll); CHKERRQ(ierr); 
      ierr = SLESSetOperators(sles_ll,A_LL,
			      A_LL,SAME_NONZERO_PATTERN); CHKERRQ(ierr); 
      ierr = SLESGetKSP(sles_ll,&ksp_ll); CHKERRQ(ierr); 
      ierr = SLESGetPC(sles_ll,&pc_ll); CHKERRQ(ierr); 

      ierr = KSPSetType(ksp_ll,KSPPREONLY); CHKERRQ(ierr); 
      ierr = PCSetType(pc_ll,PCLU); CHKERRQ(ierr); 
      // printf("setting pc_lu_fill = %f\n",pc_lu_fill);
      // ierr = PCLUSetFill(pc_ll,pc_lu_fill); CHKERRQ(ierr); 
      ierr = PCLUSetUseInPlace(pc_ll); CHKERRQ(ierr);

    }
#if 0 // To print the Schur matrix by columns
    for (int kk=1; kk<=2; kk++) {
      for (j = 0; j < n_int_tot; j++) {
	scal = 0.;
	ierr = VecSet(&scal,x_i); 
	CHKERRQ(ierr); 
	scal = 1.;
	ierr = VecSetValues(x_i,1,&j,&scal,INSERT_VALUES);
	CHKERRQ(ierr); 
      
	ierr = MatMult(A,x_i,res_i);
	PetscPrintf(comm,"For j = %d, column:\n",j);

	ierr = VecView(res_i,VIEWER_STDOUT_SELF);
	CHKERRQ(ierr); 
      }
    }
#endif

    ierr = VecCreateMPI(comm,
			n_loc,PETSC_DETERMINE,&res_loc); CHKERRQ(ierr); 
    ierr = VecDuplicate(res_loc,&x_loc); CHKERRQ(ierr); 

    // This could be done with a scatter
    // res -> (res_loc, res_i)
    ierr = VecGetArray(res,&res_a); CHKERRQ(ierr); 
    ierr = VecGetArray(res_loc,&res_loc_a); CHKERRQ(ierr); 
    ierr = VecGetArray(res_i,&res_i_a); CHKERRQ(ierr); 
    
    for (j=0; j<neqp; j++) {
      int dof = dofs_proc[j];
      jj = map_dof[dof];
      if (jj < n_loc_tot) {
	jj -= n_locp;
	res_loc_a[jj] = res_a[j];
      } else {
	jj -= n_intp;
	res_i_a[jj] = res_a[j];
      }
    }

    ierr = VecRestoreArray(res,&res_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(res_loc,&res_loc_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(res_i,&res_i_a); CHKERRQ(ierr); 

    // Solves system for `x_loc':
    // `x_loc   <-   - A_LL \ res_loc'
    if (local_solver == PETSc) {
      local_solve(x_loc,res_loc,0,-1.);
    } else {
      local_solve_SLU(x_loc,res_loc,0,-1.);
    }
    ierr = MatMultAdd(A_IL,x_loc,res_i,res_i);

    // Solves the interface problem (iteratively)
    ierr = SLESSolve(sles,res_i,x_i,&itss); CHKERRQ(ierr); 
    
    ierr = VecDuplicate(res_loc,&res_loc_i); CHKERRQ(ierr); 

    ierr = MatMult(A_LI,x_i,res_loc_i); CHKERRQ(ierr);

    scal = -1.;
    ierr = VecAXPY(&scal,res_loc_i,res_loc); CHKERRQ(ierr);
    
    if (local_solver == PETSc) {
      local_solve(x_loc,res_loc);
    } else {
      local_solve_SLU(x_loc,res_loc);
    }

    // Again, this could be done with a scatter
    ierr = VecGetArray(dx,&dx_a); CHKERRQ(ierr); 
    ierr = VecGetArray(x_loc,&x_loc_a); CHKERRQ(ierr); 
    ierr = VecGetArray(x_i,&x_i_a); CHKERRQ(ierr); 
    
    for (j=0; j<neqp; j++) {
      int dof = dofs_proc[j];
      jj = map_dof[dof];
      if (jj < n_loc_tot) {
	jj -= n_locp;
	dx_a[j] = x_loc_a[jj];
      } else {
	jj -= n_intp;
	dx_a[j] = x_i_a[jj];
      }
    }

    ierr = VecRestoreArray(dx,&dx_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_loc,&x_loc_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_i,&x_i_a); CHKERRQ(ierr); 

#ifdef DEBUG_IISD    
    ierr = ViewerSetFormat(matlab,
			   VIEWER_FORMAT_ASCII_MATLAB,"dxiisd"); CHKERRQ(ierr);
    ierr = VecView(dx,matlab);
    ierr = ViewerDestroy(matlab);
#endif

    ierr = VecDestroy(res_i); CHKERRQ(ierr); 
    ierr = VecDestroy(x_i); CHKERRQ(ierr); 
    ierr = VecDestroy(A_II_diag); CHKERRQ(ierr); 
    ierr = VecDestroy(res_loc); CHKERRQ(ierr); 
    ierr = VecDestroy(x_loc); CHKERRQ(ierr); 
    ierr = VecDestroy(res_loc_i); CHKERRQ(ierr); 

  } else {  // if (n_int_tot == 0 )
    
    assert(!factored);
    ierr = VecGetArray(res,&res_a); CHKERRQ(ierr); 

    scal=0.;
    ierr = VecSet(&scal,y_loc_seq); CHKERRQ(ierr);
    for (int j = 0; j < neqp; j++) {
      int dof = dofs_proc[j];
      kloc = map_dof[dof] - n_locp;
      // y_loc_seq_a[kloc] = res_a[k1+j];
      ierr = VecSetValues(y_loc_seq,1,&kloc,&res_a[dof],
			  INSERT_VALUES); CHKERRQ(ierr);
    }
    
    ierr = VecRestoreArray(res,&res_a); CHKERRQ(ierr); 

    if (n_loc > 0) {

      if (local_solver == PETSc) {
	SLES sles_lll;
	KSP ksp_lll;
	PC pc_lll;

	ierr = SLESCreate(PETSC_COMM_SELF,&sles_lll); CHKERRQ(ierr); 
	ierr = SLESSetOperators(sles_lll,A_LL,
				A_LL,SAME_NONZERO_PATTERN); CHKERRQ(ierr); 
	ierr = SLESGetKSP(sles_lll,&ksp_lll); CHKERRQ(ierr); 
	ierr = SLESGetPC(sles_lll,&pc_lll); CHKERRQ(ierr); 

	ierr = KSPSetType(ksp_lll,KSPGMRES); CHKERRQ(ierr); 

	ierr = KSPSetTolerances(ksp_lll,0,0,1e10,1); CHKERRQ(ierr); 

	ierr = PCSetType(pc_lll,PCLU); CHKERRQ(ierr); 
	ierr = KSPSetMonitor(ksp_lll,petscfem_null_monitor,PETSC_NULL);

	ierr = SLESSolve(sles_lll,y_loc_seq,x_loc_seq,&itss); CHKERRQ(ierr); 

	ierr = SLESDestroy(sles_lll); CHKERRA(ierr); CHKERRQ(ierr); 

      } else { // local_solver == SuperLU

	ierr = VecGetArray(y_loc_seq,&y_loc_seq_a); CHKERRQ(ierr); 
	ierr = VecGetArray(x_loc_seq,&x_loc_seq_a); CHKERRQ(ierr); 
	A_LL_SLU.solve(y_loc_seq_a);
	for (int j = 0; j < n_loc; j++) {
	  x_loc_seq_a[j] = y_loc_seq_a[j];
	}
	ierr = VecRestoreArray(y_loc_seq,&y_loc_seq_a); CHKERRQ(ierr); 
	ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a); CHKERRQ(ierr); 

      }
    }

    ierr = VecGetArray(dx,&dx_a); CHKERRQ(ierr); 
    ierr = VecGetArray(x_loc_seq,&x_loc_seq_a); CHKERRQ(ierr); 

    for (int j = 0; j < neqp; j++) {
      int dof = dofs_proc[j];
      kloc = map_dof[dof] - n_locp;
      dx_a[dof] = x_loc_seq_a[kloc];
    }

    ierr = VecRestoreArray(dx,&dx_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a); CHKERRQ(ierr); 
  }
  return 0;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clean_factor"
int IISDMat::clean_factor() {
  int ierr;
  if (factored && local_solver == PETSc) {
    ierr = SLESDestroy(sles_ll); CHKERRQ(ierr); 
  }
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int IISDMat::warn_iisdmat=0;

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::build_sles"
int IISDMat::build_sles() {
  int ierr;
  //o Chooses the preconditioning operator. 
  TGETOPTDEF_ND(&thash,double,pc_lu_fill,5.);
  return 0;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define DEFAULT_IISD_PC "jacobi"
#undef __FUNC__
#define __FUNC__ "IISDMat::set_preco"
int IISDMat::set_preco(const string & preco_type) {
  int ierr;
  if (preco_type=="jacobi" || preco_type=="") {
    ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
    ierr = PCShellSetApply(pc,&iisd_jacobi_pc_apply,this); 
    // printf("[%d] setting apply to %p\n",MY_RANK,&iisd_jacobi_pc_apply);
  } else if (preco_type=="none" ) {
    ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  } else {
    if ( !warn_iisdmat ) {
      warn_iisdmat=1;
      PetscPrintf(comm,
		  "PETScFEM warning: IISD operator does not support any\n"
		  "preconditioning. Entered \"%s\", switching to \"PCNONE\"\n",
		  preco_type.c_str());
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "iisd_jacobi_pc_apply"
int iisd_jacobi_pc_apply(void *ctx,Vec x ,Vec y) {
  int ierr;
  PFMat *A = (PFMat *) ctx;
  IISDMat *AA;
  AA = dynamic_cast<IISDMat *> (A);
  ierr = (AA==NULL); CHKERRQ(ierr);
  AA->jacobi_pc_apply(x,y);
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::jacobi_pc_apply"
int IISDMat::jacobi_pc_apply(Vec x,Vec w) {
  int ierr;
  // Computes the componentwise division w = x/y. 
  ierr = VecPointwiseDivide(x,A_II_diag,w); CHKERRQ(ierr);  
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::jacobi_pc_apply"
void IISDMat::print(void) {
  lgraph.print();
}
#endif

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/home/mstorti/PETSC/petscfem/tools/pfcpp")
  End: 
*/
