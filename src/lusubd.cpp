//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.63 2001/12/08 12:29:56 mstorti Exp $

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
#include <src/pfmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>

enum PETScFEMErrors {
  iisdmat_set_value_out_of_range
};

PFMat * PFMat_dispatch(const char *s) {
  PFMat *A;
  IISDMat *AA;
  // IISD solver with PETSc or SuperLU local solver
  if (!strcmp(s,"iisd_superlu")) {
    AA =  new IISDMat;
    AA->local_solver = IISDMat::SuperLU;
    return AA;
  } else if (!strcmp(s,"iisd_petsc")) {
    AA =  new IISDMat;
    AA->local_solver = IISDMat::PETSc;
    return AA;
  } else if (!strcmp(s,"iisd")) {
    // local solver is chosen by default
    AA =  new IISDMat;
    return AA;
  } else if (!strcmp(s,"petsc")) {
    // PETSc (iterative) solver 
    A = new PETScMat;
    return A;
  } else if (!strcmp(s,"direct_superlu")) {
    A = new SparseDirect("SuperLU");
    return A;
  } else if (!strcmp(s,"direct") || 
	     !strcmp(s,"direct_petsc")) {
    A = new SparseDirect("PETSc");
    return A;
  } else {
    PETSCFEM_ERROR("PFMat type not known: %s\n",s);
  }
}

int PFMat::solve(Vec &res,Vec &dx) {
  int retval;
  if (!factored) {
    build_sles(&thash);
    retval = factor_and_solve(res,dx);
    factored=1;
    return retval;
  } else {
    return solve_only(res,dx);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PFMat::PFMat() : sles_was_built(0), 
  A(NULL), P(NULL), factored(0) {
#if 0
  char *s;
  int n = asprintf(&s,"PFMat matrix %p",this);
  assert(n>=0);
  string name = string(s);
  thash.register_name(name);
  free(s);
#endif
};

PFMat::~PFMat() {clear();};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::duplicate"
int PFMat::duplicate(MatDuplicateOption op,const PFMat &A) {
  PETSCFEM_ERROR("Not implemented yet!! duplicate operation \n"
		 "for \"%s\" PFMat derived type\n",
		 typeid(*this).name());
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
#undef __FUNC__
#define __FUNC__ "PETScMat::duplicate"
int PETScMat::duplicate(MatDuplicateOption op,const PFMat &B) {
  const PETScMat *BB = dynamic_cast<const PETScMat *> (&B);
  if (!BB) {
    PETSCFEM_ERROR("Not implemented yet!! duplicate operation \n"
		   "for \"%s\" PFMat derived type\n",
		   typeid(*this).name());
    return 1;
  }
  ierr = MatDuplicate(BB->A,op,&A); CHKERRA(ierr);
  P = A;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
IISDMat::~IISDMat() {
  delete part;
  part=NULL;
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
  MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  const int &neqp = dofmap->neqproc[myrank];
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
  MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  const int &neqp = dofmap->neqproc[myrank];
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
    map_dof(I->first,row_t,row_indx);
    row_indx -= n_locp;
    const Row &row = I->second;
    J1 = row.begin();
    J2 = row.end();
    for (J = J1; J != J2; J++) {
      map_dof(J->first,col_t,col_indx);
      assert (row_t == L && col_t == L);
      col_indx -= n_locp;
      if (row_indx < 0 || row_indx >= n_loc
	  || col_indx < 0 || col_indx >= n_loc) {
	printf("[%d] LL element not in this proc, \n"
	       "global/local row: %d/%d, column %d/%d,  n_loc: %d, "
	       "n_locp: %d\n",MY_RANK,I->first,row_indx,
	       J->first,col_indx,n_loc,n_locp);
	// MPI_Abort(PETSC_COMM_WORLD,iisdmat_set_value_out_of_range);
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
int IISDMat::view(Viewer viewer) {
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

  PFMat::clear();
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
#define __FUNC__ "IISDMat::map_dof"
void IISDMat::map_dof(int gdof,int &block,int &ldof) {
  ldof = map[gdof];
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

  map_dof(row,row_t,row_indx);
  map_dof(col,col_t,col_indx);
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
  if (n_int_tot > 0 ) {

#ifdef DEBUG_IISD
    if (MY_RANK==0) {
      FILE *fid = fopen("map.dat","w");
      fprintf(fid,
	      "# Created by PETSc-FEM \n# name: %s\n# type: matrix\n"
	      "# rows: %d\n# columns: %d\n","map",dofmap->neq,2);
      for (j=0; j<dofmap->neq; j++) 
	fprintf(fid,"%d %d\n",j,map[j]+1);
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
    ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			   "debug_iisd.dat",&matlab); CHKERRQ(ierr);
    ierr = ViewerSetFormat(matlab,
			   VIEWER_FORMAT_ASCII_MATLAB,"resiisd"); CHKERRQ(ierr);
    ierr = VecView(res,matlab); CHKERRQ(ierr); 
#endif
    
    ierr = VecCreateMPI(PETSC_COMM_WORLD,n_int,PETSC_DETERMINE,&res_i); CHKERRQ(ierr); 
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
	PetscPrintf(PETSC_COMM_WORLD,"For j = %d, column:\n",j);

	ierr = VecView(res_i,VIEWER_STDOUT_SELF);
	CHKERRQ(ierr); 
      }
    }
#endif

    ierr = VecCreateMPI(PETSC_COMM_WORLD,
			n_loc,PETSC_DETERMINE,&res_loc); CHKERRQ(ierr); 
    ierr = VecDuplicate(res_loc,&x_loc); CHKERRQ(ierr); 

    // This could be done with a scatter
    // res -> (res_loc, res_i)
    ierr = VecGetArray(res,&res_a); CHKERRQ(ierr); 
    ierr = VecGetArray(res_loc,&res_loc_a); CHKERRQ(ierr); 
    ierr = VecGetArray(res_i,&res_i_a); CHKERRQ(ierr); 
    
    for (j=0; j<neqp; j++) {
      jj = map[k1+j];
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
      jj = map[k1+j];
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
      kloc = map[k1+j] - n_locp;
      // y_loc_seq_a[kloc] = res_a[k1+j];
      ierr = VecSetValues(y_loc_seq,1,&kloc,&res_a[k1+j],
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
      kloc = map[k1+j] - n_locp;
      dx_a[k1+j] = x_loc_seq_a[kloc];
    }

    ierr = VecRestoreArray(dx,&dx_a); CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a); CHKERRQ(ierr); 
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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int IISDMat::warn_iisdmat=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::build_sles"
int IISDMat::build_sles(TextHashTable *thash,char *name=NULL) {
  int ierr;
  ierr = PFMat::build_sles(thash,name); CHKERRQ(ierr);
  //o Chooses the preconditioning operator. 
  TGETOPTDEF_ND(thash,double,pc_lu_fill,5.);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define DEFAULT_IISD_PC "jacobi"
#undef __FUNC__
#define __FUNC__ "PETScMat::build_sles"
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
      PetscPrintf(PETSC_COMM_WORLD,
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMatPETSc::create"
void PETScMat::create(Darray *da,const Dofmap *dofmap,
		     int debug_compute_prof=0) {
  int k,k1,k2,neqp,keq,leq,pos,sumd=0,sumdcorr=0,sumo=0,ierr,myrank;

  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);
  Node *nodep;
  // range of dofs in this processor: k1 <= dof <= k2
  k1=dofmap->startproc[myrank];
  neqp=dofmap->neqproc[myrank];
  k2=k1+neqp-1;

  // vectors for dimensioning the PETSc matrix
  int *d_nnz,*o_nnz,diag_ok;
  d_nnz = new int[neqp];
  o_nnz = new int[neqp];
  //loop over local dof's
  for (k=0;k<neqp;k++) {
    // initialize vector
    d_nnz[k]=0;
    o_nnz[k]=0;
    // global dof number
    keq=k1+k;
    if (debug_compute_prof) printf("-------- keq = %d: ",keq);
    // dofs connected to a particular dof are stored in the Darray as
    // a single linked list of integers. `pos' is a cursor to the
    // cells
    // The list for dof `keq' starts at position `keq' and is linked
    // by the `next' field. 
    // fixme:= we could have only the headers for the local dof's
    pos=keq;
    // PETSc doc says that you have to add room for the diagonal entry
    // even if it doesn't exist. But apparently it isn't needed. 
    // diag_ok=0;
    while (1) {
      // Recover cell
      nodep = (Node *)da_ref(da,pos);
      // Is the terminator cell?
      if (nodep->next==-1) break;
      // connected dof
      leq = nodep->val;
      if (debug_compute_prof) printf("%d ",leq);
      // if (leq==keq) diag_ok=1;
      if (k1<=leq && leq<=k2) {
	// Count in `diagonal' part (i.e. `leq' in the same processor
	// than `keq')
	d_nnz[k]++;
      } else {
	// Count in `off-diagonal' part (i.e. `leq' NOT in the same
	// processor than `keq')
	o_nnz[k]++;
      }	
      // Follow link
      pos = nodep->next;
    }
    if (debug_compute_prof) 
      printf("     --    d_nnz %d   o_nnz %d\n",d_nnz[k],o_nnz[k]);
    //d_nnz[k] += 1;
    //o_nnz[k] += 1;
    // To add room for the diagonal entry, even if it doesn't exist
    // if (!diag_ok) {
    //      printf("No diagonal element for dof %d in proc %d\n",keq,myrank);
    // fixme:= Do not add room for the diagonal term. 
    // d_nnz[k]+=1; 
    // }
    // For statistics
    sumd += d_nnz[k];
    sumdcorr += d_nnz[k];
    sumo += o_nnz[k];
  }

  // deallocate darray
  da_destroy(da);

  // Print statistics
  double avo,avd,avdcorr;
  avo = double(sumo)/double(neqp); // Average off-diag
  avd = double(sumd)/double(neqp); // Average diag
  avdcorr = double(sumdcorr)/double(neqp); // Average corrected
					   // (Used no more)
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "On processor %d,\n"
			  "       diagonal block terms: %d, (%f av.)\n"
			  // Corrected does not make sense anymore
			  // since it seems that PETSc does not need
			  // the diagonal terms. 
			  // "                (corrected): %d, (%f av.)\n"
			  "   off diagonal block terms: %d, (%f av.)\n",
			  // myrank,sumd,avd,sumdcorr,avdcorr,sumo,avo);
			  myrank,sumd,avd,sumo,avo);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  
  // Create matrices
  int neq=dofmap->neq;
  ierr =  MatCreateMPIAIJ(PETSC_COMM_WORLD,dofmap->neqproc[myrank],
			  dofmap->neqproc[myrank],neq,neq,
			  PETSC_NULL,d_nnz,PETSC_NULL,o_nnz,&A); 
  ierr =  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR);
  // P and A are pointers (in PETSc), otherwise this may be somewhat risky
  P=A;
  PETSCFEM_ASSERT0(ierr==0,"Error creating PETSc matrix\n");
  delete[] d_nnz;
  delete[] o_nnz;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clear"
void PFMat::clear() {
  if (sles_was_built) {
    int ierr = SLESDestroy(sles); 
    PETSCFEM_ASSERT0(ierr==0,"Error destroying SLES\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::clear"
void PETScMat::clear() {
  PFMat::clear();
  // P is not destroyed, since P points to A
  if (A) {
    int ierr = MatDestroy(A); 
    PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix \"A\"\n");
  }
#if 0    // Should us destroy P ?
  if (P) {
    int ierr = MatDestroy(P); 
    PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix \"P\"\n");
  }
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::build_sles"
int PFMat::build_sles(TextHashTable *thash,char *name=NULL) {

  int ierr;
  //o Absolute tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PFMAT(thash,double,atol,1e-6);
  //o Relative tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PFMAT(thash,double,rtol,1e-3);
  //o Divergence tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PFMAT(thash,double,dtol,1e+3);
  //o Krylov space dimension in solving the monolithic linear
  // system (Newton linear subiteration) by GMRES.
  TGETOPTDEF_ND_PFMAT(thash,int,Krylov_dim,50);
  //o Maximum iteration number in solving the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND_PFMAT(thash,int,maxits,Krylov_dim);
  //o Prints convergence in the solution of the GMRES iteration. 
  TGETOPTDEF_ND_PFMAT(thash,int,print_internal_loop_conv,0);
  //o Defines the KSP method
  TGETOPTDEF_S_ND_PFMAT(thash,string,KSP_method,gmres);
  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S_PFMAT(thash,string,preco_type,jacobi);

  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRQ(ierr);
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

  ierr = KSPSetMonitor(ksp,PFMat_default_monitor,this);
  sles_was_built = 1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::set_preco"
int PFMat::set_preco(const string & preco_type) {
  // warning:= avoiding `const' restriction!!
  int ierr = PCSetType(pc,(char *)preco_type.c_str()); CHKERRQ(ierr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::destroy_sles"
int PFMat::destroy_sles() {
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
#define __FUNC__ "PFMat::destroy_sles"
int IISDMat::destroy_sles() {
  int ierr;
  clean_factor();
  if (sles_was_built) {
    PFMat::destroy_sles();
    if (local_solver == PETSc) {
      ierr = MatDestroy(A_LL); CHKERRQ(ierr); 
      A_LL = NULL;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat_default_monitor"
int PFMat_default_monitor(KSP ksp,int n,double rnorm,void *A_) {
  PFMat *A = (PFMat *)A_;
  return A->monitor(ksp,n,rnorm);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::monitor"
int PFMat::monitor(KSP ksp,int n,double rnorm) {
  int ierr;
  if (print_internal_loop_conv) {
    if (n==0) PetscPrintf(PETSC_COMM_WORLD,
			  " Begin internal iterations "
			  "--------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD,
		"iteration %d KSP Residual_norm = %14.12e \n",n,rnorm);
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::factor_and_solve"
int PETScMat::factor_and_solve(Vec &res,Vec &dx) {
  int ierr = SLESSolve(sles,res,dx,&its_); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::solve_only"
int PETScMat::solve_only(Vec &res,Vec &dx) {
  return factor_and_solve(res,dx);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::view"
int PETScMat::view(Viewer viewer) {
  ierr = MatView(A,viewer); CHKERRQ(ierr); 
//    ierr = SLESView(sles,VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::zero_entries"
int PETScMat::zero_entries() {
  ierr=MatZeroEntries(A); CHKERRQ(ierr);
  return 0;
};

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/home/mstorti/PETSC/petscfem/tools/pfcpp")
  End: 
*/
