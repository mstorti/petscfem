//__INSERT_LICENSE__
//$Id: iisdmat.cpp,v 1.33 2002/09/16 00:16:29 mstorti Exp $
// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

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
#include <petscmat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>

#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/graph.h>
#include <src/distmap2.h>
#include <src/distcont2.h>

DofPartitioner::~DofPartitioner() {}

enum PETScFEMErrors {
  iisdmat_set_value_out_of_range
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PFPETScMat::~PFPETScMat() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Adjacency graph classes
// based on STL map+set, demands too much memory, CPU time OK
#define STORE_GRAPH_1 0

// Based on dynamic vector of pair of indices with resorting,
// demands too much CPU time, RAM is OK
#define GRAPH_DV 1

// For each vertex wee keep a linked list of cells containing the
// adjacent nodes. Each insertion is O(m^2) where `m' is the average
// number of adjacent vertices. This seems to be optimal for
// FEM connectivities.
#define LINK_GRAPH 2

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PFPETScMat::PFPETScMat(int MM,const DofPartitioner &pp,MPI_Comm comm_) 
  : sles(NULL), comm(comm_), part(pp), pf_part(part), 
  lgraph1(MM,&part,comm_), 
  lgraph_dv(MM,&part,comm_), 
  lgraph_lkg(0,&part,comm_), 
  lgraph(&lgraph_lkg), 
  // lgraph(&lgraph_dv), 
  A(NULL), P(NULL), factored(0), mat_size(MM) { 
  //o Choice representation of the profile graph. Possible values are:
  // 0) Adjacency graph classes
  // based on STL map+set, demands too much memory, CPU time OK.
  // 1) Based on dynamic vector of pair of indices with resorting,
  // demands too much CPU time, RAM is OK
  // 2) For each vertex wee keep a linked list of cells containing the
  // adjacent nodes. Each insertion is $O(m^2)$ where $m$ is the average
  // number of adjacent vertices. This seems to be optimal for
  // FEM connectivities.
  TGETOPTDEF(GLOBAL_OPTIONS,int,use_compact_profile,LINK_GRAPH);
  //o Size of chunk for the dynamic vector used in computing the
  // mstrix profile. 
  TGETOPTDEF(GLOBAL_OPTIONS,int,compact_profile_graph_chunk_size,0);
  if (use_compact_profile==GRAPH_DV) {

    lgraph = &lgraph_dv;
    if (compact_profile_graph_chunk_size>0)
      lgraph_dv.set_chunk_size(compact_profile_graph_chunk_size);

  } else if (use_compact_profile==LINK_GRAPH) {
    lgraph = &lgraph_lkg;
    if (compact_profile_graph_chunk_size>0)
      lgraph_lkg.set_chunk_size(compact_profile_graph_chunk_size);
    lgraph_lkg.init(MM);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PFPETScMat::clean_prof_a() { 
  lgraph->clear(); 
  if(mat_size>0) lgraph->init(mat_size); 
  return 0; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::duplicate_a"
int PFPETScMat::duplicate_a(MatDuplicateOption op,const PFMat &A) {
  PETSCFEM_ERROR("Not implemented yet!! duplicate operation \n"
		 "for \"%s\" PFMat derived type\n",
		 typeid(*this).name());
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::build_sles"
int PFPETScMat::build_sles() {

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
  //o Uses right or left preconditioning
  TGETOPTDEF_S_PF(thash,string,preco_side,right);

  if (KSP_method == "cg" && preco_side == "right") {
    PetscPrintf(PETSC_COMM_WORLD,__FUNC__ 
		": can't choose \"right\" preconditioning with KSP CG\n");
    preco_side = "left";
  }

  ierr = SLESDestroy_maybe(sles); CHKERRQ(ierr);
  ierr = SLESCreate(comm,&sles); CHKERRQ(ierr);
  ierr = SLESSetOperators(sles,A,
			  P,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRQ(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRQ(ierr);

  set_preco(preco_type);

  // warning:= avoiding `const' restriction!!
  ierr = KSPSetType(ksp,(char *)KSP_method.c_str()); CHKERRQ(ierr);
  if (KSP_method=="gmres") {
    int (*fcn )(KSP,int);
    //o Orthogonalization method used in conjunction with GMRES. 
    // May be  #unmodified_gram_schmidt#,
    // #modified_gram_schmidt# or #ir_orthog# (default). (Iterative refinement).
    // See PETSc documentation. 
    TGETOPTDEF_S_PF(thash,string,gmres_orthogonalization,ir_orthog);

#define SETORTH(key,fun) if (gmres_orthogonalization==key) fcn = &fun
    SETORTH("ir_orthog",KSPGMRESIROrthogonalization);
    else SETORTH("unmodified_gram_schmidt",
		 KSPGMRESUnmodifiedGramSchmidtOrthogonalization);
    else SETORTH("modified_gram_schmidt",
		 KSPGMRESModifiedGramSchmidtOrthogonalization);
    else PETSCFEM_ERROR("PFPETScMat::build_sles():: "
			"Bad \"gmres_orthogonalization\": %s\n",
			gmres_orthogonalization.c_str());  
    ierr = KSPGMRESSetOrthogonalization(ksp,fcn);
    CHKERRQ(ierr);
  }
  if (preco_side == "right")
    ierr = KSPSetPreconditionerSide(ksp,PC_RIGHT);
  else if (preco_side == "left") {}
  else PetscPrintf(PETSC_COMM_WORLD,
		   "PFPETScMat::build_sles: bad \"preco_side\" option: %s\n",
		   preco_side.c_str());
    
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  ierr = KSPGMRESSetRestart(ksp,Krylov_dim); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr); 

  ierr = KSPSetMonitor(ksp,PFPETScMat_default_monitor,this,NULL);
  CHKERRQ(ierr); 
  // sles_was_built = 1; // included in `factored'
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_preco"
int PFPETScMat::set_preco(const string & preco_type) {
  // warning:= avoiding `const' restriction!!
  int ierr = PCSetType(pc,(char *)preco_type.c_str()); CHKERRQ(ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFPETScMat::clean_factor"
int PFPETScMat::clean_factor_a() {
  ierr = SLESDestroy_maybe(sles); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clean_factor_a"
int IISDMat::clean_factor_a() {
  ierr = PFPETScMat::clean_factor_a(); CHKERRQ(ierr); 
  ierr = SLESDestroy_maybe(sles_ll); CHKERRQ(ierr); 
  ierr = SLESDestroy_maybe(sles_ii); CHKERRQ(ierr); 
  ierr = MatDestroy_maybe(A_LL); CHKERRQ(ierr); 
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

  ierr = VecRestoreArray(y_loc,&a); CHKERRQ(ierr); 
  ierr = VecRestoreArray(y_loc_seq,&aa); CHKERRQ(ierr); 

  // Solve local system: x_loc_seq <- XL
  if (trans) {
    ierr = SLESSolveTranspose(sles_ll,y_loc_seq,x_loc_seq,&its_);
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

  // x_loc <- A_IL' * XI
  ierr = MatMultTranspose(A_IL,x,x_loc); CHKERRQ(ierr); 
  if (local_solver == PETSc) {
    ierr = local_solve(x_loc,x_loc,1,-1.); CHKERRQ(ierr); 
  } else {
    ierr = local_solve_SLU(x_loc,x_loc,1,-1.); CHKERRQ(ierr); 
  }
  ierr = MatMultTranspose(A_II,x,y); CHKERRQ(ierr); 
  ierr = MatMultTransposeAdd(A_LI,x_loc,y,y); CHKERRQ(ierr); 
  return 0;
}

//#define PF_CHKERRQ(ierr) assert(ierr)
#define PF_CHKERRQ(ierr) CHKERRQ(ierr)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::assembly_begin_a"
int IISDMat::assembly_begin_a(MatAssemblyType type) {
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

  ierr = MatAssemblyBegin(A_LI,type); PF_CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_II,type); PF_CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_IL,type); PF_CHKERRQ(ierr);
  if (local_solver == PETSc) {
    ierr = MatAssemblyBegin(A_LL,type); PF_CHKERRQ(ierr);
  }
  A_LL_other->scatter();
#if 0
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] t1 %f, t2 %f, t3 %f, scattered %d, sr %d\n",
			  MY_RANK,t1,t2,t3,scattered,sr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::assembly_end_a"
int IISDMat::assembly_end_a(MatAssemblyType type) {

#if 0
  // This prints the time elapsed in 
  double beg,li,ii,il,ll;
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] %d %d %d %d\n",
			  MY_RANK,N_SET[0],N_SET[1],N_SET[2],N_SET[3]);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  PetscFinalize();
  exit(0);
  beg = chrono.elapsed(); chrono.start();
  ierr = MatAssemblyEnd(A_LI,type); PF_CHKERRQ(ierr);
  li  = chrono.elapsed(); chrono.start();
  ierr = MatAssemblyEnd(A_II,type); PF_CHKERRQ(ierr);
  ii  = chrono.elapsed(); chrono.start(); 
  ierr = MatAssemblyEnd(A_IL,type); PF_CHKERRQ(ierr);
  il  = chrono.elapsed(); chrono.start();
  if (local_solver == PETSc) {
    ierr = MatAssemblyEnd(A_LL,type); PF_CHKERRQ(ierr);
  }
  ll  = chrono.elapsed(); chrono.start();
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			  "[%d] iisdmat-assembly-end beg-li-ii-il-ll-tot: "
			  "%f %f %f %f %f %f\n",
			  MY_RANK,beg,li,ii,il,ll,beg+li+ii+il+ll);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#else
  ierr = MatAssemblyEnd(A_LI,type); PF_CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_II,type); PF_CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_IL,type); PF_CHKERRQ(ierr);
  if (local_solver == PETSc) {
    ierr = MatAssemblyEnd(A_LL,type); PF_CHKERRQ(ierr);
  }
#endif
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clean_mat_a"
int IISDMat::clean_mat_a() {

  if (local_solver == PETSc) {
    ierr = MatDestroy_maybe(A_LL); PF_CHKERRQ(ierr); 
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_loc,n_loc,PETSC_NULL,
			   d_nnz_LL.begin(),&A_LL); PF_CHKERRQ(ierr); 
    ierr =  MatSetOption(A_LL, MAT_NEW_NONZERO_ALLOCATION_ERR);
    CHKERRQ(ierr); 
    ierr=MatZeroEntries(A_LL); PF_CHKERRQ(ierr);
  } else if (local_solver == SuperLU) {
    A_LL_SLU.clear().resize(n_loc,n_loc);
  } else assert(0);

  ierr=MatZeroEntries(A_IL); PF_CHKERRQ(ierr);
  ierr=MatZeroEntries(A_LI); PF_CHKERRQ(ierr);
  ierr=MatZeroEntries(A_II); PF_CHKERRQ(ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::view"
int IISDMat::view(PetscViewer viewer=PETSC_VIEWER_STDOUT_WORLD) {
  int ierr;
  PetscViewer matlab;
  if (local_solver == PETSc) {
    for (int rank=0; rank<SIZE; rank++) {
      char f[10];
      sprintf(f,"a_ll_%03d",MY_RANK);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,
			     f,&matlab); PF_CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(matlab,
			     PETSC_VIEWER_ASCII_MATLAB); PF_CHKERRQ(ierr);
      ierr = MatView(A_LL,matlab); PF_CHKERRQ(ierr);
      ierr = PetscViewerDestroy(matlab);
    }
  } else {
    A_LL_SLU.print("L-L part");
  }
  ierr = MatView(A_LI,viewer); PF_CHKERRQ(ierr);
  ierr = MatView(A_IL,viewer); PF_CHKERRQ(ierr);
  ierr = MatView(A_II,viewer); PF_CHKERRQ(ierr);

  return 0;
//    ViewerASCIIPrintf(viewer,"% IISD SLES\n");
//    ierr =  SLESView(sles,viewer);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::clean_prof_a"
int IISDMat::clean_prof_a() {
  // P is not destroyed, since P points to A

  ierr = PFPETScMat::clean_prof_a(); CHKERRQ(ierr); 
  if (local_solver == PETSc) {
    int ierr = MatDestroy_maybe(A_LL); CHKERRQ(ierr); 
  } else {
    A_LL_SLU.clear();
  }
  ierr = MatDestroy_maybe(A_LI); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_LI (loc-int)\n");

  ierr = MatDestroy_maybe(A_IL); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_IL (int-loc)\n");

  ierr = MatDestroy_maybe(A_II); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_II (int-int)\n");

  ierr = MatDestroy_maybe(A); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix shell A\n");

  delete A_LL_other;
  A_LL_other = NULL;
  return 0;
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
#define __FUNC__ "IISDMat::set_value_a"
int IISDMat::set_value_a(int row,int col,PetscScalar value,
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
      return 0;
    } 
    if (local_solver == SuperLU) {
      if (mode==ADD_VALUES) {
	val = A_LL_SLU.get(row_indx,col_indx) + value;
      } else {
	val = value;
      }
      A_LL_SLU.set(row_indx,col_indx,val);
      return 0;
    }
  } 
#ifndef DEBUG_IISD_DONT_SET_VALUES
  ierr = MatSetValues(*(AA[row_t][col_t]),
		      1,&row_indx,1,&col_indx,&value,mode);
  CHKERRQ(ierr);
  return 0;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define PETSC_OBJECT_DESTROY_MAYBE(type)	\
int type##Destroy_maybe(type &v) {		\
  if (v) {					\
    int ierr = type##Destroy(v); CHKERRQ(ierr);	\
    v = NULL;					\
  }						\
  return 0;					\
}

PETSC_OBJECT_DESTROY_MAYBE(Vec);
PETSC_OBJECT_DESTROY_MAYBE(Mat);
PETSC_OBJECT_DESTROY_MAYBE(SLES);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::maybe_factor_and_solve"
int IISDMat::maybe_factor_and_solve(Vec &res,Vec &dx,int factored=0) {

  int ierr,kloc,itss,j,jj;
  PetscViewer matlab;
  double *res_a,*res_i_a,*res_loc_a,*y_loc_seq_a,
    *x_loc_seq_a,*x_loc_a,*dx_a,scal,*x_a,*x_i_a;
  Vec res_i=NULL,x_i=NULL,res_loc=NULL,x_loc=NULL,res_loc_i=NULL;

  if (!factored) build_sles();

  if (n_int_tot > 0 ) {
    
    ierr = VecCreateMPI(comm,n_int,PETSC_DETERMINE,&res_i); PF_CHKERRQ(ierr); 
    ierr = VecDuplicate(res_i,&x_i); PF_CHKERRQ(ierr); 

    // Get diagonal part of A_II matrix for preconditioning
    ierr = VecDuplicate(res_i,&A_II_diag); PF_CHKERRQ(ierr); 
    ierr = MatGetDiagonal(A_II,A_II_diag); PF_CHKERRQ(ierr);

    if (!factored && local_solver == PETSc) {
    
      ierr = SLESDestroy_maybe(sles_ll); PF_CHKERRQ(ierr); 
      ierr = SLESCreate(PETSC_COMM_SELF,&sles_ll); PF_CHKERRQ(ierr); 
      ierr = SLESSetOperators(sles_ll,A_LL,
			      A_LL,SAME_NONZERO_PATTERN); PF_CHKERRQ(ierr); 
      ierr = SLESGetKSP(sles_ll,&ksp_ll); PF_CHKERRQ(ierr); 
      ierr = SLESGetPC(sles_ll,&pc_ll); PF_CHKERRQ(ierr); 

      ierr = KSPSetType(ksp_ll,KSPPREONLY); PF_CHKERRQ(ierr); 
      ierr = PCSetType(pc_ll,PCLU); PF_CHKERRQ(ierr); 
      ierr = PCLUSetFill(pc_ll,pc_lu_fill); PF_CHKERRQ(ierr); 

      if (use_interface_full_preco) {
	ierr = SLESDestroy_maybe(sles_ii); PF_CHKERRQ(ierr); 
	ierr = SLESCreate(comm,&sles_ii); PF_CHKERRQ(ierr); 
	ierr = SLESSetOperators(sles_ii,A_II,
				A_II,SAME_NONZERO_PATTERN); PF_CHKERRQ(ierr); 
	ierr = SLESGetKSP(sles_ii,&ksp_ii); PF_CHKERRQ(ierr); 
	ierr = SLESGetPC(sles_ii,&pc_ii); PF_CHKERRQ(ierr); 
	// ierr = KSPSetType(ksp_ii,KSPGMRES); PF_CHKERRQ(ierr); 
	ierr = KSPSetType(ksp_ii,KSPRICHARDSON); PF_CHKERRQ(ierr); 
	ierr = KSPRichardsonSetScale(ksp_ii,interface_full_preco_relax_factor);
	if(print_interface_full_preco_conv) {
	  ierr = KSPSetMonitor(ksp_ii,KSPDefaultMonitor,NULL,NULL);
	  PF_CHKERRQ(ierr); 
	}
	ierr = KSPSetTolerances(ksp_ii,0.,0.,1.e10,
				interface_full_preco_maxits); PF_CHKERRQ(ierr); 
	ierr = PCSetType(pc_ii,PCJACOBI); PF_CHKERRQ(ierr); 
      }
    }

    if (print_Schur_matrix) {
      // To print the Schur matrix by columns
      for (j = 0; j < n_int_tot; j++) {
	scal = 0.;
	ierr = VecSet(&scal,x_i); 
	PF_CHKERRQ(ierr); 
	scal = 1.;
	ierr = VecSetValues(x_i,1,&j,&scal,INSERT_VALUES);
	PF_CHKERRQ(ierr); 
      
	ierr = MatMult(A,x_i,res_i);
	PetscPrintf(comm,"For j = %d, column:\n",j);

	ierr = VecView(res_i,PETSC_VIEWER_STDOUT_SELF);
	PF_CHKERRQ(ierr); 
      }
    }

    ierr = VecCreateMPI(comm,
			n_loc,PETSC_DETERMINE,&res_loc); PF_CHKERRQ(ierr); 
    ierr = VecDuplicate(res_loc,&x_loc); PF_CHKERRQ(ierr); 

    // This could be done with a scatter
    // res -> (res_loc, res_i)
    ierr = VecGetArray(res,&res_a); PF_CHKERRQ(ierr); 
    ierr = VecGetArray(res_loc,&res_loc_a); PF_CHKERRQ(ierr); 
    ierr = VecGetArray(res_i,&res_i_a); PF_CHKERRQ(ierr); 
    
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

    ierr = VecRestoreArray(res,&res_a); PF_CHKERRQ(ierr); 
    ierr = VecRestoreArray(res_loc,&res_loc_a); PF_CHKERRQ(ierr); 
    ierr = VecRestoreArray(res_i,&res_i_a); PF_CHKERRQ(ierr); 

    // Solves system for `x_loc':
    // `x_loc   <-   - A_LL \ res_loc'
    if (local_solver == PETSc) {
      local_solve(x_loc,res_loc,0,-1.);
    } else {
      local_solve_SLU(x_loc,res_loc,0,-1.);
    }

    ierr = MatMultAdd(A_IL,x_loc,res_i,res_i);


    // Solves the interface problem (iteratively)
    ierr = SLESSolve(sles,res_i,x_i,&itss); PF_CHKERRQ(ierr); 
    
    ierr = VecDuplicate(res_loc,&res_loc_i); PF_CHKERRQ(ierr); 

    ierr = MatMult(A_LI,x_i,res_loc_i); PF_CHKERRQ(ierr);

    scal = -1.;
    ierr = VecAXPY(&scal,res_loc_i,res_loc); PF_CHKERRQ(ierr);
    
    if (local_solver == PETSc) {
      local_solve(x_loc,res_loc);
    } else {
      local_solve_SLU(x_loc,res_loc);
    }

    // Again, this could be done with a scatter
    ierr = VecGetArray(dx,&dx_a); PF_CHKERRQ(ierr); 
    ierr = VecGetArray(x_loc,&x_loc_a); PF_CHKERRQ(ierr); 
    ierr = VecGetArray(x_i,&x_i_a); PF_CHKERRQ(ierr); 
    
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

    ierr = VecRestoreArray(dx,&dx_a); PF_CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_loc,&x_loc_a); PF_CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_i,&x_i_a); PF_CHKERRQ(ierr); 

#ifdef DEBUG_IISD    
    ierr = ViewerSetFormat(matlab,
			   PETSC_VIEWER_ASCII_MATLAB,"dxiisd"); PF_CHKERRQ(ierr);
    ierr = VecView(dx,matlab);
    ierr = ViewerDestroy(matlab);
#endif
    ierr = VecDestroy_maybe(res_i); PF_CHKERRQ(ierr); 
    ierr = VecDestroy_maybe(x_i); PF_CHKERRQ(ierr); 
    ierr = VecDestroy_maybe(A_II_diag); PF_CHKERRQ(ierr); 
    ierr = VecDestroy_maybe(res_loc); PF_CHKERRQ(ierr); 
    ierr = VecDestroy_maybe(x_loc); PF_CHKERRQ(ierr); 
    ierr = VecDestroy_maybe(res_loc_i); PF_CHKERRQ(ierr); 

  } else {  // if (n_int_tot == 0 )
    
    ierr = VecGetArray(res,&res_a); PF_CHKERRQ(ierr); 

    scal=0.;
    ierr = VecSet(&scal,y_loc_seq); PF_CHKERRQ(ierr);
    for (int j = 0; j < neqp; j++) {
      int dof = dofs_proc[j];
      kloc = map_dof[dof] - n_locp;
      // y_loc_seq_a[kloc] = res_a[k1+j];
      ierr = VecSetValues(y_loc_seq,1,&kloc,&res_a[dof],
			  INSERT_VALUES); PF_CHKERRQ(ierr);
    }
    
    ierr = VecRestoreArray(res,&res_a); PF_CHKERRQ(ierr); 

    if (n_loc > 0) {

      if (local_solver == PETSc) {
	SLES sles_lll;
	KSP ksp_lll;
	PC pc_lll;

	ierr = SLESCreate(PETSC_COMM_SELF,&sles_lll); PF_CHKERRQ(ierr); 
	ierr = SLESSetOperators(sles_lll,A_LL,
				A_LL,SAME_NONZERO_PATTERN); PF_CHKERRQ(ierr); 
	ierr = SLESGetKSP(sles_lll,&ksp_lll); PF_CHKERRQ(ierr); 
	ierr = SLESGetPC(sles_lll,&pc_lll); PF_CHKERRQ(ierr); 

	ierr = KSPSetTolerances(ksp_lll,0,0,1e10,1); PF_CHKERRQ(ierr); 

	ierr = PCSetType(pc_lll,PCLU); PF_CHKERRQ(ierr); 
	ierr = KSPSetMonitor(ksp_lll,petscfem_null_monitor,PETSC_NULL,NULL);

	ierr = SLESSolve(sles_lll,y_loc_seq,x_loc_seq,&itss); PF_CHKERRQ(ierr); 

	ierr = SLESDestroy(sles_lll); CHKERRA(ierr); PF_CHKERRQ(ierr); 

      } else { // local_solver == SuperLU

	ierr = VecGetArray(y_loc_seq,&y_loc_seq_a); PF_CHKERRQ(ierr); 
	ierr = VecGetArray(x_loc_seq,&x_loc_seq_a); PF_CHKERRQ(ierr); 
	A_LL_SLU.solve(y_loc_seq_a);
	for (int j = 0; j < n_loc; j++) {
	  x_loc_seq_a[j] = y_loc_seq_a[j];
	}
	ierr = VecRestoreArray(y_loc_seq,&y_loc_seq_a); PF_CHKERRQ(ierr); 
	ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a); PF_CHKERRQ(ierr); 

      }
    }

    ierr = VecGetArray(dx,&dx_a); PF_CHKERRQ(ierr); 
    ierr = VecGetArray(x_loc_seq,&x_loc_seq_a); PF_CHKERRQ(ierr); 

    for (int j = 0; j < neqp; j++) {
      int dof = dofs_proc[j];
      kloc = map_dof[dof] - n_locp;
      dx_a[dof] = x_loc_seq_a[kloc];
    }

    ierr = VecRestoreArray(dx,&dx_a); PF_CHKERRQ(ierr); 
    ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a); PF_CHKERRQ(ierr); 
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::factor_and_solve_a"
int IISDMat::factor_and_solve_a(Vec &res,Vec &dx) {
  int ierr = maybe_factor_and_solve(res,dx,0); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::solve_only_a"
int IISDMat::solve_only_a(Vec &res,Vec &dx) {
  int ierr = maybe_factor_and_solve(res,dx,1); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int IISDMat::warn_iisdmat=0;

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
		  "PETScFEM warning: IISD solver only support any\n"
		  "preconditioning. Entered \"%s\", switching to \"PCNONE\"\n",
		  preco_type.c_str());
    }
  }
  return 0;
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
  return 0;
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::jacobi_pc_apply"
int IISDMat::jacobi_pc_apply(Vec x,Vec w) {
  int ierr;
  if (use_interface_full_preco) {
    int its;
    // Solves `w = A_II \ x' iteratively. 
    ierr = SLESSolve(sles_ii,x,w,&its);
    CHKERRQ(ierr);  
    
  } else {
    // Computes the componentwise division w = x/y. 
    ierr = VecPointwiseDivide(x,A_II_diag,w); CHKERRQ(ierr);  
  }
  return 0;
}
