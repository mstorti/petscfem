//__INSERT_LICENSE__
//$Id: petscmat.cpp,v 1.1.2.4 2002/01/09 20:33:09 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;
#include <typeinfo>
#include <mat.h>

#include <src/libretto.h>

#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/petscmat.h>
#include <src/graph.h>

PETScMat::~PETScMat() {clear(*gu);};

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
  ierr = MatDuplicate(BB->A,op,&A); CHKERRQ(ierr);
  P = A;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::create_a"
int PETScMat::create_a() {
  int k,neqp,keq,leq,pos,sumd=0,sumdcorr=0,sumo=0,ierr,myrank,
    debug_compute_prof=0;
  set<int> ngbrs_v;
  set<int>::iterator q,qe;

  MPI_Comm_rank (comm, &myrank);

  // Scatter the profile graph
  lgraph.scatter();

  const int &neq = M;
  
  // number of dofs in this processor
  neqp = 0;
  dofs_proc.clear();
  proc2glob.clear();
  for (k=0; k<neq; k++) {
    if (part.processor(k)==myrank) {
      dofs_proc.push_back(k);
      proc2glob[k] = neqp++;
    }
  }
  dofs_proc_v = dofs_proc.begin();

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
    keq = dofs_proc_v[k];
    if (debug_compute_prof) printf("-------- keq = %d: ",keq);
    ngbrs_v.clear();
    lgraph.set_ngbrs(keq,ngbrs_v);
    // PETSc doc says that you have to add room for the diagonal entry
    // even if it doesn't exist. But apparently it isn't needed. 
    // diag_ok=0;
    qe = ngbrs_v.end();
    for (q=ngbrs_v.begin(); q!=qe; q++) {
      // leq:= number of dof connected to `keq'
      leq = *q;
      if (debug_compute_prof) printf("%d ",leq);
      // if (leq==keq) diag_ok=1;
      if (part.processor(leq) == myrank) {
	// Count in `diagonal' part (i.e. `leq' in the same processor
	// than `keq')
	d_nnz[k]++;
      } else {
	// Count in `off-diagonal' part (i.e. `leq' NOT in the same
	// processor than `keq')
	o_nnz[k]++;
      }	
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

  // Print statistics
  double avo,avd,avdcorr;
  avo = double(sumo)/double(neqp); // Average off-diag
  avd = double(sumd)/double(neqp); // Average diag
  avdcorr = double(sumdcorr)/double(neqp); // Average corrected
					   // (Used no more)
  PetscSynchronizedPrintf(comm,
			  "On processor %d,\n"
			  "       diagonal block terms: %d, (%f av.)\n"
			  // Corrected does not make sense anymore
			  // since it seems that PETSc does not need
			  // the diagonal terms. 
			  // "                (corrected): %d, (%f av.)\n"
			  "   off diagonal block terms: %d, (%f av.)\n",
			  // myrank,sumd,avd,sumdcorr,avdcorr,sumo,avo);
			  myrank,sumd,avd,sumo,avo);
  PetscSynchronizedFlush(comm);
  
  // Create matrices
  ierr =  MatCreateMPIAIJ(comm,neqp,neqp,neq,neq,
			  PETSC_NULL,d_nnz,PETSC_NULL,o_nnz,&A);
  CHKERRQ(ierr); 
  ierr =  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR);
  CHKERRQ(ierr); 
  // P and A are pointers (in PETSc), otherwise this may be somewhat risky
  P=A;
  delete[] d_nnz;
  delete[] o_nnz;
  return 0;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::clear"
void PETScMat::clear() {
  ierr = PFPETScMat::clear(); CHKERRQ(ierr); 
  // P is not destroyed, since P points to A
  int ierr = MatDestroy_maybe(A); CHKERRQ(ierr); 
  CHKERRQ(ierr); 
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::factor_and_solve"
int PETScMat::factor_and_solve_a(Vec &res,Vec &dx) {
  ierr = build_sles(); CHKERRQ(ierr); 
  ierr = solve_only_a(res,dx); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::solve_only"
int PETScMat::solve_only_a(Vec &res,Vec &dx) {
  int ierr = SLESSolve(sles,res,dx,&its_); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::view"
int PETScMat::view(Viewer viewer=VIEWER_STDOUT_WORLD) {
  ierr = MatView(A,viewer); CHKERRQ(ierr); 
//    ierr = SLESView(sles,VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PETScMat::clean_mat_a"
int PETScMat::clean_mat_a() {
  ierr=MatZeroEntries(A); CHKERRQ(ierr);
  return 0;
};

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/home/mstorti/PETSC/petscfem/tools/pfcpp")
  End: 
*/
