//__INSERT_LICENSE__
//$Id: petscmat.cpp,v 1.1.2.1 2001/12/27 19:55:47 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;
#include <typeinfo>
#include <mat.h>

#include <src/libretto.h>

#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/petscmat.h>
#include <src/graph.h>

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
#undef __FUNC__
#define __FUNC__ "PETScMat::create"
void PETScMat::create() {
  int k,k1,k2,neqp,keq,leq,pos,sumd=0,sumdcorr=0,sumo=0,ierr,myrank;

  MPI_Comm_rank (comm, &myrank);
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
  int neq=dofmap->neq;
  ierr =  MatCreateMPIAIJ(comm,dofmap->neqproc[myrank],
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
#define __FUNC__ "PETScMat::clear"
void PETScMat::clear() {
  PFMat::clear();
  // P is not destroyed, since P points to A
  if (A) {
    int ierr = MatDestroy(A); 
    PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix \"A\"\n");
  }
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
