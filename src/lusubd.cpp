//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.2 2001/07/12 20:10:26 mstorti Exp $

#include <vector>

#include "fem.h"
#include "utils.h"

class PFMat {
public:
  virtual void set_profile(DArray *da)=0;
  virtual ~PFMat()=0;
}

class PFMatLu {
  Dofmap *dofmap;
public:
  PFMat(Dofmap *dofmap_,DArray *da);
}

PFMat::PFMat(Dofmap *dofmap_,DArray *da) {
  dofmap = dofmap_;
  int k,k1,k2;
  Node *nodep;
  // neqp:= k2:= k1:= unknowns in this processor are thos in the range
  // k1 <= k <= k2 = k1+neqp-1
  k1 = dofmap->startproc[myrank];
  neqp = dofmap->neqproc[myrank];
  k2=k1+neqp-1;

  // these are the vectors for use with the PETSc matrix constructors
  // falg := will contain `0' for local dof's and `1' for `interface' dof's
  vector<int> map(neqp,0),d_nnz_ll,d_nnz_li,d_nnz_il,d_nnz_ii,
    o_nnz_ll,o_nnz_li,o_nnz_il,o_nnz_ii,flag(neq,0);

  // First, decide which dof's are marked as `interface' and which as `local'. 

  // A dof in processor `k' is marked as interface if it is connected
  // (has a non zero matrix value) with a dof in a processor with a lower index "k'<k"
  // n_int := number of nodes in this processor that are `interface'
  // n_loc := number of nodes in this processor that are `local'
  int n_int=0, n_loc=0;
  for (k=0;k<neqp;k++) {
    keq=k1+k;
    pos=keq;
    diag_ok=0;
    while (1) {
      nodep = (Node *)da_ref(da,pos);
      if (nodep->next==-1) break;
      leq = nodep->val;
      if (leq<k1) {
      pos = nodep->next;
    }
    if (debug_compute_prof) printf("     --    d_nnz %d   o_nnz %d\n",d_nnz[k],o_nnz[k]);
    //d_nnz[k] += 1;
    //o_nnz[k] += 1;
    // To add room for the diagonal entry, even if it doesn't exist
    sumd += d_nnz[k];
    if (!diag_ok) {
      //      printf("No diagonal element for dof %d in proc %d\n",keq,myrank);
      // fixme:= No hago correccion para elementos diagonales
      //d_nnz[k]+=1; 
    }
    sumdcorr += d_nnz[k];
    sumo += o_nnz[k];
  }

  // deallocate darray
  // wait_from_console("antes de destruir da ");
  da_destroy(da);
  // wait_from_console("despues de destruir da ");

  double avo,avd,avdcorr;
  avo = double(sumo)/double(neqp);
  avd = double(sumd)/double(neqp);
  avdcorr = double(sumdcorr)/double(neqp);
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
			  PETSC_NULL,d_nnz,PETSC_NULL,o_nnz,A); CHKERRA(ierr);
  delete[] d_nnz;
  delete[] o_nnz;
  return 0;
}
