//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.3 2001/07/13 03:48:23 mstorti Exp $

#include "lusubd.h"

PFMat::~PFMat() {};

PFMatLU::PFMatLU(Dofmap *dofmap_,Darray *da) {
  int myrank;
  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);

  dofmap = dofmap_;
  int k,k1,k2,neqp,pos,keq,leq,nloc,nint,jj;
  int &neq=dofmap->neq;
  Node *nodep;
  // neqp:= k2:= k1:= unknowns in this processor are thos in the range
  // k1 <= k <= k2 = k1+neqp-1
  k1 = dofmap->startproc[myrank];
  neqp = dofmap->neqproc[myrank];
  k2=k1+neqp-1;

  // these are the vectors for use with the PETSc matrix constructors
  // falg := will contain `0' for local dof's and `1' for `interface' dof's
  vector<int> d_nnz_ll,d_nnz_li,d_nnz_il,d_nnz_ii,
    o_nnz_ll,o_nnz_li,o_nnz_il,o_nnz_ii,flag0(dofmap->neq,0),flag(dofmap->neq,0);

  // First, decide which dof's are marked as `interface' and which as `local'. 

  // A dof in processor `k' is marked as interface if it is connected
  // (has a non zero matrix value) with a dof in a processor with a lower index "k'<k"
  for (k=0;k<neqp;k++) {
    // keq:= number of dof
    keq=k1+k;
    pos=keq;
    while (1) {
      nodep = (Node *)da_ref(da,pos);
      if (nodep->next==-1) break;
      // leq:= number of dof connected to `keq'
      leq = nodep->val;
      // if leq is in other processor, then either `keq' of `leq' are
      // interface. This depends on `leq<keq' or `leq>keq'
      if (leq<k1) {
	printf("[%d] marking node %d\n",keq);
	flag0[keq]=1;
      } else if (leq>k2) {
	flag0[leq]=1;
	// printf("[%d] marking node %d\n",myrank,leq);
      }
      pos = nodep->next;
    }
  }
  
  // Each processor marks dof's in that belongs to them and to other
  // so that, after, we have to combine all them with an Allreduce
  MPI_Allreduce(flag0.begin(), flag.begin(), neq, MPI_INT, 
		MPI_MAX, PETSC_COMM_WORLD);
  flag0.clear();
//    if (myrank==0) {
//      for (k=0; k<neq; k++) 
//        printf("dof %d: flag %d\n",k1+k,flag[k]);
//    }
//    PetscFinalize();
//    exit(0);
  // map:= map[k] is the location of dof `k1+k' in reordering such
  // that the `local' dofs are the first `nloc' and the `interface'
  // dofs are the last nint.
  // n_int := number of nodes in this processor that are `interface'
  // n_loc := number of nodes in this processor that are `local'
  map.resize(neqp,0);
  nloc=0;
  for (k=0;k<neqp;k++) if(flag[k]==0) map[k] = nloc++;
  jj=nloc;
  for (k=0;k<neqp;k++) if(flag[k]==1) map[k] = jj++;
  assert(jj==neqp);

  nint = neqp-nloc;
  flag.clear();

  // Now we have to construct the `d_nnz' and `o_nnz' vectors
  d_nnz_ll.resize(n_loc,0);
    d_nnz_li = 
    d_nnz_il = 
    d_nnz_ii = 
    o_nnz_ll = 
    o_nnz_li = 
    o_nnz_il = 
    o_nnz_ii = 
  
}
