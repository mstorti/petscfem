//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.4 2001/07/14 10:53:24 mstorti Exp $

#include "lusubd.h"

PFMat::~PFMat() {};

PFMatLU::PFMatLU(Dofmap *dofmap_,Darray *da) {
  int myrank,size;
  int k,k1,k2,neqp,pos,keq,leq,jj,row,row_t,col_t,od,
    d_nz,o_nz,nrows,ierr;
  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);
  MPI_Comm_size (PETSC_COMM_WORLD, &size);

  dofmap = dofmap_;
  
  int &neq=dofmap->neq;
  Node *nodep;
  // neqp:= k2:= k1:= unknowns in this processor are thos in the range
  // k1 <= k <= k2 = k1+neqp-1
  k1 = dofmap->startproc[myrank];
  neqp = dofmap->neqproc[myrank];
  k2=k1+neqp-1;

  // these are the vectors for use with the PETSc matrix constructors
  // flag:= will contain `0' for local dof's and `1' for `interface'
  // dof's. 
  // nnz:= eight vectors that contain the d_nnz, o_nnz for
  // dimensioning PETSc matrices
  // nnz[diagonal/off-diagonal][row-block-type][column-block-type]
  // For instance nnz[D][L][L] = d_nn_z for the 'local-local' block. 
  const int D=0,O=1,L=0,I=1;
  vector<int> nnz[2][2][2],flag0(dofmap->neq,0),flag(dofmap->neq,0),
    *d_nnz,*o_nnz;

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
	// 	printf("[%d] marking node %d\n",keq);
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
  // that the `local' dofs are the first `n_loc' and the `interface'
  // dofs are the last n_int.
  // n_int := number of nodes in this processor that are `interface'
  // n_loc := number of nodes in this processor that are `local'
  map.resize(neqp,0);
  n_loc=0;
  for (k = 0; k < neqp; k++) if(flag[k1+k] == 0) map[k] = n_loc++;
  jj=n_loc;
  for (k = 0; k < neqp; k++) if(flag[k1+k] == 1) map[k] = jj++;
  assert(jj == neqp);

  n_int = neqp - n_loc;

  // Now we have to construct the `d_nnz' and `o_nnz' vectors
  // od:= may be `D' (0) or `I' (1). Diagonal or off-diagonal (in the
  // PETSc sense)
  // row_t:= col_t:= may be local (`L=0') or `interface ('I=1') 
  // is a a block index, when decomposing the dof's
  // at each processor as `local' and `interface'. 
  for (od = 0; od < 2; od++) {
    for (row_t = 0; row_t < 2; row_t++) {
      for (col_t = 0; col_t < 2; col_t++) {
	// The size of the `d_nnz' and `o_nnz' vectors is that of
	// the `row' type index, i.e. the second index, for instance
	// the size of the nnz[O][L][I] index is that of the `L'
	// block, i.e. `n_loc'. 
	nnz[od][row_t][col_t].resize((row_t == L ? n_loc : n_int),0);
      }
    }
  }
    
  // For each dof in this processor we scan all connected dof's and
  // add the corresponding element in the `d_nnz' or `o_nnz' vectors. 
  for (k = 0; k < neqp; k++) {
    // keq:= number of dof
    keq = k1 + k;
    // type of dof (local or interface)
    row_t = (flag[keq] ? I : L);
    // Index in the PETSc matrices (maped index)
    row = map[k];
    // Correct for interface dof's
    if (row > n_loc) row = row - n_loc;
    // loop over the connected dof's
    pos = keq;
    while (1) {
      nodep = (Node *)da_ref(da,pos);
      if (nodep->next==-1) break;
      // leq:= number of dof connected to `keq' i.e. `A(keq,leq) != 0' 
      leq = nodep->val;
      // type of dof
      col_t = (flag[leq] ? I : L);
      // diagonal or off-diagonal (in PETSc sense)
      od = ((leq < k1 || leq > k2) ? O : D);
      // By construction, the local-local block should be in the
      // diagonal part
      assert(!(od==O && row_t==L && col_t==L));
      // count 
      nnz[od][row_t][col_t][row]++;
      // next link
      pos = nodep->next;
    }
  }

#if 0
  // Prints d_nnz, o_nnz for block LL, IL, IL and II in turn
  // For each block prints the d_nnz and o_nnz in turn
  for (int row_t=0; row_t<2; row_t++) {
    for (int col_t=0; col_t<2; col_t++) {
      PetscPrintf(PETSC_COMM_WORLD,"[%s]-[%s] block\n",
		  (row_t==0? "LOC" : "INT"),(col_t==0? "LOC" :
					     "INT"));
      // number of rows in this block in this processor
      nrows=(row_t==L ? n_loc : n_int);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "%d/%d rows on processor [%d]\n",
			      nrows,neqp,myrank);
      for (row=0; row<nrows; row++) {
	d_nz = nnz[D][row_t][col_t][row];
	o_nz = nnz[O][row_t][col_t][row];
	if (d_nz|| o_nz) 
	  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				  "row=%d, d_nnz=%d, o_nnz=%d\n",row,
				  d_nz,o_nz);
      }
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
    }
  }
#endif
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_loc,n_loc,PETSC_NULL,
			 nnz[D][L][L].begin(),&A_LL); 
  PETSCFEM_ASSERT0(ierr==0,"Error creating loc-loc matrix"); 

  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_loc,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][L][I].begin(),
			 PETSC_NULL,nnz[O][L][I].begin(),
			 &A_LI);
  PETSCFEM_ASSERT0(ierr==0,"Error creating loc-int matrix"); 
    
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_int,n_loc,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][L].begin(),
			 PETSC_NULL,nnz[O][I][L].begin(),
			 &A_IL);
  PETSCFEM_ASSERT0(ierr==0,"Error creating int-loc matrix"); 
    
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_int,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][I].begin(),
			 PETSC_NULL,nnz[O][I][I].begin(),
			 &A_II);
  PETSCFEM_ASSERT0(ierr==0,"Error creating int-int matrix"); 
    
  PetscFinalize();
  exit(0);
}
