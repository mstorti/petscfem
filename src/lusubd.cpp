//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.6 2001/07/16 20:27:48 mstorti Exp $

#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include "mat.h"

#include "fem.h"
#include "dofmap.h"
#include "elemset.h"
#include "pfmat.h"

PFMat::~PFMat() {};

#if 0
void LUsubdMat::create(Darray *da,const Dofmap *dofmap_,
		int debug_compute_prof) {
  int myrank,size;
  int k,k1,k2,neqp,pos,keq,leq,jj,row,row_t,col_t,od,
    d_nz,o_nz,nrows,ierr;
  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);
  MPI_Comm_size (PETSC_COMM_WORLD, &size);

  dofmap = dofmap_;
  
  const int &neq=dofmap->neq;
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
    
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int PFMatPETSc::create(Darray *,Dofmap *,int =0)"
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
  // P and A are pointers (in PETSc), otherwise this may be somewhat risky
  P=A;
  PETSCFEM_ASSERT0(ierr==0,"Error creating PETSc matrix\n");
  delete[] d_nnz;
  delete[] o_nnz;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void PFMat::clear()"
void PFMat::clear() {
  if (sles_was_built) {
    int ierr = SLESDestroy(sles); CHKERRA(ierr);
    PETSCFEM_ASSERT0(ierr==0,"Error destroying SLES\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void PETScMat::clear()"
void PETScMat::clear() {
  // P is not destroyed, since P points to A
  int ierr = MatDestroy(A); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int PETScMat::build_sles(TextHashTable *,char *=NULL)"
void PFMat::build_sles(TextHashTable *thash,char *name=NULL) {

  int ierr;
  //o Absolute tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND(thash,double,atol,1e-6);
  //o Relative tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND(thash,double,rtol,1e-3);
  //o Divergence tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND(thash,double,dtol,1e+3);
  //o Krylov space dimension in solving the monolithic linear
  // system (Newton linear subiteration) by GMRES.
  TGETOPTDEF_ND(thash,int,Krylov_dim,50);
  //o Maximum iteration number in solving the monolithic linear
  // system (Newton linear subiteration).
  TGETOPTDEF_ND(thash,int,maxits,Krylov_dim);
  //o Prints convergence in the solution of the GMRES iteration. 
  TGETOPTDEF_ND(thash,int,print_internal_loop_conv,0);
  //o After computing the analytic Jacobian, Computes the
  // Jacobian in order to verify the analytic one. 

  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S(thash,string,preco_type,jacobi);
  // I had to do this since `c_str()' returns `const char *'
  char *preco_type_ = new char[preco_type.size()+1];
  strcpy(preco_type_,preco_type.c_str());

  ierr = SLESCreate(PETSC_COMM_WORLD,&sles); CHKERRA(ierr);
  ierr = SLESSetOperators(sles,A,
			  P,SAME_NONZERO_PATTERN); CHKERRA(ierr);
  ierr = SLESGetKSP(sles,&ksp); CHKERRA(ierr);
  ierr = SLESGetPC(sles,&pc); CHKERRA(ierr);

  ierr = KSPSetType(ksp,KSPGMRES); CHKERRA(ierr);

  ierr = KSPSetPreconditionerSide(ksp,PC_RIGHT);
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  ierr = KSPGMRESSetRestart(ksp,Krylov_dim); CHKERRA(ierr);
  ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits);

  ierr = PCSetType(pc,preco_type_); CHKERRA(ierr);
  delete[] preco_type_;
  ierr = KSPSetMonitor(ksp,PFMat_default_monitor,this);
  sles_was_built = 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int PFMat_default_monitor(KSP ksp,int n,double rnorm,void *A_) {
  PFMat *A = (PFMat *)A_;
  return A->monitor(ksp,n,rnorm);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
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
#define __FUNC__ "void PETScMat::solve(Vec res,Vec dx)"
void PFMat::solve(Vec res,Vec dx) {
  int ierr = SLESSolve(sles,res,dx,&its_); 
  PETSCFEM_ASSERT0(ierr==0,"Error solving linear system"); 
}

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "~/PETSC/petscfem/tools/pfcpp")
  End: 
*/
