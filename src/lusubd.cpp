//__INSERT_LICENSE__
//$Id: lusubd.cpp,v 1.8 2001/07/18 13:59:32 mstorti Exp $

#include <typeinfo>
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

const int IISDMat::D=0;
const int IISDMat::O=1;
const int IISDMat::L=0;
const int IISDMat::I=1;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int IISD_mult(Mat A,Vec x,Vec y)"
int IISD_mult(Mat A,Vec x,Vec y) {
  void *ctx;
  IISDMat *pfA;
  int ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr); 
  pfA = (IISDMat *) ctx;
  pfA->mult(x,y);
  return 0;
}

void IISDMat::create(Darray *da,const Dofmap *dofmap_,
		int debug_compute_prof) {
  int myrank,size;
  int k,pos,keq,leq,jj,row,row_t,col_t,od,
    d_nz,o_nz,nrows,ierr,n_loc_h,n_int_h,k1h,k2h,rank;
  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);
  MPI_Comm_size (PETSC_COMM_WORLD, &size);

  dofmap = dofmap_;
  const int &neq = dofmap->neq;
  
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

  // map:= map[k] is the location of dof `k1+k' in reordering such
  // that the `local' dofs are the first `n_loc' and the `interface'
  // dofs are the last n_int.
  // n_int := number of nodes in this processor that are `interface'
  // n_loc := number of nodes in this processor that are `local'
  map.resize(neq,0);
  n_int_v.resize(size+1,0);
  n_loc_v.resize(size+1,0);

  // number all `loc' dof's
  n_loc_tot = 0;
  n_loc_v[0] = 0;
  for (rank = 0; rank < size; rank++) {
    // PetscPrintf(PETSC_COMM_WORLD,"startproc %d, neqproc %d\n",
    // dofmap->startproc[rank],dofmap->neqproc[rank]);
    k1h = dofmap->startproc[rank];
    for (k = 0; k < dofmap->neqproc[rank]; k++) {
      // PetscPrintf(PETSC_COMM_WORLD,"dof %d, flag %d\n",k1h+k,flag[k1h+k]);
      if(flag[k1h+k] == 0) map[k1h+k] = n_loc_tot++;
    }
    n_loc_v[rank+1] = n_loc_tot;
  }

  // number all `int' dof's
  n_int_tot = n_loc_tot;
  n_int_v[0] = n_loc_tot;
  for (rank = 0; rank < size; rank++) {
    k1h = dofmap->startproc[rank];
    for (k = 0; k < dofmap->neqproc[rank]; k++) 
      if(flag[k1h+k] == 1) map[k1h+k] = n_int_tot++;
    n_int_v[rank+1] = n_int_tot;
  }
  PETSCFEM_ASSERT0(n_int_tot == neq,"Failed to count all dof's in "
		   "int-loc partitioning");   
  // Correct number of total `int' dof's
  n_int_tot = n_int_tot - n_loc_tot;
  n_int = n_int_v[myrank+1]-n_int_v[myrank];
  n_loc = n_loc_v[myrank+1]-n_loc_v[myrank];
  n_locp = n_loc_v[myrank];

#if 0 // Debug
  int ldof,type;
  if (myrank==0) {
    for (rank=0; rank<size; rank++) 
      printf("[%d] loc: %d-%d, int: %d-%d\n",
	     rank,n_loc_v[rank],n_loc_v[rank+1],
	     n_int_v[rank],n_int_v[rank+1]);
    for (rank=0; rank<size; rank++) {
      printf("In processor [%d]: loc %d, int %d\n",rank,
	     n_loc_v[rank+1]-n_loc_v[rank],
	     n_int_v[rank+1]-n_int_v[rank]);
      for (k=0; k<dofmap->neqproc[rank]; k++) {
	keq = dofmap->startproc[rank]+k;
	ldof = map[keq];
	if (ldof < n_loc_tot) {
	  type = L;
	  ldof = ldof - n_loc_v[rank];
	} else {
	  type = I;
	  ldof = ldof - n_int_v[rank];
	}
	printf("keq: %d, type: %s, local: %d\n",
	       keq,(type==L ? "L" : "I"),ldof);
      }
    }
  }
  PetscFinalize();
  exit(0);
#endif

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
    // row_t = (flag[keq] ? I : L);
    // Index in the local PETSc matrices (maped index)
    row = map[keq];
    row_t = (row < n_loc_tot ? L : I);
    
    // Correct dof's
    if (row >= n_loc_tot) row = row - n_int_v[myrank];
    row = row - n_loc_v[myrank];
    // loop over the connected dof's
    pos = keq;
    while (1) {
      nodep = (Node *)da_ref(da,pos);
      if (nodep->next==-1) break;
      // leq:= number of dof connected to `keq' i.e. `A(keq,leq) != 0' 
      leq = nodep->val;
      // type of dof
      // col_t = (flag[leq] ? I : L);
      col_t = (map[leq] < n_loc_tot ? L : I);
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
  PetscFinalize();
  exit(0);
#endif
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_loc,n_loc,PETSC_NULL,
			 nnz[D][L][L].begin(),&A_LL); 
  PETSCFEM_ASSERT0(ierr==0,"Error creating loc-loc matrix\n"); 

  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_loc,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][L][I].begin(),
			 PETSC_NULL,nnz[O][L][I].begin(),
			 &A_LI);
  PETSCFEM_ASSERT0(ierr==0,"Error creating loc-int matrix\n"); 
    
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_int,n_loc,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][L].begin(),
			 PETSC_NULL,nnz[O][I][L].begin(),
			 &A_IL);
  PETSCFEM_ASSERT0(ierr==0,"Error creating int-loc matrix\n"); 
    
  ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,n_int,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][I].begin(),
			 PETSC_NULL,nnz[O][I][I].begin(),
			 &A_II);
  PETSCFEM_ASSERT0(ierr==0,"Error creating int-int matrix\n"); 
  
  // extern int mult(Mat,Vec,Vec);
  ierr = MatCreateShell(PETSC_COMM_WORLD,n_int,n_int,
			PETSC_DETERMINE,PETSC_DETERMINE,this,&A);
  PETSCFEM_ASSERT0(ierr==0,"Error creating shell matrix\n"); 
  P=A;
  MatShellSetOperation(A,MATOP_MULT,(void *)(&IISD_mult));
  ierr = VecCreateMPI(PETSC_COMM_WORLD,n_loc,PETSC_DETERMINE,&x_loc);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc' vector\n"); 
  ierr = VecCreateSeq(PETSC_COMM_SELF,n_loc,&y_loc_seq);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `y_loc_seq' vector\n"); 
  ierr = VecDuplicate(y_loc_seq,&x_loc_seq);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc_seq' vector\n"); 

  ierr = SLESCreate(PETSC_COMM_SELF,&sles_ll); 
  PETSCFEM_ASSERT0(ierr==0,"Error creating SLES.\n"); 
  ierr = SLESSetOperators(sles_ll,A_LL,
			  A_LL,SAME_NONZERO_PATTERN);
  PETSCFEM_ASSERT0(ierr==0,"Error with `SLESSetOperators'.\n"); 
  ierr = SLESGetKSP(sles_ll,&ksp_ll); 
  PETSCFEM_ASSERT0(ierr==0,"Error with `SLESGetKSP'.\n"); 
  ierr = SLESGetPC(sles_ll,&pc_ll);
  PETSCFEM_ASSERT0(ierr==0,"Error with `SLESGetPC'.\n"); 

  ierr = KSPSetType(ksp_ll,KSPGMRES);
  PETSCFEM_ASSERT0(ierr==0,"Error with `KSPSetType'.\n"); 

  ierr = KSPSetTolerances(ksp_ll,0,0,1e10,1);
  PETSCFEM_ASSERT0(ierr==0,"Error setting tolerances.\n"); 

  ierr = PCSetType(pc_ll,PCLU);
  PETSCFEM_ASSERT0(ierr==0,"Error setting PC type.\n"); 
  ierr = KSPSetMonitor(ksp_ll,PETSC_NULL,PETSC_NULL);

  // Shortcuts
  AA[L][L] = &A_LL;
  AA[L][I] = &A_LI;
  AA[I][L] = &A_IL;
  AA[I][I] = &A_II;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void IISDMat::mult(Vec x,Vec y)"
void IISDMat::mult(Vec x,Vec y) {

  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  const int &neqp = dofmap->neqproc[myrank];
  // We are solving the kernel of the linear system
  //
  // ALL XL + ALI XI = 0
  // AIL XL + AII XI = RI
  //
  // XI comes in x (interface nodes) and we have to compute y = RI 
  
  int j,ierr,its_;
  double *a,*aa;;

  // x_loc <- A_LI * XI
  ierr = MatMult(A_LI,x,x_loc);
  PETSCFEM_ASSERT0(ierr==0,"Error performing `A_LI*x' product.\n"); 

  // Localize on procesor and make y_loc_seq <- - A_LI * XI
  ierr = VecGetArray(x_loc,&a);
  PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray'.\n"); 

  ierr = VecGetArray(y_loc_seq,&aa);
  PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray'.\n"); 

  for (j=0; j<neqp; j++) aa[j] = -a[j];

  ierr = VecRestoreArray(y_loc_seq,&aa);

  // Solve local system: x_loc_seq <- XL
  ierr = SLESSolve(sles_ll,y_loc_seq,x_loc_seq,&its_); 
  PETSCFEM_ASSERT0(ierr==0,"Error solving subdomain problem.\n"); 
  
  // Pass to global vector: x_loc <- XL
  ierr = VecGetArray(x_loc_seq,&aa);
  PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray'.\n"); 

  for (j=0; j<neqp; j++) a[j] = aa[j];

  ierr = VecRestoreArray(x_loc_seq,&aa);
  ierr = VecRestoreArray(x_loc,&a);

  // 
  ierr = MatMult(A_II,x,y);
  PETSCFEM_ASSERT0(ierr==0,"Error performing `A_II*x' product.\n"); 

  ierr = MatMultAdd(A_IL,x_loc,y,y);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int IISDMat::assembly_begin(MatAssemblyType type)"
int IISDMat::assembly_begin(MatAssemblyType type) {
  int ierr;
  ierr = MatAssemblyBegin(A_LL,type); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_IL,type); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_LI,type); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A_II,type); CHKERRQ(ierr);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int IISDMat::assembly_end(MatAssemblyType type)"
int IISDMat::assembly_end(MatAssemblyType type) {
  int ierr;
  ierr = MatAssemblyEnd(A_LL,type); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_IL,type); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_LI,type); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A_II,type); CHKERRQ(ierr);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int IISDMat::zero_entries()"
int IISDMat::zero_entries() {
  int ierr;
  ierr=MatZeroEntries(A_LL); CHKERRQ(ierr);
  ierr=MatZeroEntries(A_IL); CHKERRQ(ierr);
  ierr=MatZeroEntries(A_LI); CHKERRQ(ierr);
  ierr=MatZeroEntries(A_II); CHKERRQ(ierr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int IISDMat::view()"
int IISDMat::view(Viewer viewer) {
  int ierr;
  PetscPrintf(PETSC_COMM_WORLD," IISD matrix - printing [LOC][LOC] block\n");
  ierr = MatView(A_LL,viewer); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD," IISD matrix - printing [LOC][INT] block\n");
  ierr = MatView(A_LI,viewer); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD," IISD matrix - printing [INT][LOC] block\n");
  ierr = MatView(A_IL,viewer); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD," IISD matrix - printing [INT][INT] block\n");
  ierr = MatView(A_II,viewer); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD," IISD matrix - printing [INT][INT] block\n");

  ierr =  SLESView(sles,VIEWER_STDOUT_SELF);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void IISDMat::clear()"
void IISDMat::clear() {
  // P is not destroyed, since P points to A
  PFMat::clear();
  int ierr = MatDestroy(A_LL); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_LL (loc-loc)\n");
  ierr = MatDestroy(A_LI); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_LI (loc-int)\n");
  ierr = MatDestroy(A_IL); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_IL (int-loc)\n");
  ierr = MatDestroy(A_II); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix A_II (int-int)\n");

  ierr = MatDestroy(A); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix shell A\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void IISDMat::map_dof(int gdof,int &block,int &ldof)"
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
#define __FUNC__ ""
void IISDMat::set_value(int row,int col,Scalar value,
			InsertMode mode=ADD_VALUES) {
  int row_indx,col_indx,row_t,col_t;
  map_dof(row,row_t,row_indx);
  map_dof(col,col_t,col_indx);
  MatSetValues(*(AA[row_t][col_t]),1,&row_indx,1,&col_indx,&value,mode);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void IISDMat::solve(Vec res,Vec dx)"
void IISDMat::solve(Vec res,Vec dx) {
  int ierr,kloc;
  double *res_a,*y_loc_seq_a,*x_loc_seq_a,*dx_a;
  if (n_int_tot > 0 ) {
    ierr = SLESSolve(sles,res,dx,&its_); 
    PETSCFEM_ASSERT0(ierr==0,"Error solving linear system"); 
  } else {
    ierr = VecGetArray(res,&res_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray' on `res'.\n"); 
#if 0
    ierr = VecGetArray(y_loc_seq,&y_loc_seq_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray' on `y_loc_seq'.\n"); 
#endif

    for (int j = 0; j < neqp; j++) {
      kloc = map[k1+j] - n_locp;
      y_loc_seq_a[kloc] = res_a[k1+j];
    }
    
    ierr = VecRestoreArray(res,&res_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecRestoreArray' on `res'.\n"); 
    ierr = VecRestoreArray(y_loc_seq,&y_loc_seq_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecRestoreArray' on `y_loc_seq'.\n"); 
    
    if (n_loc > 0) {
      ierr = SLESSolve(sles_ll,y_loc_seq,x_loc_seq,&its_); 
      PETSCFEM_ASSERT0(ierr==0,"Error solving local system"); 
    }
#if 1
    ierr = VecView(y_loc_seq,VIEWER_STDOUT_SELF);
    ierr = VecView(x_loc_seq,VIEWER_STDOUT_SELF);
    ierr = MatView(A_LL,VIEWER_STDOUT_SELF);

    PetscFinalize();
    exit(0);
#endif

    ierr = VecGetArray(dx,&dx_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray' on `dx'.\n"); 
    ierr = VecGetArray(x_loc_seq,&x_loc_seq_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecGetArray' on `x_loc_seq'.\n"); 

    for (int j = 0; j < neqp; j++) {
      kloc = map[k1+j] - n_locp;
      dx_a[k1+j] = x_loc_seq_a[kloc];
    }

    ierr = VecRestoreArray(dx,&dx_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecRestoreArray' on `dx'.\n"); 
    ierr = VecRestoreArray(x_loc_seq,&x_loc_seq_a);
    PETSCFEM_ASSERT0(ierr==0,"Error performing `VecRestoreArray' on `x_loc_seq'.\n"); 
  }
}

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
  PFMat::clear();
  // P is not destroyed, since P points to A
  int ierr = MatDestroy(A); 
  PETSCFEM_ASSERT0(ierr==0,"Error destroying PETSc matrix\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int PETScMat::build_sles(TextHashTable *,char *=NULL)"
void PFMat::build_sles(TextHashTable *thash,char *name=NULL) {

  static int warn_iisdmat=0;
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

  if (!warn_iisdmat && typeid(*this)==typeid(IISDMat) 
      && preco_type != "none") {
    warn_iisdmat=1;
    PetscPrintf(PETSC_COMM_WORLD,
		"PETScFEM warning: IISD operator does not support any\n"
		"preconditioning. Switching to \"none\"\n");
    preco_type = "none";
  }

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int view(Viewer viewer)"
int PETScMat::view(Viewer viewer) {
  ierr = MatView(A,viewer); CHKERRQ(ierr); 
  ierr = SLESView(sles,VIEWER_STDOUT_SELF); CHKERRQ(ierr); 
};

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "~/PETSC/petscfem/tools/pfcpp")
  End: 
*/
