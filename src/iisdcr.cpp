//__INSERT_LICENSE__
//$Id: iisdcr.cpp,v 1.6 2001/11/27 14:34:45 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

//  #define DEBUG_IISD
//  #define DEBUG_IISD_DONT_SET_VALUES

#include <src/debug.h>
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class IISDGraph : public Graph {
private:
  Node *nodep;
  int vrtxf;
public:
  const int *flag;
  int k1,k2;
  /// Auxiliary functions
  const int *dof2loc,*loc2dof;
  /// Libretto dynamic array that contains the graph adjacency matrix
  Darray *da;
  /// callback user function to return the neighbors for a 
  void set_ngbrs(int vrtx_f,vector<int> &ngbrs_v);
  /// Callback user function for the user to set the weight of a given fine vertex. 
  double weight(int vrtx_f);
  /// Clean all memory related 
  ~IISDGraph() {clear();}
  /// Constructor
  IISDGraph() : Graph() {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void IISDGraph::set_ngbrs(int loc1,vector<int> &ngbrs_v) {
  int pos,loc2,dof2;
  pos = loc2dof[loc1]+k1;
  while (1) {
    nodep = (Node *)da_ref(da,pos);
    if (nodep->next==-1) break;
    // loc2:= number of global dof connected to `'
    dof2 = nodep->val;
    if (k1<=dof2 && dof2<=k2 && !flag[dof2] ) {
      loc2 = dof2loc[dof2-k1];
      ngbrs_v.push_back(loc2);
    }
    pos = nodep->next;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double IISDGraph::weight(int elem) {
  return 1.;
}

//---:---<*>---:---<*>---:a---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISD_mult"
int IISD_mult(Mat A,Vec x,Vec y) {
  void *ctx;
  IISDMat *pfA;
  int ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr); 
  pfA = (IISDMat *) ctx;
  ierr = pfA->mult(x,y); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISD_mult_trans"
int IISD_mult_trans(Mat A,Vec x,Vec y) {
  void *ctx;
  IISDMat *pfA;
  int ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr); 
  pfA = (IISDMat *) ctx;
  ierr = pfA->mult_trans(x,y); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::create"
void IISDMat::create(Darray *da,const Dofmap *dofmap_,
		int debug_compute_prof) {

  int myrank,size,max_partgraph_vertices;
  int k,pos,keq,leq,jj,row,row_t,col_t,od,
    d_nz,o_nz,nrows,ierr,n_loc_h,n_int_h,k1h,k2h,rank,
    n_loc_pre,loc,dof,subdoj,subdok,vrtx_k;
  vector<int> dof2loc,loc2dof,ngbrs_v;
  vector<int>::iterator q,qe;
  IISDGraph graph;

  MPI_Comm_rank (PETSC_COMM_WORLD, &myrank);
  MPI_Comm_size (PETSC_COMM_WORLD, &size);

  dofmap = dofmap_;
  part = new DofmapPartitioner(dofmap);
  A_LL_other = new DistMatrix(part,PETSC_COMM_WORLD);
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
  vector<int> nnz[2][2][2],flag0,flag;
  flag0.resize(dofmap->neq,0);
  flag.resize(dofmap->neq,0);

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
      // if leq is in other processor, then either `keq' or `leq' are
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
  
  // Each processor marks as interface dof's that belongs to them and
  // to other so that, after, we have to combine all them with an
  // Allreduce
  MPI_Allreduce(flag0.begin(), flag.begin(), neq, MPI_INT, 
		MPI_MAX, PETSC_COMM_WORLD);

#if 1
  // SUBPARTITIONING. //---:---<*>---:---<*>---:---<*>---:---<*>---:
  // We partition the graph of those vertices (dof's) that are local
  // to this processor. Those that result in internal interfaces are
  // then marked as interface. So that the total effect of this step
  // is to add some local nodes to the `interface' set.

  // Build the `dof2loc' and `loc2dof' maps. 
  // n_loc_pre:= the number of local nodes (so far). After this step
  // the number of local vertices is incremented.
  // dof2loc:= maps dofs in this processor (0<dof<neqp) to local dof's
  // (0<loc<n_loc_pre)
  // loc2dof:= the inverse of dof2loc
  n_loc_pre = 0;
  dof2loc.resize(neqp,-1);
  for (k=0;k<neqp;k++) {
    if (flag[k1+k]==0) dof2loc[k] = n_loc_pre++;
  }
  loc2dof.resize(n_loc_pre);
  for (dof = 0; dof < neqp; dof++) {
    loc = dof2loc[dof];
    if (loc>=0) loc2dof[loc] = dof;
  }
  
  // Fill Graph class members 
  graph.da = da;
  graph.dof2loc = dof2loc.begin();
  graph.loc2dof = loc2dof.begin();
  graph.k1 = k1;
  graph.k2 = k2;
  graph.flag = flag.begin();

#define INF INT_MAX
  //o The maximum number of vertices in the coarse mesh. 
  TGETOPTDEF_ND_PFMAT(&thash,int,max_partgraph_vertices,INF);
#undef INF
  TGETOPTDEF_ND_PFMAT(&thash,int,iisd_subpart,1);
  // assert(iisd_subpart!=1);

  graph.part(n_loc_pre,max_partgraph_vertices,iisd_subpart);
  // Mark those local dofs that are connected to a local dof in a
  // subdomain with lower index in the subpartitioning as interface.
  for (k=0; k<n_loc_pre; k++) {
    subdoj = graph.vrtx_part(k);
    ngbrs_v.clear();
    graph.set_ngbrs(k,ngbrs_v);
    qe = ngbrs_v.end();
    for (q=ngbrs_v.begin(); q!=qe; q++) {
      if (graph.vrtx_part(*q)<subdoj) {
	printf("[%d] marking %d as interface\n",myrank,k1+loc2dof[k]);
#if 0 // debug:=
	int kkk = k1+loc2dof[k];
	assert(0<= kkk &&kkk<neq);
#endif
	flag[k1+loc2dof[k]] = 1;
	break;
      }
    }
  }

  debug.trace("after remarking flag");
  graph.clear();
  dof2loc.clear();
  loc2dof.clear();

  // We have to combine all them again with an Allreduce
  debug.trace("before allreduce");
  MPI_Allreduce(flag.begin(), flag0.begin(), neq, MPI_INT, 
		MPI_MAX, PETSC_COMM_WORLD);
  debug.trace("after allreduce");
  // recopy on `flag'...
  memcpy(flag.begin(),flag0.begin(),neq*sizeof(int));
  flag0.clear();
#endif

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
  n_intp = n_int_v[myrank];

#ifdef DEBUG_IISD // Debug
  int ldof,type;
  if (myrank==0) {
    for (rank=0; rank<size; rank++) 
      printf("[%d] loc: %d-%d, int: %d-%d\n",
	     rank,n_loc_v[rank],n_loc_v[rank+1],
	     n_int_v[rank],n_int_v[rank+1]);
#if 0
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
	printf("keq: %d, map: %d, type: %s, local: %d\n",
	       keq,map[keq],(type==L ? "L" : "I"),ldof);
      }
    }
#endif
  }
#endif

  debug.trace("trace 0");
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
    
  debug.trace("trace 0.1");
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
    row -= (row < n_loc_tot ? n_locp : n_intp);
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
      if (!(row>=0 && row < nnz[od][row_t][col_t].size())) {
	printf("row %d, size: %d\n",row,nnz[od][row_t][col_t].size());
	MPI_Abort(PETSC_COMM_WORLD,1);
      }
      nnz[od][row_t][col_t][row]++;
      // next link
      pos = nodep->next;
    }
  }

  debug.trace("trace 1");
  // deallocate Libretto dynamic darray
  da_destroy(da);

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

  debug.trace("trace 2");
//    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,n_loc,n_loc,PETSC_NULL,
//   			 nnz[D][L][L].begin(),&A_LL); 
//    PETSCFEM_ASSERT0(ierr==0,"Error creating loc-loc matrix\n"); 
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
  
  debug.trace("trace 3");
  // extern int mult(Mat,Vec,Vec);
  ierr = MatCreateShell(PETSC_COMM_WORLD,n_int,n_int,
			PETSC_DETERMINE,PETSC_DETERMINE,this,&A);
  PETSCFEM_ASSERT0(ierr==0,"Error creating shell matrix\n"); 
  P=A;

  MatShellSetOperation(A,MATOP_MULT,(void *)(&IISD_mult));
  MatShellSetOperation(A,MATOP_MULT_TRANS,(void *)(&IISD_mult_trans));

  ierr = VecCreateMPI(PETSC_COMM_WORLD,n_loc,PETSC_DETERMINE,&x_loc);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc' vector\n"); 
  ierr = VecCreateSeq(PETSC_COMM_SELF,n_loc,&y_loc_seq);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `y_loc_seq' vector\n"); 
  ierr = VecDuplicate(y_loc_seq,&x_loc_seq);
  PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc_seq' vector\n"); 

  // Shortcuts
  AA[L][L] = &A_LL;
  AA[L][I] = &A_LI;
  AA[I][L] = &A_IL;
  AA[I][I] = &A_II;

  // Save copy of d_nnz_LL for use later when recreating A_LL
  d_nnz_LL = nnz[D][L][L];

#if 0
  // Build the interface preco stuff
  int_layers.clear();
  set<int> queue;
  for (k=k1; k<k2; k++) {
    row_t = map_dof(I->first,row_t,row_indx);
    if (row_t == L) ;
  }
#endif
}
