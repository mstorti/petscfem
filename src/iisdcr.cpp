//__INSERT_LICENSE__
//$Id: iisdcr.cpp,v 1.8.4.14 2002/01/09 16:31:07 mstorti Exp $

// fixme:= this may not work in all applications
extern int MY_RANK,SIZE;

#include <src/debug.h>
#include <typeinfo>
#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <mat.h>

#include <src/debug.h>
#include <src/fem.h>
#include <src/utils.h>
#include <src/elemset.h>
#include <src/pfmat.h>
#include <src/iisdmat.h>
#include <src/iisdgraph.h>

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

/** This is an internal auxiliar class for
    sub-partitioning inside each processor
*/
class LocalGraph : public Graph {
public:
  StoreGraph *lgraph;
  const int *loc2dof,*dof2loc,*dofs_proc,*flag;  
  map<int,int> *proc2glob;
  const DofPartitioner *partit;
  int myrank;
  void set_ngbrs(int vrtx_f,set<int> &ngbrs_v);
  /// Auixiliary set ot ngbrs
  set<int> ngbrs_v_aux;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void LocalGraph::set_ngbrs(int loc1,set<int> &ngbrs_v) {
  set<int>::const_iterator q,qe;
  int keq,dof2,loc2;
  // loc1 is a `local' dof numbering. Map to global.
  keq = dofs_proc[ loc2dof[loc1] ];
  // Clear set
  ngbrs_v_aux.clear();
  // Get set of ngbrs (in global numbering)
  lgraph->set_ngbrs(keq,ngbrs_v_aux);
  qe = ngbrs_v_aux.end();
  for (q=ngbrs_v_aux.begin(); q!=qe; q++) {
    // loc2:= number of global dof connected to `'
    dof2 = *q;
    if (partit->processor(dof2)==myrank && !flag[dof2] ) {
      loc2 = dof2loc[ (*proc2glob)[dof2] ];
      ngbrs_v.insert(loc2);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::create"
int IISDMat::create_a() {

  int myrank,size,max_partgraph_vertices_proc,proc_l;
  int k,pos,keq,leq,jj,row,row_t,col_t,od,
    d_nz,o_nz,nrows,ierr,n_loc_h,n_int_h,k1h,k2h,rank,
    n_loc_pre,loc,dof,subdoj,subdok,vrtx_k;
  vector<int> dof2loc,loc2dof;
  set<int> ngbrs_v;
  set<int>::iterator q,qe;
  LocalGraph local_graph;

  // this is a trick to avoid the collision of `local_solver' both
  // as member and as string-option here
  string local_solver_s;
  { string &local_solver = local_solver_s;
  //o Chooses the local solver (may be "PETSc" or "SuperLU")
  TGETOPTDEF_S_ND_PF(thash,string,local_solver,PETSc);
  }
  if (local_solver_s == string("PETSc")) local_solver = PETSc;
  else if (local_solver_s == string("SuperLU")) local_solver = SuperLU;
  else assert(0);

  //o PETSc parameter related to the efficiency in growing
  //   the factored profile.
  TGETOPTDEF_ND_PF(thash,double,pc_lu_fill,5.);

  //o Print the Schur matrix (don't try this for big problems).
  TGETOPTDEF_ND_PF(thash,int,print_Schur_matrix,0);

  //o Print Finite State Machine transitions
  TGETOPTDEF_ND_PF(thash,int,print_fsm_transition_info,0);

  // Scatter the profile graph
  lgraph.scatter();

  MPI_Comm_rank (comm, &myrank);
  MPI_Comm_size (comm, &size);

  A_LL_other = new DistMatrix(&pf_part,comm);
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

  // these are the vectors for use with the PETSc matrix constructors
  // flag:= will contain `0' for local dof's and `1' for `interface'
  // dof's. 
  // nnz:= eight vectors that contain the d_nnz, o_nnz for
  // dimensioning PETSc matrices
  // nnz[diagonal/off-diagonal][row-block-type][column-block-type]
  // For instance nnz[D][L][L] = d_nn_z for the 'local-local' block. 
  vector<int> nnz[2][2][2],flag0,flag;
  flag0.resize(neq,0);
  flag.resize(neq,0);

  // First, decide which dof's are marked as `interface' and which as `local'. 

  // A dof in processor `k' is marked as interface if it is connected
  // (has a non zero matrix value) with a dof in a processor with a lower index "k'<k"
  for (k=0; k < neqp; k++) {
    // keq:= number of dof
    keq = dofs_proc_v[k];
    ngbrs_v.clear();
    lgraph.set_ngbrs(keq,ngbrs_v);
    qe = ngbrs_v.end();
    for (q=ngbrs_v.begin(); q!=qe; q++) {
      // leq:= number of dof connected to `keq'
      leq = *q;
      proc_l = part.processor(leq);
      // if leq is in other processor, then either `keq' or `leq' are
      // interface. This depends on `leq<keq' or `leq>keq'
      if (proc_l < myrank) {
	// 	printf("[%d] marking node %d\n",keq);
	flag0[keq]=1;
      } else if (proc_l > myrank) {
	flag0[leq]=1;
	// printf("[%d] marking node %d\n",myrank,leq);
      }
    }
  }
  
  // Each processor marks as interface dof's that belongs to them and
  // to other so that, after, we have to combine all them with an
  // Allreduce
  MPI_Allreduce(flag0.begin(), flag.begin(), neq, MPI_INT, 
		MPI_MAX, comm);

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
    if (flag[ dofs_proc_v[k] ]==0) dof2loc[k] = n_loc_pre++;
  }
  loc2dof.resize(n_loc_pre);
  for (dof = 0; dof < neqp; dof++) {
    loc = dof2loc[dof];
    if (loc>=0) loc2dof[loc] = dof;
  }
  
#define INF INT_MAX
  //o The maximum number of vertices in the coarse mesh for
  // sub-partitioning the dof graph in the IISD matrix. 
  TGETOPTDEF_ND_PFMAT(&thash,int,max_partgraph_vertices_proc,INF);
#undef INF
  TGETOPTDEF_ND_PFMAT(&thash,int,iisd_subpart,1);

  local_graph.lgraph = &lgraph;
  local_graph.init(n_loc_pre);
  local_graph.loc2dof = loc2dof.begin();
  local_graph.dof2loc = dof2loc.begin();
  local_graph.dofs_proc = dofs_proc.begin();
  local_graph.proc2glob = &proc2glob;
  local_graph.partit = &part;
  local_graph.myrank = myrank;
  local_graph.flag = flag.begin();

  local_graph.part(max_partgraph_vertices_proc,iisd_subpart);
  // Mark those local dofs that are connected to a local dof in a
  // subdomain with lower index in the subpartitioning as interface.
  for (k=0; k<n_loc_pre; k++) {
    subdoj = local_graph.vrtx_part(k);
    ngbrs_v.clear();
    local_graph.set_ngbrs(k,ngbrs_v);
    qe = ngbrs_v.end();
    for (q=ngbrs_v.begin(); q!=qe; q++) {
      if (local_graph.vrtx_part(*q)<subdoj) {
	flag[dofs_proc_v[ loc2dof[k] ] ] = 1;
	break;
      }
    }
  }
  ngbrs_v.clear();

  local_graph.clear();
  dof2loc.clear();
  loc2dof.clear();

  // We have to combine all them again with an Allreduce
  MPI_Allreduce(flag.begin(), flag0.begin(), neq, MPI_INT, 
		MPI_MAX, comm);
  // recopy on `flag'...
  memcpy(flag.begin(),flag0.begin(),neq*sizeof(int));
  flag0.clear();

  // map:= map[k] is the location of dof `k1+k' in reordering such
  // that the `local' dofs are the first `n_loc' and the `interface'
  // dofs are the last n_int.
  // n_int := number of nodes in this processor that are `interface'
  // n_loc := number of nodes in this processor that are `local'
  map_dof.resize(neq,0);
  n_int_v.resize(size+1,0);
  n_loc_v.resize(size+1,0);

  // number all `loc' dof's
  n_loc_tot = 0;
  n_loc_v[0] = 0;
  // Number first all loc dof's in processor 0, then on processor 1,
  // and so on
  for (rank = 0; rank < size; rank++) {
    // This is no scalable on the number of processors, but is very
    // fast so that it may cause problems only for a large number of
    // processors. 
    for (keq=0; keq<neq; keq++)
      if(part.processor(keq)==rank 
	 && flag[keq] == 0) map_dof[keq] = n_loc_tot++;
    n_loc_v[rank+1] = n_loc_tot;
  }

  // number all `int' dof's
  n_int_tot = n_loc_tot;
  n_int_v[0] = n_loc_tot;
  for (rank = 0; rank < size; rank++) {
    for (keq=0; keq<neq; keq++)
      if(part.processor(keq)==rank &&
	 flag[keq] == 1) map_dof[keq] = n_int_tot++;
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
    keq = dofs_proc_v[k];
    // type of dof (local or interface)
    // row_t = (flag[keq] ? I : L);
    // Index in the local PETSc matrices (maped index)
    row = map_dof[keq];
    row_t = (row < n_loc_tot ? L : I);
    
    // Correct dof's
    row -= (row < n_loc_tot ? n_locp : n_intp);
    // loop over the connected dof's
    ngbrs_v.clear();
    lgraph.set_ngbrs(keq,ngbrs_v);

    qe = ngbrs_v.end();
    for (q=ngbrs_v.begin(); q!=qe; q++) {
      // leq:= number of dof connected to `keq' i.e. `A(keq,leq) != 0' 
      leq = *q;
      // type of dof
      // col_t = (flag[leq] ? I : L);
      col_t = (map_dof[leq] < n_loc_tot ? L : I);
      // diagonal or off-diagonal (in PETSc sense)
      od = ( (part.processor(leq)!= myrank) ? O : D);
      // By construction, the local-local block should be in the
      // diagonal part
      assert(!(od==O && row_t==L && col_t==L));
      // count 
      if (!(row>=0 && row < nnz[od][row_t][col_t].size())) {
	printf("row %d, size: %d\n",row,nnz[od][row_t][col_t].size());
	MPI_Abort(comm,1);
      }
      nnz[od][row_t][col_t][row]++;
    }

  }

  // deallocate profile (graph)
  lgraph.clear();

#if 0
  // Prints d_nnz, o_nnz for block LL, IL, IL and II in turn
  // For each block prints the d_nnz and o_nnz in turn
  for (int row_t=0; row_t<2; row_t++) {
    for (int col_t=0; col_t<2; col_t++) {
      PetscPrintf(comm,"[%s]-[%s] block\n",
		  (row_t==0? "LOC" : "INT"),(col_t==0? "LOC" :
					     "INT"));
      // number of rows in this block in this processor
      nrows=(row_t==L ? n_loc : n_int);
      PetscSynchronizedPrintf(comm,
			      "%d/%d rows on processor [%d]\n",
			      nrows,neqp,myrank);
      for (row=0; row<nrows; row++) {
	d_nz = nnz[D][row_t][col_t][row];
	o_nz = nnz[O][row_t][col_t][row];
	if (d_nz|| o_nz) 
	  PetscSynchronizedPrintf(comm,
				  "row=%d, d_nnz=%d, o_nnz=%d\n",row,
				  d_nz,o_nz);
      }
      PetscSynchronizedFlush(comm); 
    }
  }
  PetscFinalize();
  exit(0);
#endif

  ierr = MatCreateMPIAIJ(comm,n_loc,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][L][I].begin(),
			 PETSC_NULL,nnz[O][L][I].begin(),
			 &A_LI); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating loc-int matrix\n"); 
    
  ierr = MatCreateMPIAIJ(comm,n_int,n_loc,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][L].begin(),
			 PETSC_NULL,nnz[O][I][L].begin(),
			 &A_IL); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating int-loc matrix\n"); 
    
  ierr = MatCreateMPIAIJ(comm,n_int,n_int,
			 PETSC_DETERMINE,PETSC_DETERMINE,
			 PETSC_NULL,nnz[D][I][I].begin(),
			 PETSC_NULL,nnz[O][I][I].begin(),
			 &A_II); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating int-int matrix\n"); 
  
  // extern int mult(Mat,Vec,Vec);
  ierr = MatCreateShell(comm,n_int,n_int,
			PETSC_DETERMINE,PETSC_DETERMINE,this,&A);
  CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating shell matrix\n"); 
  P=A;

  MatShellSetOperation(A,MATOP_MULT,(void *)(&IISD_mult));
  MatShellSetOperation(A,MATOP_MULT_TRANS,(void *)(&IISD_mult_trans));

  ierr = VecCreateMPI(comm,n_loc,PETSC_DETERMINE,&x_loc);
  CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc' vector\n"); 

  ierr = VecCreateSeq(PETSC_COMM_SELF,n_loc,&y_loc_seq);
  CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating `y_loc_seq' vector\n"); 

  ierr = VecDuplicate(y_loc_seq,&x_loc_seq); CHKERRQ(ierr); 
  // PETSCFEM_ASSERT0(ierr==0,"Error creating `x_loc_seq' vector\n"); 

  // Shortcuts
  AA[L][L] = &A_LL;
  AA[L][I] = &A_LI;
  AA[I][L] = &A_IL;
  AA[I][I] = &A_II;

  // Save copy of d_nnz_LL for use later when recreating A_LL
  d_nnz_LL = nnz[D][L][L];

  ierr = clean_mat_a(); CHKERRQ(ierr);

  return 0;
}
