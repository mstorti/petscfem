//__INSERT_LICENSE__
//$Id: elemset.cpp,v 1.26 2001/10/15 01:31:02 mstorti Exp $

#include "fem.h"
#include <vector>
#include <set>
#include "utils.h"
#include "getprop.h"
#include "elemset.h"
#include "idmap.h"
#include "dofmap.h"
#include "arglist.h"
#include "readmesh.h"
#include "pfmat.h"

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ICONE(j,k) VEC2(icone,j,k,nel)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "compute_this_elem" 
int compute_this_elem(const int & iele,const Elemset *elemset,const int & myrank,
		      int iter_mode) {
  int flag;
  flag = (elemset->epart[iele] == myrank+1);
  if (flag) return 1;
  if (iter_mode == INCLUDE_GHOST_ELEMS) {
    int is_ghost_elem = da_bsearch(elemset->ghost_elems,&iele,int_cmp,NULL);
    return (is_ghost_elem>=0);
  } else {
    return 0;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int_cmp"
int int_cmp (const void *left,const void *right, void *args) {
  int l = *(int *) left;
  int r = *(int *) right;
  if (l > r) return +1;
  if (l < r) return -1;
  return 0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void node_insert(Darray *da,int j,int k) {
  Node *nodep,nodeq;
  int posit=j,newpos;

  //  printf(" Inserting %d %d, pasing by (",j,k);
  while (1) {
    nodep = (Node *)da_ref(da,posit);
    if (nodep->next == -1) {
      nodeq = Node(-1,0);
      newpos = da_append(da,&nodeq);
      //      printf("). Appended at position %d\n",posit);
      nodeq = Node(newpos,k);
      da_set(da,posit,&nodeq);
      break;
    } else if (nodep->val==k) {
      //      printf("). Found at position %d\n",posit);
      break;
    } else {
      //      printf(" %d",nodep->val);
      posit = nodep->next;
    }
  }
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Elemset::download_vector(int nel,int ndof,Dofmap *dofmap,
			     arg_data &argd,
			     int myrank,int el_start,int el_last,
			     int iter_mode,const TimeData *time_data=NULL) {

  int iele,iele_here,kloc,node,kdof,locdof;
  
  double *locst = argd.locst;
  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {

    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;

    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=1; kdof<=ndof; kdof++) {
	dofmap->get_nodal_value(node,kdof,argd.sstate,argd.ghost_vals,
				time_data,
				LOCST(iele_here,kloc,kdof-1));
	// printf("setting node %d,kdof %d to %f\n",node,
	// kdof,LOCST(iele,kloc,kdof-1));
      }
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "upload_vector"
int Elemset::upload_vector(int nel,int ndof,Dofmap *dofmap,
		  int options,arg_data &argd,int myrank,
		  int el_start,int el_last,int iter_mode,
		  int klocc=0,int kdofc=0) {

  int iele,kloc,node,kdof,locdof,lloc,nodel,ldof,locdofl,ierr,
    load_vec,load_mat,load_mat_col,comp_prof,comp_mat,iele_here,
    pfmat;

  double *retval = argd.retval;

  int neq;
  //  sp_entry *spe, *spel;
  double q,ql,val,vall;
  //  ierr = VecGetSize(vec,&neq); CHKERRQ(ierr); 
  neq = dofmap->neq;
  row_t::iterator entry,entryc;
  
  // static Darray *row, *rowc;
  row_t row,rowc;
  IdMapRow row_v,rowc_v;
  IdMapEntry *entry_v,*entryc_v;

  load_vec = (options & (UPLOAD_VECTOR | UPLOAD_VECTOR_LOCST));
  load_mat = (options & (UPLOAD_MATRIX | UPLOAD_PROFILE));
  load_mat_col = (options & IS_FDJ);
  pfmat = (options & PFMAT);

  // In order to compute profiles (comp_prof==1) we have to do all the
  // work in all th processors
  comp_prof= (options & UPLOAD_PROFILE); // to be defined later

  InsertMode mode = ADD_VALUES;
  if (options & UPLOAD_VECTOR_LOCST) mode = INSERT_VALUES;

  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {
    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=0; kdof<ndof; kdof++) {
	dofmap->get_row(node,kdof+1,row_v);

	for (int ientry=0; ientry<row_v.size(); ientry++) {
	  entry_v = &row_v[ientry];
	  locdof = entry_v->j;
	  if (locdof>neq) continue; // ony load free nodes
	  val = (entry_v->coef) * RETVAL(iele_here,kloc,kdof);

	  // load local values on global vector
	  if (load_vec) {
	    VecSetValue(*(argd.x),locdof-1,val,mode);
	  }

	  // load finite difference jacobian computed by differences
	  if (load_mat_col) {
	    nodel = ICONE(iele,klocc);
	    vall = RETVAL(iele_here,kloc,kdof);
	    dofmap->get_row(nodel,kdofc+1,rowc_v);

	    for (int ientryc=0; ientryc<rowc_v.size(); ientryc++) {
	      entryc_v = &rowc_v[ientryc];
	      locdofl = entryc_v->j;
	      if (locdofl>neq) continue;
	      val = entry_v->coef * entryc_v->coef * vall;

	      int kd=locdof-1,kdl=locdofl-1;
	      if (val != 0) {
		if (comp_prof) {
		  node_insert(argd.da,kd,kdl);
		  node_insert(argd.da,kdl,kd); // be sure that profile
				// is symmetric
		} else {
		  if (pfmat) {
		    argd.pfA->set_value(kd,kdl,val,ADD_VALUES); 
		  } else {
		    MatSetValue(*argd.A,kd,kdl,val,ADD_VALUES); 
		  }
		}
	      }
	    }
	  }

	  // load local values on global matrix
	  if (!load_mat) continue;
	  for (lloc=0; lloc<nel; lloc++) {
	    nodel = ICONE(iele,lloc);
	    for (ldof=0; ldof<ndof; ldof++) {
	      vall = RETVALMAT(iele_here,kloc,kdof,lloc,ldof);
	      dofmap->get_row(nodel,ldof+1,rowc_v);

	      for (int ientryc=0; ientryc<rowc_v.size(); ientryc++) {
		entryc_v = &rowc_v[ientryc];
		locdofl = entryc_v->j;
		if (locdofl>neq) continue; // only load for free dof's
		
		val = (entry_v->coef) * (entryc_v->coef) * vall;
		if (val != 0) {
		  if (comp_prof) {
		    int kd=locdof-1,kdl=locdofl-1;
		    // be sure that profile is symmetric
		    node_insert(argd.da,locdof-1,locdofl-1);
		    node_insert(argd.da,locdofl-1,locdof-1); 
		  } else {
		    // printf("(%d,%d) -> %f\n",locdof,locdofl,val);
		    if (pfmat) {
		      argd.pfA->set_value(locdof-1,locdofl-1,val,ADD_VALUES); 
		    } else {
		      MatSetValue(*argd.A,locdof-1,locdofl-1,val,ADD_VALUES);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void vector_assoc_gather(vector<int> *,int,int)" 
int vector_assoc_gather(vector<double> *vector_assoc,
			 int options,int myrank,int size) {
  // each processor passes a vector of doubles, and we have to
  // take the MAX, MIN or ADD opertions over all the processors. 
  // We do this by defining a PETSc vector and performing the
  // appropriate operation. 
  if (size==1) return 0;
  Vec one_elem_per_proc;
  int ierr = VecCreateMPI(PETSC_COMM_WORLD,1,size,&one_elem_per_proc);
  double global_val;int jj;
  for (int j=0; j<vector_assoc->size(); j++) {
    // set local value
    // printf("On proc [%d] j=%d local val %f\n",myrank,j,(*vector_assoc)[j]);
    VecSetValue(one_elem_per_proc,myrank,(*vector_assoc)[j],INSERT_VALUES);
    ierr = VecAssemblyBegin(one_elem_per_proc); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(one_elem_per_proc); CHKERRQ(ierr);
    // perform appropriate function
    if (options & VECTOR_ADD) {
      ierr = VecSum(one_elem_per_proc,&global_val);
    } else if (options & VECTOR_MAX) {
      ierr = VecMax(one_elem_per_proc,NULL,&global_val);
    } else if (options & VECTOR_MIN) {
      ierr = VecMin(one_elem_per_proc,NULL,&global_val);
    } else {
      PFEMERRQ("Internal error: invalid option.\n");
    }
    (*vector_assoc)[j]=global_val;
    // printf("On proc [%d] j=%d gathered val %f\n",
    // myrank,j,(*vector_assoc)[j]);
  }
  ierr = VecDestroy(one_elem_per_proc);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int assemble(Mesh *,arg_list ,Dofmap *,char *)"
int assemble(Mesh *mesh,arg_list argl,
	     Dofmap *dofmap,const char *jobinfo,const TimeData *time_data=NULL) {

#define ARGVJ (arg_data_v[j])

  int iele,nelem,nel,ndof,*icone,ndoft,kloc,node,kdof,locdof,
    lloc,ldof,nodel,locdofl,kk,myrank,ierr,kdoft,iele_here,k;
  int rvsize;

  Darray *ghostel;
  Darray *elemsetlist = mesh->elemsetlist;
  Nodedata *nodedata = mesh->nodedata;
  Chrono chrono;

  //o Debug the process of building the matrix profile. 
  TGETOPTDEF(mesh->global_options,int,debug_compute_prof,0);

  // This is the argument list to be passed to the element routine
  int narg = argl.size();
  arg_data_list arg_data_v(narg);

  // pref:= Local values (reference state for finite difference jacobian).
  double *pref;

  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  // max weight (processor speed)
  float w_max=dofmap->tpwgts[0]; 
  for (int jproc=1; jproc<dofmap->size; jproc++)
    if (dofmap->tpwgts[jproc] > w_max) w_max = dofmap->tpwgts[jproc];

  int nghost_dofs = dofmap->ghost_dofs->size();

  int j;
  int iter_mode = NOT_INCLUDE_GHOST_ELEMS;
  // any_fdj:= (boolean) flags if the call includes a "perturbed
  // vector". This implies the calaculation of matrices by finite
  // difference approximation to jacobian. 
  // j_pert:= this points to which argument is the vector to be
  // perturbed. 
  int any_fdj = 0,j_pert;
  for (j=0; j<narg; j++) {
    /// Copy the options to the arg_data_v structure. 
    ARGVJ.options = argl[j].options; 
    // PetscPrintf(PETSC_COMM_WORLD,"Argument %d\n",j);
    if (argl[j].options & DOWNLOAD_VECTOR) {

      Vec *x = (Vec *) (argl[j].arg);
      Vec *ghost_vec = new Vec;
      ierr = dofmap->create_MPI_ghost_vector(*ghost_vec); CHKERRQ(ierr);
      ARGVJ.x = x;
      ARGVJ.ghost_vec = ghost_vec;
      
      ierr = VecScatterBegin(*x,*ghost_vec,INSERT_VALUES,
			     SCATTER_FORWARD,*dofmap->ghost_scatter); CHKERRA(ierr); 
      ierr = VecScatterEnd(*x,*ghost_vec,INSERT_VALUES,
			   SCATTER_FORWARD,*dofmap->ghost_scatter); CHKERRA(ierr); 
      ierr = VecGetArray(*ghost_vec,
			 &(ARGVJ.ghost_vals)); CHKERRQ(ierr);
      ierr = VecGetArray(*x,&(ARGVJ.sstate)); CHKERRQ(ierr);

    } 

    if (argl[j].options & UPLOAD_VECTOR) 
      ARGVJ.x = (Vec *) (argl[j].arg);

    if (argl[j].options & (UPLOAD_MATRIX | IS_FDJ_MATRIX)) {
      if (argl[j].options & PFMAT) {
	ARGVJ.pfA = (PFMat *) (argl[j].arg);
      } else {
	ARGVJ.A = (Mat *) (argl[j].arg);
      }
    }

    if (argl[j].options & UPLOAD_PROFILE ) {
      iter_mode = INCLUDE_GHOST_ELEMS;
      ARGVJ.da= da_create_len(sizeof(Node),dofmap->neq);
      Node nodeq = Node(-1,0);
      for (k=0; k < dofmap->neq; k++) {
	da_set(ARGVJ.da,k,&nodeq);
      }
    }

    if (argl[j].options & IS_PERT_VECTOR ) {
      if (any_fdj) PFEMERRQ("More than one perturbed vector does not make sense!!");
      any_fdj = 1;
      j_pert=j;
    }

    if (argl[j].options & VECTOR_ASSOC ) {
      arg_data_v[j].vector_assoc = (vector<double> *)(argl[j].arg);
      arg_data_v[j].was_set = 0;
    }

    if (argl[j].options & USER_DATA ) 
      arg_data_v[j].user_data = (void *)(argl[j].arg);
  }

  for (int ielset=0; ielset<da_length(elemsetlist); ielset++) {
    
    Elemset *elemset  = *(Elemset **)da_ref(elemsetlist,ielset);

    // Skip this elemset if indicated.
    int skip_elemset;
    elemset->ask(jobinfo,skip_elemset);
    if (skip_elemset) continue;

    int nel_here=elemset->nelem_here;
    
    // Initialize chronometer, in order to perform load balance 
    chrono.start();
 
    // Put the state in a locst per element and pass
    // to the elemset routine. 

    icone = elemset->icone;
    nelem = elemset->nelem;
    nel   = elemset->nel;
    ndof  = elemset->ndof;
    ghostel = elemset->ghost_elems;

    // ndoft:= size to be passed per element
    ndoft = nel*ndof;
    
    //o Chunk size for the elemset. 
    TGETOPTDEF(elemset->thash,int,chunk_size,ELEM_CHUNK_SIZE);
    //o The increment in the variables in order to
    // compute the finite difference approximation to the
    // Jacobian. Should be order epsilon=sqrt(precision)*(typical
    // magnitude of the variable). Normally, precision=1e-15 so
    // that epsilon=1e-7*(typical magnitude of the
    // variable)
    TGETOPTDEF(elemset->thash,double,epsilon_fdj,EPSILON_FDJ);

    // int chunk_size = ELEM_CHUNK_SIZE;
    // ierr = get_int(elemset->thash,"chunk_size",&chunk_size,1);

    int local_chunk_size;
    // scaled chunk_size in order to balance processors 
    local_chunk_size = (int)(chunk_size*dofmap->tpwgts[myrank]/w_max) +1 ;

    CHKERRQ(ierr);
    // This fixes the 'bug100' bug: If the number of ghost_elements is
    // lower than the numberof local elements in this elemset, then we
    // have to allocate, at less room for the ghost elements (if we
    // have to iterate on the ghost_elems also). 
    int max_chunk_size = nel_here;
    if (iter_mode == INCLUDE_GHOST_ELEMS) 
      // I'm not clear about this, but using the `max' version was
      // wrong because in some cases it performed two chunks when only
      // needed one. 
      // old version:
      // max_chunk_size = maxi(2,max_chunk_size,da_length(ghostel));
      // new version:
      max_chunk_size = max_chunk_size + da_length(ghostel);
    chunk_size = mini(2,local_chunk_size,max_chunk_size);
      
    for (j=0; j<narg; j++) {
      if (argl[j].options & DOWNLOAD_VECTOR) 
	ARGVJ.locst = new double[chunk_size*ndoft];
      if (argl[j].options & UPLOAD_VECTOR_LOCST ) {
	ARGVJ.locst = new double[chunk_size*ndoft];
	ARGVJ.retval = ARGVJ.locst;
      }
      if (argl[j].options & UPLOAD_VECTOR) 
	ARGVJ.retval = new double[chunk_size*ndoft];
      if (argl[j].options & (IS_FDJ)) {
	ARGVJ.retval = new double[chunk_size*ndoft];
	ARGVJ.refres = new double[chunk_size*ndoft];
      }
      if (argl[j].options & ALLOC_MATRIX) 
	ARGVJ.retval = new double[chunk_size*ndoft*ndoft];
    }

    if (any_fdj) pref = new double[chunk_size*ndoft];

    int el_start = 0, chunk = 0, el_last, last_chunk=0;

    //#define DEBUG_CHUNK_PROCESSING
#ifdef DEBUG_CHUNK_PROCESSING
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			    "[%d] %d chunks (approx.)\n"
			    ,myrank,(nel_here-1)/local_chunk_size+1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

    while (1) { // Loop over chunks. 
      chunk++; // number of the current chunk

      iele_here = -1;
      for (iele=el_start; iele<nelem; iele++) {
	if (!compute_this_elem(iele,elemset,myrank,iter_mode)) continue;
	iele_here++;
	if (iele_here == chunk_size-1) break;
      }

#ifdef DEBUG_CHUNK_PROCESSING
      PetscPrintf(PETSC_COMM_WORLD,"chunk %d\n",chunk);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d]     elements in chunk: %d\n",
			      myrank,
			      iele_here+1);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

      el_last = iele;
      if (el_last >=nelem) el_last = nelem-1;
      if (el_last==nelem-1) last_chunk=1;
      // printf("[%d] jobinfo %s, chunk %d, chunk_size %d, here %d,range %d-%d\n",
      // myrank,jobinfo,chunk,chunk_size,iele_here+1,el_start,el_last);

      for (j=0; j<narg; j++) {
	if (argl[j].options & DOWNLOAD_VECTOR) 
	  elemset->download_vector(nel,ndof,dofmap,ARGVJ,
				   myrank,el_start,el_last,iter_mode,
				   time_data);
      }
      
#if 0
      if (!strcmp(jobinfo,"comp_res")) {
	Chrono chrono;
	chrono.start();
	int N=10;
	for (int jjj=0; jjj<N; jjj++) {
	  printf("[loop iter %d]\n",jjj+1);
	  elemset->assemble(arg_data_v,nodedata,dofmap,
			    jobinfo,myrank,el_start,el_last,iter_mode,
			    time_data);
	}
	double cpu = chrono.elapsed();
	printf("total %f, ntimes %d, nelems %d, rate %f [sec/1000/elems/iter]\n",
	       cpu, N, iele_here+1, cpu/(N*(iele_here+1))*1000);
	PetscFinalize();
	exit(0);
      }
#endif 

      if (iele_here > -1) {
	// if (1) {
	elemset->assemble(arg_data_v,nodedata,dofmap,
			  jobinfo,myrank,el_start,el_last,iter_mode,
			  time_data);
      } else {
	// printf("[%d] not processing because no elements...\n",myrank);
      }

      // Upload return values
      for (j=0; j<narg; j++) 
	if (argl[j].options & UPLOAD_RETVAL) 
	  elemset->upload_vector(nel,ndof,dofmap,argl[j].options,ARGVJ,
				 myrank,el_start,el_last,iter_mode);

      // compute columns of jacobian matrices by perturbing each
      // local degree of freedom
      if (any_fdj) {

	// copy to reference state
	memcpy (pref,arg_data_v[j_pert].locst,
		sizeof(double)*chunk_size*ndoft);
	for (j=0; j<narg; j++) 
	  if (argl[j].options & IS_FDJ) 
	    memcpy (ARGVJ.refres,ARGVJ.retval,
		    sizeof(double)*chunk_size*ndoft);

	//	double epsilon=EPSILON_FDJ;
	for (kloc=0; kloc<nel; kloc++) {
	  for (kdof=0; kdof<ndof; kdof++) {

	    // copy reference state on perturbed  state
	    memcpy (arg_data_v[j_pert].locst,pref,
		    sizeof(double)*chunk_size*ndoft);
	  
#define PSTAT(iele,j,k) VEC3(arg_data_v[j_pert].locst,iele,j,nel,k,ndof)
	    // perturb state 
	    iele_here=-1;
	    for (iele=el_start; iele <= el_last; iele++) {
	      if (!compute_this_elem(iele,elemset,myrank,iter_mode)) continue;
	      iele_here++;
	      PSTAT(iele_here,kloc,kdof) += epsilon_fdj;
	    }

	    // compute perturbed residual
	    elemset->assemble(arg_data_v,nodedata,dofmap,
			      jobinfo,myrank,el_start,el_last,iter_mode,
			      time_data);

#define RETVALT(iele,jk) VEC2(ARGVJ.retval,iele,jk,ndoft)
#define REFREST(iele,jk) VEC2(ARGVJ.refres,iele,jk,ndoft)
	    // make finite difference aprox. to the jacobian
	    for (j=0; j<narg; j++) {
	      if (argl[j].options & IS_FDJ) {
		iele_here=-1;
		for (int iele=el_start; iele <= el_last;iele++) {
		  if (!compute_this_elem(iele,elemset,myrank,iter_mode)) continue;
		  iele_here++;
		  for (kdoft=0; kdoft<ndoft; kdoft++) {
		    RETVALT(iele_here,kdoft) =
		      -(RETVALT(iele_here,kdoft)-REFREST(iele_here,kdoft))/epsilon_fdj;
		  }
		}

		// load on Petsc matrix
		elemset->upload_vector(nel,ndof,dofmap,argl[j].options,ARGVJ,
				       myrank,el_start,el_last,iter_mode,
				       kloc,kdof);

	      }
	    }
	  }
	}
      }

      for (j=0; j<narg; j++) {
	if (argl[j].options & ASSEMBLY_MATRIX) {
	  if (argl[j].options & PFMAT) {
	    ierr = (ARGVJ.pfA)
	      ->assembly_begin(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = (ARGVJ.pfA)
	      ->assembly_end(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  } else {
	    ierr = MatAssemblyBegin(*(ARGVJ.A),
				    MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*(ARGVJ.A),
				  MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  }
	}
      }
      // Has finished processing chunks this processor?
      int local_has_finished= (el_last==nelem-1);

      // Globally has finished to process chunks if
      // all processors have finished
      int global_has_finished;
      ierr = MPI_Allreduce((void *)&local_has_finished,
			(void *)&global_has_finished,1,MPI_INT,
			MPI_LAND,PETSC_COMM_WORLD);

#ifdef DEBUG_CHUNK_PROCESSING
      PetscPrintf(PETSC_COMM_WORLD,"chunk %d\n",chunk);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d] local finished:  %s\n"
			      "     global finished: %s\n",
			      myrank,
			      (local_has_finished ? "YES" : "NO"),
			      (global_has_finished ? "YES" : "NO")),
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

      el_start = el_last+1;
      if (global_has_finished) break;
    } // end loop over chunks

    //o Report consumed time for the elemset. Useful for building
    // the table of weights per processor. 
    TGETOPTDEF(elemset->thash,int,report_consumed_time,0);
    if (report_consumed_time) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Performance report for elemset \"%s\" task \"%s\"\n"
		  "[proc] - total[sec] - rate[sec/Kelement]\n",
		  elemset->type,jobinfo);
      double elapsed;
      elapsed=chrono.elapsed();
      double rate=1000.0*elapsed/elemset->nelem_here;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[proc %d]   %g   %g\n",myrank,elapsed,rate);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }

    // To be done for each elemset
    if (any_fdj) delete[] pref;
    for (j=0; j<narg; j++) {
      if (argl[j].options & DOWNLOAD_VECTOR) {
	delete[] ARGVJ.locst;
      }

      if (argl[j].options & DELETE_RETVAL) 
	delete[] ARGVJ.retval;

      if (argl[j].options & IS_FDJ) {
	delete[] ARGVJ.retval;
	delete[] ARGVJ.refres;
      }
    }
  }

  // To be done after processing all elemsets
  for (j=0; j<narg; j++) {
    if (argl[j].options & DOWNLOAD_VECTOR) {
      ierr = VecRestoreArray(*(ARGVJ.ghost_vec),
			     &(ARGVJ.ghost_vals)); CHKERRQ(ierr); 
      ierr = VecDestroy(*(ARGVJ.ghost_vec));
      ierr = VecRestoreArray(*(ARGVJ.x),
			     &(ARGVJ.sstate)); CHKERRQ(ierr); 
      delete ARGVJ.ghost_vec;
    }

    if (argl[j].options & ASSEMBLY_MATRIX) {
      if (argl[j].options & PFMAT) {
	ierr = (ARGVJ.pfA)->assembly_begin(MAT_FINAL_ASSEMBLY); 
	CHKERRQ(ierr);
	ierr = (ARGVJ.pfA)->assembly_end(MAT_FINAL_ASSEMBLY); 
	CHKERRQ(ierr);
      } else {
	ierr = MatAssemblyBegin(*(ARGVJ.A),
				MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*(ARGVJ.A),
			      MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      }
    }

    if (argl[j].options & ASSEMBLY_VECTOR) {
      ierr = VecAssemblyBegin(*(ARGVJ.x)); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(*(ARGVJ.x)); CHKERRQ(ierr);
    }

    if (argl[j].options & UPLOAD_PROFILE) {
      if (argl[j].options & PFMAT) {
	((PFMat *)(argl[j].arg))->create(ARGVJ.da,dofmap,debug_compute_prof);
      } else {
        ierr = compute_prof(ARGVJ.da,dofmap,
			    myrank,(Mat *)(argl[j].arg),
			    debug_compute_prof);
      }
    }

    if (argl[j].options & VECTOR_ASSOC ) {
      // Perform MIN, MAX and ADD gather operations. 
      ierr = vector_assoc_gather(arg_data_v[j].vector_assoc,
				 argl[j].options,myrank,dofmap->size); CHKERRQ(ierr);
    }
  }

  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::print() {
  printf("========================\n");
  SHVS(nelem,d);
  SHVS(nel,d);
  SHVS(ndof,d);
  thash->print();
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
const char * Elemset::name() {
  const char *name_r,*name_d="__ANONYMOUS__";
  thash->get_entry("name",name_r);
  if (name_r) {
    return name_r;
  } else {
    return name_d;
  }
}

void NewElemset::get_prop(Property &prop,const char *prop_name,int n=1) const {
  // looks int the properties-per-element table
  props_hash_entry *phe = (props_hash_entry *)
    g_hash_table_lookup(elem_prop_names,(void *)prop_name);
  if(phe!=NULL) {
  // If the name is found in the per element properties table
  // (the line props in the text_hash_table) then we store the
  // position in the table and the value will be loaded with
  // 'load_props' for each element
    prop.indx = phe->position;
    prop.length = phe->width;
  } else {
    // If it was not found in the per element properties table then it
    // should be found as a general property in the text_hash_table. 
    int ierr = get_vec_double(prop_name,prop.val,0);
    if (ierr) {
      prop.ptr = NULL;
      prop.length = 0;
    } else {
      prop.ptr = prop.val.begin();
      prop.length = prop.val.size();
    }
  }
}

double NewElemset::prop_val(ElementIterator &element,Property &prop) const {
  return *prop_array(element,prop);
}

const double *NewElemset::prop_array(ElementIterator &element,
				     Property &prop) const {
  if (prop.indx<0) {
    return prop.ptr;
  } else {
    return (element_props(element)+prop.indx);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
   "int NewElemset::get_vec_double(const char *,vector<double> &,int=0) const"
int NewElemset::get_vec_double(const char *name,
		   vector<double> &retval,int defval=0) const {

  const char *value;
  if (!defval) retval.clear();
  get_entry(name,value);
  if (!defval & !value) return 1;// Either we have a default value or
				// the user enters an entry
  if (!value ) return 0;
  read_double_array(retval,value);
  return 0;
}
