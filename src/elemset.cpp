//__INSERT_LICENSE__
//$Id: elemset.cpp,v 1.92 2004/09/25 09:37:57 mstorti Exp $

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <vector>
#include <set>

#include <src/fem.h>
#include <src/utils.h>
#include <src/getprop.h>
#include <src/elemset.h>
#include <src/idmap.h>
#include <src/dofmap.h>
#include <src/arglist.h>
#include <src/readmesh.h>
#include <src/pfmat.h>
#include <src/timestat.h>
#include <src/util3.h>
#include <src/autostr.h>

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1
extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define MASK(j,k,p,q) VEC4(argd.profile,j,k,ndof,p,nel,q,ndof)
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
int Elemset::real_nodes(int iele,const int *&nodes) {
  nodes = &ICONE(iele,0);
  return nel;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Elemset::download_vector(int nel,int ndof,Dofmap *dofmap,
			     arg_data &argd,
			     int myrank,int el_start,int el_last,
			     int iter_mode,const TimeData *time_data) {

  int iele,iele_here,kloc,node,kdof;
  
  // If the vector passed has a time_data value, then
  // use it as time
  const TimeData *time_d = time_data;
  if (argd.options & USE_TIME_DATA) {
    assert(argd.time_data);
    time_d = argd.time_data;
  }

  double *locst = argd.locst;
  int nen = nel*ndof;
  iele_here=-1;
  for (iele=el_start; iele<=el_last; iele++) {

    if (!compute_this_elem(iele,this,myrank,iter_mode)) continue;
    iele_here++;

#if 0
    // Old slow version
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(iele,kloc);
      for (kdof=1; kdof<=ndof; kdof++) {
	dofmap->get_nodal_value(node,kdof,argd.sstate,argd.ghost_vals,
				time_d,LOCST(iele_here,kloc,kdof-1));
	// printf("setting node %d,kdof %d to %f\n",node,
	// kdof,LOCST(iele,kloc,kdof-1));
      }
    }
#else
    // Fast new version
    int *nodep = icone+iele*nel;
    int *nodep_end = nodep + nel;
    double *w = locst + iele_here * nen;
    while (nodep < nodep_end) {
      for (kdof=1; kdof<=ndof; kdof++) {
 	dofmap->get_nodal_value(*nodep,kdof,argd.sstate,
 				argd.ghost_vals,time_d,*w++);
      }
      nodep++;
    }
#endif
  }
  return 0;
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
  double global_val;
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
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Stat {
public:
  double vmin,vmax,sum;
  int count,initialized;
  Stat() {initialized=0;}
  void reset() {initialized=0;}
  void add(double val) {
    if (!initialized) {
      initialized = 1;
      vmin = val;
      vmax = val;
      count = 1;
      sum = val;
    } else {
      if (val<vmin) vmin = val;
      if (val>vmax) vmax = val;
      sum += val;
      count++;
    }
  }
  double avrg() { 
    assert(initialized);
    return sum/double(count);
  }
  double min() { 
    assert(initialized);
    return vmin;
  }
  double max() { 
    assert(initialized);
    return vmax;
  }
  int n() {
    assert(initialized);
    return count;
  }
  double total() {
    assert(initialized);
    return sum;
  }
  void print_stat(char * s= NULL) {
    if (s) PetscPrintf(PETSC_COMM_WORLD,
		       "Event %s ------------------------\n",s);
    if (initialized) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, 
			      "[%d] total: %g, max: %g, min: "
			      "%g, avrg: %g, count: %d\n",
			      MY_RANK, total(), max(), min(), 
			      avrg(), n());
    } else {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[not initialized]\n");
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::clear_error() { error_code_m=0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::set_error(int error_code_a) { error_code_m = error_code_a; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::check_error() {
  int error;
  int ierr = MPI_Allreduce(&error_code_m,&error,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
  error_code_m = error;
  handle_error(error_code_m);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::handle_error(int error) {
  if (error) PETSCFEM_ERROR("Elemset name \"%s\", ptr %p, set error %d\n",
			    name(),this,error);  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Elemset::error_code() { return error_code_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int assemble(Mesh *,arg_list ,Dofmap *,char *)"
int assemble(Mesh *mesh,arg_list argl,
	     Dofmap *dofmap,const char *jobinfo,const TimeData *time_data) {

#define ARGVJ (arg_data_v[j])

  int iele,nelem,nel,ndof,*icone,ndoft,kloc,kdof,
    myrank,ierr,kdoft,iele_here,k;

  Darray *ghostel;
  Darray *elemsetlist = mesh->elemsetlist;
  Nodedata *nodedata = mesh->nodedata;
  HPChrono hpchrono,hpc2,hpc3,hpcassmbl;
  Stat out_of_loop, in_loop, wait;
  // If `delayed_flush=1' then we compute new element values before
  // doing the flush assembly. This can be more efficient. 
  int delayed_flush = 1;

  // Time statistics
  double total, compt, upload, download,
    tot_s,tot_e,compt_s,bus_e,upl,upl_s,
    assmbly, assmbly_s;
  hpc2.start();
  hpchrono.start();

  //o Debug the process of building the matrix profile. 
  TGETOPTDEF(mesh->global_options,int,debug_compute_prof,0);
  //o Debug the process of building the matrix profile. 
  TGETOPTDEF(mesh->global_options,int,report_assembly_time,0);
  //o Print the local chunk size used for each elemset in each
  // processor for each chunk. 
  TGETOPTDEF(mesh->global_options,int,print_local_chunk_size,0);

  // This is the argument list to be passed to the element routine
  int narg = argl.size();
  arg_data_list arg_data_v(narg);

  // pref:= Local values (reference state for finite difference jacobian).
  double *pref,fdj;

  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  // max weight (processor speed)
  float w_max=dofmap->tpwgts[0]; 
  for (int jproc=1; jproc<dofmap->size; jproc++)
    if (dofmap->tpwgts[jproc] > w_max) w_max = dofmap->tpwgts[jproc];

  int j;
  int iter_mode = NOT_INCLUDE_GHOST_ELEMS;
  // any_fdj:= (boolean) flags if the call includes a "perturbed
  // vector". This implies the calculation of matrices by finite
  // difference approximation to jacobian. 
  // j_pert:= this points to which argument is the vector to be
  // perturbed. 
  int any_fdj = 0,j_pert;  
  // any_include_ghost_elems:= any_not_include_ghost_elems=0:=
  // Flag whether any argument corresponds to that iteration mode. 
  // Iteration modes are mutually exclusive, so that we must check
  // that `!any_include_ghost_elems || !any_not_include_ghost_elems'
  int any_include_ghost_elems=0, any_not_include_ghost_elems=0; 
  for (j=0; j<narg; j++) {
    /// Copy the options to the arg_data_v structure. 
    ARGVJ.options = argl[j].options; 
    ARGVJ.must_flush = 0;
    // PetscPrintf(PETSC_COMM_WORLD,"Argument %d\n",j);
    if (argl[j].options & DOWNLOAD_VECTOR) {
      Vec *x;
      if (argl[j].options & USE_TIME_DATA) {
	State *state = (State *) (argl[j].arg);
	x = (Vec *)&state->v();
	ARGVJ.time_data = (TimeData*) (&state->t());     
      } else {
	x = (Vec *) (argl[j].arg);
      }
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

    if (argl[j].options & UPLOAD_PROFILE) {
      if (argl[j].options & PFMAT) {
	ARGVJ.pfA = (PFMat *) (argl[j].arg);
      } else {
	iter_mode = INCLUDE_GHOST_ELEMS;
	ARGVJ.da= da_create_len(sizeof(Node),dofmap->neq);
	Node nodeq = Node(-1,0);
	for (k=0; k < dofmap->neq; k++) {
	  da_set(ARGVJ.da,k,&nodeq);
	}
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

    /// Set the appropriate flag
    if (iter_mode == INCLUDE_GHOST_ELEMS) {
      any_include_ghost_elems=1;
    } else {
      any_not_include_ghost_elems=1;
    }

  }

  // Check that iteration modes are not mixed
  assert(!any_include_ghost_elems || !any_not_include_ghost_elems);

  for (int ielset=0; ielset<da_length(elemsetlist); ielset++) {
    
    // Initialize time accumulators
    compt = 0.0;
    upload = 0.0;
    download = 0.0;
    assmbly = 0.0;
    tot_s = MPI_Wtime();

    Elemset *elemset  = *(Elemset **)da_ref(elemsetlist,ielset);

    // Skip this elemset if indicated.
    int skip_elemset;
    elemset->ask(jobinfo,skip_elemset);
    if (skip_elemset) continue;

    int nel_here=elemset->nelem_here;
    
    // Initialize chronometer, in order to perform load balance 
    hpc3.start();
 
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
    //o Report consumed time for the elemset. Useful for building
    // the table of weights per processor. 
    TGETOPTDEF(elemset->thash,int,report_consumed_time,0);
    //o Print statistics about time spent in communication and residual evaluation
    TGETOPTDEF(mesh->global_options,int,report_consumed_time_stat,0);

    int local_chunk_size;
    // scaled chunk_size in order to balance processors 
    local_chunk_size = 
      (int)(chunk_size*dofmap->tpwgts[myrank]/w_max) +1;
    if (print_local_chunk_size) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d] type %s, chunk_size %d, local_chunk_size %d\n",
			      MY_RANK,elemset->type,chunk_size,local_chunk_size);
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
    }

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
	ARGVJ.retval  = new double[chunk_size*ndoft*ndoft];
      if (argl[j].options & (ASSEMBLY_MATRIX | UPLOAD_PROFILE)) {
	ARGVJ.profile = new double[ndoft*ndoft];
	for (int kk=0; kk<ndoft*ndoft; kk++) ARGVJ.profile[kk] = 1.;
      }
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

    // here application writers can perform tasks before the
    // chunk loop calling `assemble'
    elemset->before_assemble(arg_data_v,nodedata,dofmap,
			     jobinfo,myrank,el_start,el_last,iter_mode,
			     time_data);
    
    out_of_loop.add(hpchrono.elapsed());

    while (1) { // Loop over chunks. 
      hpchrono.start();
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
	if (argl[j].options & DOWNLOAD_VECTOR) {
	  upl_s = MPI_Wtime();
	  elemset->download_vector(nel,ndof,dofmap,ARGVJ,
				   myrank,el_start,el_last,iter_mode,
				   time_data);
	  download += MPI_Wtime() - upl_s;
	}
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

      
      elemset->clear_error();
      if (iele_here > -1) {
	// if (1) {
	compt_s = MPI_Wtime();
	elemset->assemble(arg_data_v,nodedata,dofmap,
			  jobinfo,myrank,el_start,el_last,iter_mode,
			  time_data);
	compt += MPI_Wtime() - compt_s;
      } else {
	// printf("[%d] not processing because no elements...\n",myrank);
      }
      elemset->check_error();

      // Upload return values
      for (j=0; j<narg; j++) {
	if (report_assembly_time) hpcassmbl.start();
	// Do flush assmbly before to upload new values in matrix
	if ((argl[j].options & ASSEMBLY_MATRIX) 
	    && ARGVJ.must_flush) {
	  assmbly_s = MPI_Wtime();
	  ierr = (ARGVJ.pfA)
	    ->assembly_end(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ARGVJ.must_flush = 0;
	  assmbly += MPI_Wtime() - assmbly_s;
	}
	if (argl[j].options & UPLOAD_RETVAL) { 
	  upl_s = MPI_Wtime();
	  elemset->upload_vector(nel,ndof,dofmap,argl[j].options,ARGVJ,
				 myrank,el_start,el_last,iter_mode);
	  upload += MPI_Wtime() - upl_s;
	}
	if (report_assembly_time) {
	  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				  "[%d] Upload time %f secs.\n",
				  MY_RANK,hpcassmbl.elapsed());
	  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
	}
      }

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

	    //#define DEBUG_ME
#ifdef DEBUG_ME
	    printf("kloc %d, kdof %d\n",kloc,kdof);
#endif
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
		    fdj = -(RETVALT(iele_here,kdoft)-
			    REFREST(iele_here,kdoft))/epsilon_fdj;
#ifdef DEBUG_ME
		    printf("ref, new, jac: %g %g %g\n",
			   REFREST(iele_here,kdoft),
			   RETVALT(iele_here,kdoft),fdj);
#endif
		    RETVALT(iele_here,kdoft) = fdj;
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

      in_loop.add(hpchrono.elapsed());
      hpchrono.start();
  
      for (j=0; j<narg; j++) {
	
	assmbly_s = MPI_Wtime();
	if (argl[j].options & ASSEMBLY_MATRIX) {
	  if (report_assembly_time) hpcassmbl.start();
	  if (argl[j].options & PFMAT) {
	    if (ARGVJ.must_flush) {
	      ierr = (ARGVJ.pfA)
		->assembly_end(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	      ARGVJ.must_flush = 0;
	    }
	    ierr = (ARGVJ.pfA)
	      ->assembly_begin(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ARGVJ.must_flush = 1;
	    if (!delayed_flush && ARGVJ.must_flush) {
	      ierr = (ARGVJ.pfA)
		->assembly_end(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	      ARGVJ.must_flush = 0;
	    }
	  } else {
	    ierr = MatAssemblyBegin(*(ARGVJ.A),
				    MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*(ARGVJ.A),
				  MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  }
	  if (report_assembly_time) {
	    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				    "[%d] Assembly time \"%s\"/\"%s\" %f secs.\n",
				    MY_RANK,elemset->type,jobinfo,hpcassmbl.elapsed());
	    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
	  }
	}
	assmbly += MPI_Wtime() - assmbly_s;
      }

      wait.add(hpchrono.elapsed());
      hpchrono.start();

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
      in_loop.add(hpchrono.elapsed());
      hpchrono.start();
      if (global_has_finished) break;
    } // end loop over chunks

    elemset->after_assemble(jobinfo);

    if (report_consumed_time) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Performance report elemset \"%s\" task \"%s\"\n"
		  "[proc] - elems - compt[sec](rate[sec/Ke]) - upl/dwl[sec]"
		  " - assmbly[sec]  - other[sec]\n",
		  elemset->type,jobinfo);

      double rate=1000.0*compt/elemset->nelem_here;
      total = MPI_Wtime() - tot_s;
      double upd = upload+download;
      double other = total-(compt+upd+assmbly);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[%d]  %7d  %7.2g/%3.0f%%(%7.2g)  "
			      "%7.2g/%3.0f%%  %7.2g/%3.0f%%  %7.2g/%3.0f%%\n",
			      myrank,elemset->nelem_here,
			      compt,100.0*compt/total,rate,
			      upd,100.0*upd/total,
			      assmbly,100.0*assmbly/total,
			      other,100.0*other/total);
#if 0
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "upload %10.3gsecs, (%10.3g secs/Ke)\n",
			      upload,1000.0*upload/elemset->nelem_here);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "download %10.3gsecs, (%10.3g secs/Ke)\n",
			      download,1000.0*download/elemset->nelem_here);
#endif
      PetscSynchronizedFlush(PETSC_COMM_WORLD);

      PetscPrintf(PETSC_COMM_WORLD,"Total element compt. %10.3gsecs\n",total);
    }

    // To be done for each elemset
    if (any_fdj) delete[] pref;
    for (j=0; j<narg; j++) {
      if (argl[j].options & ASSEMBLY_MATRIX
	  && (argl[j].options & PFMAT)
	  && ARGVJ.must_flush) {
	    ierr = (ARGVJ.pfA)
	      ->assembly_end(MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ARGVJ.must_flush = 0;
      }
      if (argl[j].options & DOWNLOAD_VECTOR) {
	delete[] ARGVJ.locst;
      }

      if (argl[j].options & DELETE_RETVAL) 
	delete[] ARGVJ.retval;

      if (argl[j].options & IS_FDJ) {
	delete[] ARGVJ.retval;
	delete[] ARGVJ.refres;
      }
      delete[] ARGVJ.profile;
      ARGVJ.profile = NULL;
    }

    if (report_consumed_time_stat) {
      out_of_loop.add(hpchrono.elapsed());
      out_of_loop.print_stat("Out of loop");
      in_loop.print_stat("In loop");
      wait.print_stat("Wait");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] Total in assemble: %g\n",
			      MY_RANK,hpc2.elapsed());
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
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
	// ((PFMat *)(argl[j].arg))->create(ARGVJ.da,dofmap,debug_compute_prof);
	((PFMat *)(argl[j].arg))->create();
      } else {
        ierr = compute_prof(ARGVJ.da,dofmap,
			    myrank,(Mat *)(argl[j].arg),
			    debug_compute_prof);
	// fixme:= should this go here also ???
	// ierr =  MatSetOption((Mat *)(argl[j].arg),
	// MAT_NEW_NONZERO_ALLOCATION_ERR);

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
int Elemset::assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
		      const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		      const TimeData *time_data) {
  printf("assemble: not known Elemset\n"); exit(1);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Elemset::weight() {
  int ierr;
  //o Element weight for the processor
  TGETOPTDEF(thash,int,element_weight,1);
  return element_weight;
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
const char * Elemset::name() {
//    static string name_r("__ANONYMOUS__");
//    ::get_string(thash,"name",name_r,1);
  return name_m.c_str();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void *& Elemset::local_store_address(int k) {
  int kk = epart2[k];
  assert(e1 <= kk < e2);
  return local_store[kk-e1];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double NewElemset::prop_val(ElementIterator &element,
			    Property &prop,double t) const {
  const double *val_p = prop_array(element,prop);
  if (prop.eval_fun) 
    return prop.eval_fun(t,(val_p ? *val_p : 0.),prop.fun_data);
  else if (val_p) return *val_p;
  else return 0.;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double NewElemset::prop_val(ElementIterator &element,
			    Property &prop) const {
  return *prop_array(element,prop);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const double *NewElemset::prop_array(ElementIterator &element,
				     Property &prop) const {
  if (prop.indx<0) {
    return prop.ptr;
  } else {
    return element_props(element) + prop.indx;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ \
   "int NewElemset::get_vec_double(const char *,vector<double> &,int=0) const"
int NewElemset::get_vec_double(const char *name,
		   vector<double> &retval,int defval) const {

  const char *value;
  if (!defval) retval.clear();
  get_entry(name,value);
  if (!defval & !value) return 1;// Either we have a default value or
				// the user enters an entry
  if (!value ) return 0;
  read_double_array(retval,value);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Mesh::Mesh() : elemsetlist(NULL), nodedata(NULL), 
  global_options(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Elemset *Mesh::find(const string &name) {
  int nelemsets = da_length(elemsetlist);
  Elemset *vol_elem = NULL;
  for (int k=0; k<nelemsets; k++) {
    Elemset *e = *(Elemset **) da_ref(elemsetlist,k);
    if (!strcmp(e->name(),name.c_str())) {
      vol_elem = e;
      break;
    }
  }
  return vol_elem;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Mesh::~Mesh() {
  if (nodedata) delete nodedata;
  nodedata=NULL;
  for (int j=0; j<da_length(elemsetlist); j++) {
    Elemset *elemset  = *(Elemset **)da_ref(elemsetlist,j);
    delete elemset;
  }
  da_destroy(elemsetlist);
  elemsetlist=NULL;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
map<string,Elemset *> Elemset::elemset_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
string Elemset::anon("__ANONYMOUS__");

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::register_name(const string &name_a,const char *type) {
  // Assign a name to the elemset.If name has tnot been entered
  // then try with `type', `type_0', `type_1' and so on, until a
  // unique name is found. In order to avoid infinite loops, assign
  // finally `type_ptr' where `ptr' is th pointer to the elemset. 
  AutoString ename;
  name_m = name_a;
  if (name_m == anon) {
#define MAX_ELEMSET_SFX 1000
    int j;
    for (j=-1; j<MAX_ELEMSET_SFX; j++) {
      // int Nbuf = asprintf(&ename,"%s_%d",type,j);
      // assert(Nbuf>=0);
      if (j<0) ename.set(type);
      else ename.sprintf("%s_%d",type,j);

      if (elemset_table.find(ename.str())
	  ==elemset_table.end()) break;
    }
    if (j==MAX_ELEMSET_SFX) ename.sprintf("%s_%p",type,this);
    name_m = local_copy(ename.str());
    PETSCFEM_ASSERT0(j!=MAX_ELEMSET_SFX,
		     "Couldn't generate automatic name for this  elemset!!\n");
  }
  elemset_table[name_m] = this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Elemset::size() { return nelem; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::read(FileStack *fstack) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Elemset::Elemset() : type(NULL), icone(NULL), elemprops(NULL),
		     elemiprops(NULL), elemprops_add(NULL),
		     elemiprops_add(NULL), thash(NULL),
                     elem_conne(NULL), error_code_m(0) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Elemset::~Elemset() {
  DELETE_VCTR(type);
  DELETE_VCTR(icone);
  DELETE_VCTR(elemprops);
  DELETE_VCTR(elemiprops);
  DELETE_VCTR(elemprops_add);
  DELETE_VCTR(elemiprops_add);
  DELETE_SCLR(thash);
  DELETE_VCTR(elem_conne);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Nodedata::Nodedata() : nodedata(NULL) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Nodedata::~Nodedata() { DELETE_VCTR(nodedata); }
