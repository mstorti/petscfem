/*__INSERT_LICENSE__*/
//$Id$

#include "fem.h"
#include <set>
#include "utils.h"
#include "getprop.h"
#include "elemset.h"
#include "idmap.h"
#include "dofmap.h"
#include "arglist.h"

// iteration modes
#define NOT_INCLUDE_GHOST_ELEMS 0
#define INCLUDE_GHOST_ELEMS 1

extern int debug_compute_prof;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
//#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
//#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)
#define ARGVJ (arg_data_v[j])

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int assemble(Mesh *,arg_list ,Dofmap *,char *)"
int assemble(Mesh *mesh,arg_list argl,
	     Dofmap *dofmap,char *jobinfo,const void *time_data=NULL) {

  int iele,nelem,nel,ndof,*icone,ndoft,kloc,node,kdof,locdof,
    lloc,ldof,nodel,locdofl,kk,myrank,ierr,kdoft,iele_here,k;
  int rvsize;

  Darray *ghostel;
  Darray *elemsetlist = mesh->elemsetlist;
  Nodedata *nodedata = mesh->nodedata;
  Chrono chrono;
  
  // This is the argument list to be passed to the element routine
  int narg = argl.size();
  arg_data_list arg_data_v(narg);

  // pref:= Local values (reference state for finite difference jacobian).
  double *pref;

  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
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

    if (argl[j].options & (UPLOAD_MATRIX | IS_FDJ_MATRIX)) 
      ARGVJ.A = (Mat *) (argl[j].arg);

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

    if (argl[j].options & VECTOR_ASSOC ) 
      arg_data_v[j].vector_assoc = argl[j].arg;
      
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
    
    int chunk_size  = ELEM_CHUNK_SIZE;
    ierr = get_int(elemset->thash,"chunk_size",&chunk_size,1);
    CHKERRQ(ierr);
    chunk_size = mini(2,chunk_size,nel_here);

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

    while (1) { // Loop over chunks. 
      chunk++;

      iele_here = -1;
      for (iele=el_start; iele<nelem; iele++) {
	if (!compute_this_elem(iele,elemset,myrank,iter_mode)) continue;
	iele_here++;
	if (iele_here==chunk_size-1) break;
      }
      el_last = iele;
      if (el_last >=nelem) el_last = nelem-1;
      if (el_last==nelem-1) last_chunk=1;

      for (j=0; j<narg; j++) {
	if (argl[j].options & DOWNLOAD_VECTOR) 
	  elemset->download_vector(nel,ndof,dofmap,ARGVJ,
				   myrank,el_start,el_last,iter_mode,
				   time_data);
      }
      
      elemset->assemble(arg_data_v,nodedata,dofmap,
  			jobinfo,myrank,el_start,el_last,iter_mode,
			time_data);

      // Upload return values
      for (j=0; j<narg; j++) 
	if (argl[j].options & UPLOAD_RETVAL) 
	  elemset->upload_vector(nel,ndof,dofmap,argl[j].options,ARGVJ,
				 myrank,el_start,el_last,iter_mode);

      // compute columns of jacobian matrices by perturbing each
      // local degree of freedom
      if (any_fdj) {

	// copy to reference state
	memcpy (pref,(void *)arg_data_v[j_pert].locst,
		sizeof(double)*chunk_size*ndoft);
	for (j=0; j<narg; j++) 
	  if (argl[j].options & IS_FDJ) 
	    memcpy ((void *)ARGVJ.refres,
		    (void *)ARGVJ.retval,
		    sizeof(double)*chunk_size*ndoft);

	// epsilon:= the increment in the variables in order to
	// compute the finite difference approximation to the
	// Jacobian. Should be order epsilon=sqrt(precision)*(typical
	// magnitude of the variable). Normally, precision=1e-15 so
	// that epsilon=1e-7*(typical magnitude of the
	// variable)
	double epsilon=EPSILON_FDJ;
	for (kloc=0; kloc<nel; kloc++) {
	  for (kdof=0; kdof<ndof; kdof++) {

	    // copy reference state on perturbed  state
	    memcpy ((void *)arg_data_v[j_pert].locst,pref,
		    sizeof(double)*chunk_size*ndoft);
	  
#define PSTAT(iele,j,k) VEC3(arg_data_v[j_pert].locst,iele,j,nel,k,ndof)
	    // perturb state 
	    iele_here=-1;
	    for (iele=el_start; iele <= el_last; iele++) {
	      if (!compute_this_elem(iele,elemset,myrank,iter_mode)) continue;
	      iele_here++;
	      PSTAT(iele_here,kloc,kdof) += epsilon;
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
		      -(RETVALT(iele_here,kdoft)-REFREST(iele_here,kdoft))/epsilon;
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
      el_start = el_last+1;
      if (last_chunk) break;
    } // end loop over chunks

    // report consumed time
    GETOPTDEF(int,report_consumed_time,0);
    // int report_consumed_time=1;
    // ierr = get_int(mesh->global_options,"report_consumed_time",
    // &report_consumed_time_size,1);
    if (report_consumed_time) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Performance report for elemset \"%s\" task \"%s\"\n"
		  "[proc] - total[sec] - rate[sec/element]\n",
		  elemset->type,jobinfo);
      double elapsed;
      elapsed=chrono.elapsed();
      double rate=elapsed/elemset->nelem_here;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,
			      "[proc %d]   %g   %g\n",myrank,elapsed,rate);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }
  }

  if (any_fdj) delete[] pref;
  for (j=0; j<narg; j++) {
    if (argl[j].options & DOWNLOAD_VECTOR) {
      delete[] ARGVJ.locst;
      ierr = VecRestoreArray(*(ARGVJ.ghost_vec),
			     &(ARGVJ.ghost_vals)); CHKERRQ(ierr); 
      ierr = VecDestroy(*(ARGVJ.ghost_vec));
      ierr = VecRestoreArray(*(ARGVJ.x),
			     &(ARGVJ.sstate)); CHKERRQ(ierr); 
      delete ARGVJ.ghost_vec;
    }

    if (argl[j].options & DELETE_RETVAL) 
      delete[] ARGVJ.retval;

    if (argl[j].options & IS_FDJ) {
      delete[] ARGVJ.retval;
      delete[] ARGVJ.refres;
    }

    if (argl[j].options & ASSEMBLY_MATRIX) {
      ierr = MatAssemblyBegin(*(ARGVJ.A),
			      MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*(ARGVJ.A),
			    MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }

    if (argl[j].options & ASSEMBLY_VECTOR) {
      ierr = VecAssemblyBegin(*(ARGVJ.x)); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(*(ARGVJ.x)); CHKERRQ(ierr);
    }

    if (argl[j].options & UPLOAD_PROFILE) {
      //      debug_compute_prof=1;
      ierr = compute_prof(ARGVJ.da,dofmap,
			  myrank,(Mat *)(argl[j].arg));
    }

    if (argl[j].options & VECTOR_ASSOC ) {
      // Perform MIN, MAX and ADD gather operations. 
      vector_assoc_gather(arg_data_v[j].vector_assoc,
 			  argl[j].options,myrank,dofmap->size);
    }
  }
  return 0;
}

