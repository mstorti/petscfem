/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

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

#if 0
    Darray *gh = elemset->ghost_elems;
    if (myrank==0) {
      printf("En el proc. 0, lista de ghost_elems\n");
      for (int kk=0; kk<da_length(gh); kk++) {
	int elem = *(int *) da_ref(gh,kk);
	printf("%d -> %d\n",kk,elem);
      }
    }
#endif

    int is_ghost_elem = da_bsearch(elemset->ghost_elems,&iele,int_cmp,NULL);
//      if (myrank==0) {
//        printf("elem %d,  is_ghost_elem = %d\n",iele,is_ghost_elem);
//      }
//      if (is_ghost_elem >=0) {
//        PetscPrintf(PETSC_COMM_WORLD,"Element %d is ghost on pocessor %d\n",
//  		  iele,myrank);
//      }
    printf("processing ghost_elem %d: ... %d\n",iele,is_ghost_elem>=0);
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
    load_vec,load_mat,load_mat_col,comp_prof,comp_mat,iele_here;

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
		} else {
		  MatSetValue(*argd.A,kd,kdl,val,ADD_VALUES); 
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
		    node_insert(argd.da,locdof-1,locdofl-1);
		  } else {
		    // printf("(%d,%d) -> %f\n",locdof,locdofl,val);
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
	     Dofmap *dofmap,char *jobinfo,const TimeData *time_data=NULL) {

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
    // int chunk_size = ELEM_CHUNK_SIZE;
    // ierr = get_int(elemset->thash,"chunk_size",&chunk_size,1);

    int local_chunk_size;
    // scaled chunk_size in order to balance processors 
    local_chunk_size = (int)(chunk_size*dofmap->tpwgts[myrank]/w_max) +1 ;

    CHKERRQ(ierr);
    chunk_size = mini(2,local_chunk_size,nel_here);

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

      for (j=0; j<narg; j++) {
	if (argl[j].options & ASSEMBLY_MATRIX) {
	  ierr = MatAssemblyBegin(*(ARGVJ.A),
				  MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(*(ARGVJ.A),
				MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
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
    GETOPTDEF(int,report_consumed_time,0);
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
      // wait_from_console("antes de assembly FINAL en FINAL");  
      ierr = MatAssemblyBegin(*(ARGVJ.A),
			      MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*(ARGVJ.A),
			    MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      // wait_from_console("despues de assembly FINAL en FINAL");  
    }

    if (argl[j].options & ASSEMBLY_VECTOR) {
      ierr = VecAssemblyBegin(*(ARGVJ.x)); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(*(ARGVJ.x)); CHKERRQ(ierr);
    }

    if (argl[j].options & UPLOAD_PROFILE) {
      ierr = compute_prof(ARGVJ.da,dofmap,
			  myrank,(Mat *)(argl[j].arg),
			  debug_compute_prof);
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

