//__INSERT_LICENSE__
//$Id: pfmat.cpp,v 1.13 2003/08/28 22:43:16 mstorti Exp $

#include <petscmat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>
// #pragma implementation "PFMat"
#include <src/pfmat.h>
#include <src/pfptscmat.h>
#include <src/iisdmat.h>
#include <src/petscmat.h>
#include <src/spdirect.h>

#define PF_ACTION_DECL(action) void action() 

#define PF_ACTION_DEF(action)			\
void pfmatFSMContext::action() {		\
  matrix_p->ierr = matrix_p->action();		\
  if (matrix_p->ierr)				\
    printf("pfmatFSMContext::action ierr=%d\n",	\
	   matrix_p->ierr);			\
}

#define PF_ACTION_LIST				\
  PF_ACTION(clean_prof_a);			\
  PF_ACTION(clean_mat_a);			\
  PF_ACTION(clean_factor_a);			\
  PF_ACTION(factor_and_solve_A);		\
  PF_ACTION(solve_only_A);

#define PF_ACTION(name) PF_ACTION_DECL(name) 
class pfmatFSMContext {
public:
  PFMat * matrix_p;
  // pfmatFSMContext(PFMat *p) : matrix_p(p) {};
  PF_ACTION_LIST;
  void FSMError(const char *e,const char *s) { 
    printf("PFMat: Not valid event \"%s\" in state \"%s\"\n",e,s);
  }
};
#undef PF_ACTION

#define PF_ACTION(name) PF_ACTION_DEF(name) 
PF_ACTION_LIST;
#undef PF_ACTION

class PFMat;

#include "pfmatFSM.h"

const char *from_g = "START", *event_g = "START", *to_g = "START";
int transition_counter = 0;

#define DEBUG_FSM_TRANSITION
#ifdef DEBUG_FSM_TRANSITION
#define PRINT_FSM_TRANSITION_INFO(from,event,to)	\
  int mode = s.matrix_p->print_fsm_transition_info_f();	\
  report_transition_event(from,event,to,mode)
#endif

void report_transition_event(const char *from, const char *event,
			     const char *to,int mode) {
  if (mode==1) {
      printf("[%d] from: \"%s\", event \"%s\", to \"%s\"\n",
	     MY_RANK,from_g,event_g,to_g);
  } else if (mode==2) {
    if (!(!strcmp(from,from_g) && !strcmp(event,event_g) && !strcmp(to,to_g))) {
      if (transition_counter) 
	printf("[%d] from: \"%s\", event \"%s\", to \"%s\"",
	       MY_RANK,from_g,event_g,to_g);
      if (transition_counter>1) printf(" (%d times)",transition_counter);
      printf("\n");
      transition_counter = 0;
      from_g = from;
      to_g = to;
      event_g = event;
    }
    transition_counter++;
  }
}

#include "pfmatFSM.cpp"

PFMat::PFMat() : ierr(0), 
  print_fsm_transition_info(0), fsm(new pfmatFSM) { 
  fsm->matrix_p = this; }

PFMat::~PFMat() { delete fsm; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_profile"
int PFMat::set_profile(int row,int col) {
  fsm->set_profile();
  CHKERRQ(ierr); 
  
  ierr = set_profile_a(row,col);
  CHKERRQ(ierr); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::create"
int PFMat::create() {
  fsm->create();
  CHKERRQ(ierr); 
  
  ierr = create_a();
  CHKERRQ(ierr); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_value"
int PFMat::set_value(int row,int col,PetscScalar value,
		     InsertMode mode) {
  fsm->set_value();
  CHKERRQ(ierr); 
  
  ierr = set_value_a(row,col,value,mode); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_values"
int PFMat::set_values(int nrows,int *idxr,int ncols,int *idxc,
		      PetscScalar *values, InsertMode mode) { 
  fsm->set_value();
  CHKERRQ(ierr); 
  
  ierr = set_values_a(nrows,idxr,ncols,idxc,values,mode); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::set_values_a"
int PFMat::set_values_a(int nrows,int *idxr,int ncols,int *idxc,
			PetscScalar *values, InsertMode mode) { 
  int row, ierr=0;
  for (int j=0; j<nrows; j++) {
    row = idxr[j];
    for (int k=0; k<ncols; k++) {
      ierr = set_value_a(row,idxc[k],values[j*ncols+k]);
      if (ierr) break;
    }
    if (ierr) break;
  }
  return ierr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly_begin"
int PFMat::assembly_begin(MatAssemblyType type) {
  fsm->assembly_begin();
  CHKERRQ(ierr); 
  
  ierr = assembly_begin_a(type); CHKERRQ(ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly_end"
int PFMat::assembly_end(MatAssemblyType type) {
  fsm->assembly_end();
  CHKERRQ(ierr); 
  
  ierr = assembly_end_a(type); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::assembly"
int PFMat::assembly(MatAssemblyType type) {
  int ierr;
  ierr = assembly_begin(type); CHKERRQ(ierr); 
  ierr = assembly_end(type); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::solve"
int PFMat::solve(Vec &res,Vec &dx) {
  res_p = &res;
  dx_p = &dx;
  fsm->solve(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::factor_and_solve"
int PFMat::factor_and_solve(Vec &res,Vec &dx) {
  res_p = &res;
  dx_p = &dx;
  fsm->factor_and_solve(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::solve_only"
int PFMat::solve_only(Vec &res,Vec &dx) {
  res_p = &res;
  dx_p = &dx;
  fsm->solve_only(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::duplicate"
int PFMat::duplicate(MatDuplicateOption op,const PFMat &A) {
  fsm->clear(); 
  duplicate_a(op,A);
  fsm->set_profile();
  fsm->set_value(); 
  return 0;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_factor"
int PFMat::clean_factor() { 
  fsm->clean_factor(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_mat"
int PFMat::clean_mat() { 
  fsm->clean_mat(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clean_prof"
int PFMat::clean_prof() { 
  fsm->clean_prof(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "PFMat::clear"
int PFMat::clear() { 
  fsm->clear(); CHKERRQ(ierr); 
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
PFMat * PFMat::dispatch(int N,DofPartitioner &part,const char *s) {
  PFMat *A;
  IISDMat *AA;
  // IISD solver with PETSc or SuperLU local solver
  if (!strcmp(s,"iisd_superlu")) {
#ifdef USE_SUPERLU
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    AA->local_solver = IISDMat::SuperLU;
    return AA;
#else
    PETSCFEM_ERROR0("Not compiled with SuperLU library!!\n");
#endif
  } else if (!strcmp(s,"iisd_petsc")) {
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    AA->local_solver = IISDMat::PETSc;
    return AA;
  } else if (!strcmp(s,"iisd")) {
    // local solver is chosen by default
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    return AA;
  } else if (!strcmp(s,"petsc")) {
    // PETSc (iterative) solver 
    A = new PETScMat(N,N,part,PETSC_COMM_WORLD);
    return A;
  } else if (!strcmp(s,"direct_superlu")) {
#ifdef USE_SUPERLU
    A = new SparseDirect(N,"SuperLU");
    return A;
#else
    PETSCFEM_ERROR0("Not compiled with SuperLU library!!\n");
#endif
  } else if (!strcmp(s,"direct") || 
	     !strcmp(s,"direct_petsc")) {
    A = new SparseDirect(N,"PETSc");
    return A;
  } else {
    PETSCFEM_ERROR("PFMat type not known: %s\n",s);
  }
}

