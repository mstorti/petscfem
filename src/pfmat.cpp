//__INSERT_LICENSE__
//$Id: pfmat.cpp,v 1.4 2002/09/05 18:23:52 mstorti Exp $

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
int PFMat::set_value(int row,int col,Scalar value,
		     InsertMode mode=ADD_VALUES) {
  fsm->set_value();
  CHKERRQ(ierr); 
  
  ierr = set_value_a(row,col,value,mode); CHKERRQ(ierr); 
  return 0;
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
    AA =  new IISDMat(N,N,part,PETSC_COMM_WORLD);
    AA->local_solver = IISDMat::SuperLU;
    return AA;
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
    A = new SparseDirect(N,"SuperLU");
    return A;
  } else if (!strcmp(s,"direct") || 
	     !strcmp(s,"direct_petsc")) {
    A = new SparseDirect(N,"PETSc");
    return A;
  } else {
    PETSCFEM_ERROR("PFMat type not known: %s\n",s);
  }
}

