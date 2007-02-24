//__INSERT_LICENSE__
//$Id: state.cpp,v 1.9.18.1 2007/02/24 02:56:42 mstorti Exp $

#include "sttfilter.h"
 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State::axpy(double,const State &)"
State & State::axpy(double alpha,const State &v) {
  int ierr = VecAXPY(&alpha,*(v.vec),*vec); 
  assert(ierr==0); ierr=0;
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State::State(const State &v)"
State::State(const State &v) 
  : vec(NULL), 
    time(v.time) {
  vec = new Vec;
  int ierr = VecDuplicate(*v.vec,vec); ierr=0;
  assert(ierr==0); ierr=0;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State::~State()" 
State::~State() {
  if (vec) {
    int ierr = VecDestroy(*vec);
    assert(ierr==0); ierr=0;
    delete vec;
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "const State & State::print_some(const char *,Dofmap *,set<int> &) const"
const State & State::print_some(const char *filename,Dofmap *dofmap,
	       set<int> & node_list) const {
  ::print_some(filename,v(),dofmap,node_list,&t());
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State & State::scale(double alpha)"
State & State::scale(double alpha) {
  double alpha_ = alpha;
  int ierr = VecScale(&alpha_,*vec);
  assert(ierr==0); ierr=0;
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "State & State::set_cnst(double a)"
State & State::set_cnst(double a) {
  double a_ = a;
  int ierr = VecSet(&a_,*vec);
  assert(ierr==0); ierr=0;
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "const State & State::print() const"
const State & State::print() const {

  printf("time: %f\n",double(t()));
  int ierr = VecView(v(),PETSC_VIEWER_STDOUT_SELF);
  assert(ierr==0); ierr=0;
  return *this;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "const State & State::print(int n) const"
const State & State::print(int n) const {

  printf("time: %f\n",double(t()));
  double *a;
  int lsize;
  int ierr = VecGetLocalSize(v(),&lsize);
  int nn = (n<=0 ||  n>lsize ? lsize : n);
  ierr = VecGetArray(v(),&a);
  assert(ierr==0); ierr=0;
  for (int j=0; j<nn; j++) {
    printf("%d %f\n",j,a[j]);
  }
  return *this;
}
