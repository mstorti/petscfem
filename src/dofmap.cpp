//__INSERT_LICENSE__
//$Id: dofmap.cpp,v 1.22.10.1 2007/02/19 20:23:56 mstorti Exp $

#include <cassert>
#include <algorithm>

#include <petsc.h> 

#include "libretto.h"
#include <libretto/darray.h>

using namespace std;

#include <newmat.h>

#include "fem.h"
#include "getprop.h"
#include "fstack.h"
#include "dofmap.h"
#include "utils.h"
#include "util2.h"
#include "fastlib.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void fixation_entry::print(void) const"
void fixation_entry::print(void) const {
  PetscPrintf(PETSCFEM_COMM_WORLD,"Temporal Amplitude function:\n");
  if (amp!=NULL) {
    amp->print();
  } else {
    PetscPrintf(PETSCFEM_COMM_WORLD,"No temporal part!\n");
    PetscPrintf(PETSCFEM_COMM_WORLD,"value: %f\n",val);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double fixation_entry::value(void *) const"
double Dofmap::value(const fixation_entry &fe, const TimeData *time_data) const {
  // double fixation_entry::value(const TimeData *time_data) const {
  if (fe.amp==NULL) {
    return fe.val;
  } else {
    assert(time_data);
    double v;
    if (!fe.amp->needs_dof_field_q()) {
      v = fe.val*fe.amp->eval(time_data);
    } else {
      int node,field;
      nodf(fe.edof,node,field);
      v = fe.val*fe.amp->eval(time_data,node,field);
    }
    return v;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fixa_entry_cmp" 
int fixa_entry_cmp (const void *left,const void *right, void *args) {
  fixa_entry *l = (fixa_entry *) left;
  fixa_entry *r = (fixa_entry *) right;
  if (l->node > r->node) return +1;
  if (l->node < r->node) return -1;
  if (l->kdof > r->kdof) return +1;
  if (l->kdof < r->kdof) return -1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_row_free" 
void Dofmap::get_row_free(int const & node,int const & kdof,row_t
			  &row) const {
  int edoff = edof(node,kdof);
  id->get_row(edoff,row);
  row_t::iterator j;
  for (j=row.begin(); j!=row.end(); j++) {
    if (j->first > neq) row.erase(j);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::row_set" 
void Dofmap::row_set(const int & node,const int & kdof,const row_t & row) {
  int edoff = edof(node,kdof);
  id->row_set(edoff,row);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_row" 
void Dofmap::get_row(int const & node,
		     int const & kdof,row_t &row) const {
  int edoff = edof(node,kdof);
  id->get_row(edoff,row);
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_row" 
void Dofmap::get_row(int const & node,
		     int const & kdof,IdMapRow &row) const {
  int edoff = edof(node,kdof);
  id->get_row(edoff,row);
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int Dofmap::edof(const int node,const int field) const"
int Dofmap::edof(const int node,const int field) const {
#if !defined(NDEBUG)
  PETSCFEM_ASSERT(node>=1 && node<=nnod,"Node out of range: %d\n",node);
  PETSCFEM_ASSERT(field>=1 && field<=ndof,"Field out of range: %d\n",field);
#endif
  return (node-1)*ndof+field;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_dofval"
double  Dofmap::get_dofval(const int & jeq,double const *sstate,
			   double const *ghost_vals,const TimeData *time_data) const {
  vector<int>::iterator it;
  if (dof1 <= jeq && jeq <= dof2) {
    return sstate[jeq-dof1];
  } else if (jeq<=neq) {
    it = lower_bound(ghost_dofs->begin(),ghost_dofs->end(),jeq-1);
    assert(*it == jeq-1);
    return *(ghost_vals+(it-ghost_dofs->begin()));
  } else {
    // printf("jeq %d -> \n",jeq);
    // fixed[jeq-neq-1].print();
    return value(fixed[jeq-neq-1],time_data);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_dofval"
double  Dofmap::get_dofval(const int & jeq,
			   const double *sstate,const TimeData *time_data) const{
  if (jeq<=neq) {
    return sstate[jeq-1];
  } else {
    // printf("jeq %d -> %f\n",jeq,fixed[jeq-neq-1].first);
    return value(fixed[jeq-neq-1],time_data);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_nodal_value" 
int Dofmap::get_nodal_value(int const & node,int const & kdof,double
			    const * sstate,double const *ghost_vals,
			    const TimeData *time_data,double & value) const {

  double w=0.;
#if 0
  // Old slow version
  IdMapRow row;
  IdMapEntry *entry;
  // row_t row;
  get_row(node,kdof,row);
  for (int k=0; k<row.size(); k++) {
    entry = &row[k];
    //for (k=row.begin(); k!=row.end(); k++) {
    int keq= entry->j;
    double coef= entry->coef;
    if (row.size()>1) {
      // printf("eq/coef: %d %g\n",keq,coef);
    }
    w += coef*get_dofval(keq,sstate,ghost_vals,time_data);
  }
  value = w;
#else
  // New fast version
  int n;
  const int *dofp, *dofp_end;
  const double *coefp;
  get_row(node,kdof,n,&dofp,&coefp);
  if (n==1) {
    value = (*coefp) * get_dofval((*dofp),sstate,ghost_vals,time_data);
  } else {
    dofp_end = dofp+n;
    while (dofp<dofp_end) {
      // printf("eq/coef: %d %g\n",*dofp,*coefp);
      w += (*coefp++)
	*get_dofval((*dofp++),sstate,ghost_vals,time_data);
    }
    value = w;
  }
#endif
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_nodal_value" 
int Dofmap::get_nodal_value(int const & node,int const & kdof,double
			    const * sstate,const TimeData *time_data,double & value
			    ) const {

  value=0.;
  row_t row;
  get_row(node,kdof,row);
  row_t::iterator k;
  for (k=row.begin(); k!=row.end(); k++) {
    int keq= k->first;
    double coef=k->second;
    value += coef*get_dofval(keq,sstate,time_data);
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "fixa_entry::print" 
void fixa_entry::print(void) {
  printf("node: %d, kdof %d, val: %f\n",node,kdof,val);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "q_entry_cmp" 
int q_entry_cmp (const void *left,const void *right, void *args) {
  q_entry *l = (q_entry *) left;
  q_entry *r = (q_entry *) right;
  if (l->node > r->node) return +1;
  if (l->node < r->node) return -1;
  if (l->kdof > r->kdof) return +1;
  if (l->kdof < r->kdof) return -1;
  if (l->keq > r->keq) return +1;
  if (l->keq < r->keq) return -1;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Dofmap::col_is_null(int j) {
  static row_t *col;
  if (col==NULL) col= new row_t;
  id->get_col(j,*col);
  return col->size()==0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Constraint::add_entry(int node,int field,double coef)"
void Constraint::add_entry(int node,int field,double coef) {
  constraint_entry ce;
  ce.node = node;
  ce.field = field;
  ce.coef = coef;
  //  ce_list.push_back(ce);
  push_back(ce);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void Constraint::empty(void)"
void Constraint::empty(void) {
  VOID_IT((*this));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void dof_range(int ,int &,int &)"
void Dofmap::dof_range(int myrank,int &dof1,int &dof2) {
  dof1=startproc[myrank];
  dof2=dof1+neqproc[myrank]-1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::create_MPI_vector(Vec &)"
int Dofmap::create_MPI_vector(Vec &v) {

  int myrank;
  MPI_Comm_rank(PETSCFEM_COMM_WORLD,&myrank);
  int ierr = VecCreateMPI(PETSCFEM_COMM_WORLD,neqproc[myrank],neq,&v);
  CHKERRQ(ierr);
  return ierr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::create_MPI_ghost_vector(Vec &)"
int Dofmap::create_MPI_ghost_vector(Vec &v) {
  
  int myrank, ierr;
  MPI_Comm_rank(PETSCFEM_COMM_WORLD,&myrank);
  ierr = VecCreateSeq(PETSC_COMM_SELF,
		      ghost_dofs->size(),&v); CHKERRQ(ierr);
  return ierr;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Dofmap::processor(int j) const {
  int size,rank;
  MPI_Comm_size(PETSCFEM_COMM_WORLD, &size);

  for (rank=0; rank<size; rank++)
    if (j>=startproc[rank] && j<startproc[rank+1]) break;
  return rank;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Dofmap::Dofmap() : 
  idmap2(NULL), special_ptr(NULL), sp_eq(NULL), coefs(NULL),
  comm(PETSCFEM_COMM_WORLD),
  nnod(0), neq(0), dof1(0), dof2(0), neqf(0), neqtot(0), ndof(0),
  ident(NULL), ghost_dofs(NULL), id(NULL),  fixa(NULL), q(NULL),
  startproc(NULL), neqproc(NULL), size(0), tpwgts(NULL), npart(NULL),
  ghost_scatter(NULL), scatter_print(NULL)
{ }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Dofmap::~Dofmap() {
  DELETE_SCLR(id);
  DELETE_VCTR(tpwgts);
}
