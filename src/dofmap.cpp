/*__INSERT_LICENSE__*/
//$Id: dofmap.cpp,v 1.4 2001/05/30 03:58:50 mstorti Exp $

#include <cassert>
#include <algorithm>

#include <petsc.h> 

#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <libretto/darray.h>

#include <newmat.h>

#include "fem.h"
#include "getprop.h"
#include "fstack.h"
#include "dofmap.h"
#include "utils.h"
#include "util2.h"
#include "fastlib.h"

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void fixation_entry::print(void) const"
void fixation_entry::print(void) const {
  PetscPrintf(PETSC_COMM_WORLD,"Temporal Amplitude function:\n");
  if (amp!=NULL) {
    amp->print();
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"No temporal part!\n");
    PetscPrintf(PETSC_COMM_WORLD,"value: %f\n",val);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double fixation_entry::value(void *) const"
double fixation_entry::value(const TimeData *time_data) const {
  if (amp==NULL) {
    return val;
  } else {
    return val*amp->eval(time_data);
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
};

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

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_column"
int Dofmap::get_column(int const & node,int const & kdof,Darray *da) {
  sp_entry spe;
  da_resize(da,0);
  int locdof = IDENT(node-1,kdof-1);
  spe = sp_entry(node,kdof,locdof,1.);
  if (locdof>0) da_append(da,&spe);
}
#endif

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
    return fixed[jeq-neq-1].value(time_data);
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
    return fixed[jeq-neq-1].value(time_data);
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_nodal_value" 
int Dofmap::get_nodal_value(int const & node,int const & kdof,double
			    const * sstate,double const *ghost_vals,
			    const TimeData *time_data,double & value) const {

  value=0.;
  row_t row;
  get_row(node,kdof,row);
  row_t::iterator k;
  for (k=row.begin(); k!=row.end(); k++) {
    int keq= k->first;
    double coef=k->second;
    value += coef*get_dofval(keq,sstate,ghost_vals,time_data);
  }
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::get_nodal_value" 
int Dofmap::get_nodal_value(int const & node,int const & kdof,double
			    const * sstate,double const *ghost_vals,
			    const TimeData *time_data,double & value) const {

  value=0.;
  IdMapRow row;
  IdMapEntry *entry;
  // row_t row;
  get_row(node,kdof,row);
  for (int k=0; k<row.size(); k++) {
    entry = &row[k];
    //for (k=row.begin(); k!=row.end(); k++) {
    int keq= entry->j;
    double coef= entry->coef;
    value += coef*get_dofval(keq,sstate,ghost_vals,time_data);
  }
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

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::add_q_entry"
void Dofmap::add_q_entry(const int node,const int kdof,const int keq,
			 const double coef) {
  q_entry qe(node,kdof,keq,coef);
  da_append(q,&qe);
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "dofmap::synchronize"
void Dofmap::synchronize() {
  // Resync fixa and q darray's
  da_sort (fixa,fixa_entry_cmp,NULL);
  da_sort (q,q_entry_cmp,NULL);

#if 0  
  int nq = da_length(q);
  for (int kk=0; kk<nq; kk++) {
    q_entry *qe = (q_entry *) da_ref(q,kk);
    printf("kk %d, node %d, kdof %d, keq %d, coef %f\n",kk,
	   qe->node,qe->kdof,qe->keq,qe->coef);
  }
#endif

  synchro=1;
}
#endif

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
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "Dofmap::print"
void Dofmap::print (char * s=NULL) {
  darray *row = da_create(sizeof(sp_entry));
  sp_entry *spe;
  PetscPrintf(PETSC_COMM_WORLD,
	      (s == NULL ? "Dofmap dump: \n" : s));
  for (int node=1; node<=nnod; node++) {
    for (int kdof=1; kdof<=ndof; kdof++) {
      get_row(node,kdof,row);
      int len=da_length(row);
      printf("row for node %d, dof %d has %d entries\n",node,kdof,len);
      for (int kk=0; kk<len; kk++) {
	spe = (sp_entry *) da_ref(row,kk);
	printf("dof %d -> coef %f\n",spe->gdof,spe->coef);
      }
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,
	      "--------- end of dofmap dump. ---------\n");
  da_destroy(row);
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Dofmap::col_is_null(int j) {
  static row_t *col;
  if (col==NULL) col= new row_t;
  id->get_col(j,*col);
  return col->size()==0;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Dofmap::set_fixation(int node,int kdof,double val) {
  int edoff = edof(ndof,kdof);
  id->set_elem(edof,edof,0);
  fixed[edoff] = pair(val,1);
}
#endif

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
#define __FUNC__ "Dofmap::set_constraint(const Constraint &constraint)"
#if 1
void Dofmap::set_constraint(const Constraint &constraint) {

  int length=constraint.size();
  assert(length>1);
  const constraint_entry *it,*it0;
  row_t row,row0;
  VOID_IT(row0);
  // Look for the maximum (abs. val.) coefficient
  int jmax=0; double cmax=fabs(constraint[0].coef);
  for (int k=1; k<length; k++) {
    double cc=fabs(constraint[k].coef);
    if (cc>cmax) {
      jmax=k;
      cmax=cc;
    }
  }
  
  it0 = &(constraint[jmax]);
  double coef0 = it0->coef;
  for (int k=0; k<length; k++) {
    if (k==jmax) continue;
    it = &(constraint[k]);
    get_row(it->node,it->field,row); 
    axpy(row0,-it->coef/coef0,row);
  }
  row_set(it0->node,it0->field,row0);
}
#else // OLD VERSION
void Dofmap::set_constraint(const Constraint &constraint) {

  int length=constraint.size();
  assert(length>1);
  const constraint_entry *it,*it0;
  row_t row,row0;
  VOID_IT(row0);
  it0 = &(constraint[0]);
  double coef0 = it0->coef;
  for (int k=1; k<length; k++) {
    it = &(constraint[k]);
    get_row(it->node,it->field,row); 
    axpy(row0,-it->coef/coef0,row);
  }
  row_set(it0->node,it0->field,row0);
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void dof_range(int ,int &,int &)"
void Dofmap::dof_range(int myrank,int &dof1,int &dof2) {
  dof1=startproc[myrank];
  dof2=dof1+neqproc[myrank]-1;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void create_ghosted_vector(Vec *v)"
int Dofmap::create_ghosted_vector (Vec *v) {

  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  int ierr = VecCreateGhost(PETSC_COMM_WORLD,neqproc[myrank],
		      neq,ghost_dofs->size(),ghost_dofs->begin(),
		      v); CHKERRQ(ierr);

}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::create_MPI_vector(Vec &)"
int Dofmap::create_MPI_vector(Vec &v) {

  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  int ierr = VecCreateMPI(PETSC_COMM_WORLD,neqproc[myrank],neq,&v);
  CHKERRQ(ierr);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::create_MPI_ghost_vector(Vec &)"
int Dofmap::create_MPI_ghost_vector(Vec &v) {
  
  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
//    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
//  	      "[%d] Dofmap::create_MPI_ghost_vector ghost_dofs->size(): %d\n",
//  	      myrank, ghost_dofs->size());
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  int ierr = VecCreateSeq(PETSC_COMM_SELF,
			  ghost_dofs->size(),&v); CHKERRQ(ierr);
}

