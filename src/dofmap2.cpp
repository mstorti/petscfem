//__INSERT_LICENSE__
//$Id: dofmap2.cpp,v 1.3 2002/12/22 19:45:42 mstorti Exp $

#include <cassert>
#include <deque>

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
#include "fastmat2.h"

using namespace std;

// Computes the inverse of `edof': given an edof value gives
// the node/field combination
void edofi(int edof, int ndof, int &node,int &field) {
  field = modulo(edof-1,ndof,&node)+1;
  node = node+1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::set_constraint(const Constraint &constraint)"
void Dofmap::set_constraint(const Constraint &constraint) {

  double tol = 1e-10;
  int length=constraint.size();
  assert(length>1);
  const constraint_entry *it,*it0;
  row_t row,row0;
  VOID_IT(row0);
  // Look for the maximum (abs. val.) coefficient
  // among all rows that have not been constrained already
  // (I'm not sure that this is the RIGHT WAY to do this). 
  int set_flag=0, jmax=0; double cmax=0.;
  for (int k=0; k<length; k++) {
    row_t row;
    it = &(constraint[k]);
    get_row(it->node,it->field,row); 
    if (row.size()!=1) continue;
    row_t::iterator q = row.begin();
    int edoff = edof(it->node,it->field);
    if (q->first!=edoff) continue;
    double cc=fabs(constraint[k].coef);
    if (!set_flag || cc>cmax) {
      set_flag = 1;
      jmax=k;
      cmax=cc;
    }
  }
  PETSCFEM_ASSERT0(set_flag,"Couldn't find an unconstrained dof for this constraint.");  

  // First insert the row `blindly'
  it0 = &(constraint[jmax]);
  double coef0 = it0->coef;
  for (int k=0; k<length; k++) {
    if (k==jmax) continue;
    it = &(constraint[k]);
    row0.insert(pair<int,double>(edof(it->node,it->field),-it->coef/coef0));
  }
  row_set(it0->node,it0->field,row0);

  // Now check for recursive contraints, i.e.  Dofs to be eliminated
  // are those connected cyclically (in the sense of a graph) to the
  // current one, i.e. edof `j' is connected to edof `i' if node `i'
  // depends on `j' and also `j' depends on `i'.  In that case we
  // build a system `Q edofs = B nedofs' where `edofs' are the dofs to
  // be eliminated and nedofs thos that aren't. Then we invert Q and
  // replace all rows corresponding to edofs by their respective rows
  // in Q^{-1} B.

  // Build set of dofs to be eliminated
  // Traverses the graph `breadth-first'
  set<int> to_elim;
  deque<int> elim_front;
  int edof0 = edof(it0->node,it0->field);
  to_elim.insert(edof0);
  elim_front.push_back(edof0);
  while (elim_front.size()>0) {
    // Take first edof in the queue
    int edof_elim = elim_front.front();
    elim_front.pop_front();
    row_t row;
    int node,field;
    edofi(edof_elim,ndof,node,field);
    get_row(node,field,row); 
    row_t::iterator l,le;
    le=row.end();
    for (l=row.begin(); l!=le; l++) {
      int edof_c = l->first;
      if (to_elim.find(edof_c)!=to_elim.end()) continue;
      row_t row_c;
      int node_c,field_c;
      edofi(edof_c,ndof,node_c,field_c);
      get_row(node_c,field_c,row_c); 
      row_t::iterator q,qe;
      qe = row_c.end();
      for (q=row_c.begin(); q!=qe; q++) {
	if (to_elim.find(q->first)==to_elim.end()) continue;
	to_elim.insert(edof_c);
	elim_front.push_back(edof_c);
      }
    }
  }

#if 0
  printf("to_elim:");
  for (set<int>::iterator q=to_elim.begin(); 
       q!=to_elim.end(); q++) printf(" %d",*q);
  printf("\n");
#endif

#if 0
  for (int k=0; k<length; k++) {
    if (k==jmax) continue;
    it = &(constraint[k]);
    get_row(it->node,it->field,row); 
    int edof1 = edof(it->node,it->field);
    row_t::iterator l,le;
    le=row.end();
    for (l=row.begin(); l!=le; l++) {
      int edoff = l->first;
      if (to_elim.find(edoff)!=to_elim.end()) to_elim.insert(edof1);
    }
  }
#endif
  
  // Number of dofs to be eliminated
  int nelim = to_elim.size();
  FastMat2 Q(2,nelim,nelim),iQ(2,nelim,nelim);
  vector<row_t *> B(nelim);

  // Costruct map from edofs to index in the set to be eliminated
  map<int,int> edof2elim;
  vector<int> elim2edof(nelim);
  int j=0;
  for (set<int>::iterator q=to_elim.begin(); 
       q!=to_elim.end(); q++) {
    edof2elim[*q] = ++j;
    elim2edof[j-1] = *q;
  }

#if 0
  // First line are the coefficients of the constraint
  B[0] = new row_t;
  for (int k=0; k<length; k++) {
    it = &(constraint[k]);
    int edof1 = edof(it->node,it->field);
    if (to_elim.find(edof1)!=to_elim.end()) {
      Q.setel(it->coef,1,edof2elim[edof1]);
    } else {
      B[0]->insert(pair<int,double>(edof1,it->coef));
    }
  }
#endif

  Q.set(0.);
  for (int r=1; r<=nelim; r++) {
    B[r-1] = new row_t;
    row_t row_o;		// the original row
    int edoff,node, field;
    edoff = elim2edof[r-1];
    edofi(edoff,ndof,node,field);
    get_row(node,field,row_o); 
    // Iterates on the original row. Columns in `to_elim'
    // are passed to the Q matrix, and the others to their respective
    // rows in `B'
    for (row_t::iterator q=row_o.begin(); 
	 q!=row_o.end(); q++) {
      int edof_c = q->first;
      if (to_elim.find(edof_c)!=to_elim.end()) {
	Q.setel(q->second,r,edof2elim[edof_c]);
      } else {
	B[r-1]->insert(*q);
      }
    }
    Q.addel(-1.,r,r);		// The row shouldn't be a `regular' row
    assert(fabs(Q.get(r,r))>tol);
  }

#if 0
  Q.print("Q:");
  for (int k=0; k<nelim; k++) {
    printf("%d:",k+1);
    print(*B[k]);
    printf("\n");
  }
#endif

  iQ.inv(Q);

  for (int k=0; k<nelim; k++) {
    row_t row,roww,row_q;
    for (int l=0; l<nelim; l++) axpy(row,-iQ.get(k+1,l+1),*B[l]);
    row_t::iterator q,qe=row.end();
    for (q=row.begin(); q!=qe; q++) {
      int node,field;
      edofi(q->first,ndof,node,field);
      get_row(node,field,row_q); 
      axpy(roww,q->second,row_q);
    }
    int edoff = elim2edof[k];
    int node,field;
    edofi(edoff,ndof,node,field);
    row_set(node,field,roww);
  }

  Q.clear();
  iQ.clear();
  // Free memory in B
  for (int k=0; k<nelim; k++) delete B[k];
}
