//__INSERT_LICENSE__
//$Id: dofmap2.cpp,v 1.6 2002/12/25 14:44:11 mstorti Exp $

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
#define DEBUG_ELIM

  // This implementation takes into account cyclic references. 
  // Cyclic reference means that edof `i' is constrained to
  // edof `j' but edof `j' has been already constrained to edof `i'.
  // So we cannot straightforwardly eliminate edof `i' in terms of `j'
  // but instead we haqve to build a system like Q u_e + B u_ne = 0
  // where `u_e' are the dofs to eliminate (i and j) and `u_ne' the
  // rest. 

  double tol = 1e-10;
  int length=constraint.size();
  assert(length>1);
  const constraint_entry *it,*it0;
  row_t row,row0;
  
  // Find edof to eliminate: that one that is not already
  // eliminated and has maximum coef


  // Look for the maximum (abs. val.) coefficient
  // among all rows that have not been constrained already
  // (I'm not sure that this is the RIGHT WAY to do this). 

  // set_flag:= flags whether a regular row has been found already or
  // not.  
  // jmax:= the index of the regular row with higher absolute value
  // coefficient
  int set_flag=0, jmax=0; double cmax=0.;
  for (int k=0; k<length; k++) {
    it = &(constraint[k]);
    int edoff = edof(it->node,it->field);
    // Look if row is regular (it doesn't point to other edof's)
    if (col_is_null(edoff)) continue;
    // Update maximum
    double cc=fabs(constraint[k].coef);
    if (!set_flag || cc>cmax) {
      set_flag = 1;
      jmax=k;
      cmax=cc;
    }
  }

  // Create row from constraint and eliminate all already eliminated
  // variables, except those that refer to this edof 
  for (int k=0; k<length; k++) {
    it = &(constraint[k]);
    int edoff = edof(it->node,it->field);
    // Look if row is regular (it doesn't point to other edof's)
    row.clear();
#ifdef DEBUG_ELIM
    printf("row item: node %d, dof %d, edof %d, coef %f\n",
	   it->node,it->field,edof(it->node,it->field),
	   it->coef);
#endif    
    if (col_is_null(edoff)) {
      get_row(it->node,it->field,row); 
    } else {
      row[edoff] = 1.;
    }
    axpy(row0,it->coef,row);
  }

#ifdef DEBUG_ELIM
  print(row0,"entered row:");
#endif

  int node, field;
  edofi(edof0,ndof,node,field);
  row0.erase(qmax);
  row.clear();
  axpy(row,-1./coef0,row0);
  row_set(node,field,row);

#ifdef DEBUG_ELIM
  print(row0,"Eliminated unknown:");
#endif

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

  // to_elim:= set of nodes to be eliminated
  set<int> to_elim;
  // elim_front:= while visiting the graph it keeps the front
  // of edof's that are inserted in `to_elim' but doesn't have
  // checked already their neighbors for cyclic reference. 
  deque<int> elim_front;
  to_elim.insert(edof0);
  elim_front.push_back(edof0);
  // This is the typical algorithm for traversing a graph
  // `breadth-first'. We keep a queue with the nodes to
  // be yet visited. The algorithm ends when the queue is empty.
  while (elim_front.size()>0) {
    // Take first edof in the queue
    int edof_elim = elim_front.front();
    elim_front.pop_front();
    // get row for this edof
    row_t row;
    int node,field;
    edofi(edof_elim,ndof,node,field);
    get_row(node,field,row); 
    // iterate over the row
    row_t::iterator l,le;
    le=row.end();
    for (l=row.begin(); l!=le; l++) {
      int edof_c = l->first;
      // edof_c:= is an edof constrained to `edof_elim' if it is
      // already marked to be eliminated then we have nothing to do at
      // this stage.xs
      if (to_elim.find(edof_c)!=to_elim.end()) continue;
      // get row for `edof_c'
      row_t row_c;
      int node_c,field_c;
      edofi(edof_c,ndof,node_c,field_c);
      get_row(node_c,field_c,row_c); 
      // iterate over the row for `edof_c' and check if some edof
      // is constrained to an edof already marked to be eliminated
      row_t::iterator q,qe;
      qe = row_c.end();
      for (q=row_c.begin(); q!=qe; q++) {
	if (to_elim.find(q->first)==to_elim.end()) continue;
	to_elim.insert(edof_c);
	elim_front.push_back(edof_c);
      }
    }
  }

#ifdef DEBUG_ELIM
  // Prints set of edofs to be eliminated for this constraint
  printf("to_elim:");
  for (set<int>::iterator q=to_elim.begin(); 
       q!=to_elim.end(); q++) printf(" %d",*q);
  printf("\n");
#endif

  // Number of dofs to be eliminated. 
  int nelim = to_elim.size();
  // Q:= see above
  // iQ:= inverse of Q
  FastMat2 Q(2,nelim,nelim),iQ(2,nelim,nelim);
  // B[j] is a pointer to row j-1
  vector<row_t *> B(nelim);

  // For each edof to be eliminated we assign an index `elim' 0<=j<nelim
  // Costruct map from edofs to index in the set to be eliminated
  map<int,int> edof2elim;
  vector<int> elim2edof(nelim);
  int elim=0;
  for (set<int>::iterator q=to_elim.begin(); 
       q!=to_elim.end(); q++) {
    edof2elim[*q] = ++elim;
    elim2edof[elim-1] = *q;
  }

  // get rows for the edofs to be elimineated and
  // build the Q and B matrices. 
  Q.set(0.);
  for (int r=1; r<=nelim; r++) {
    // recall toe free after...
    B[r-1] = new row_t;
    row_t row_o;		// the original row
    int edoff,node, field;
    // get row
    edoff = elim2edof[r-1];
    edofi(edoff,ndof,node,field);
    get_row(node,field,row_o); 
    // Iterates on the original row. Columns in `to_elim'
    // are passed to the Q matrix, and the others to their respective
    // rows in `B'
    for (row_t::iterator q=row_o.begin(); 
	 q!=row_o.end(); q++) {
      // get column number
      int edof_c = q->first;
      if (to_elim.find(edof_c)!=to_elim.end()) {
	// column number has to be eliminated also ....
	Q.setel(q->second,r,edof2elim[edof_c]);
      } else {
	// ... or not
	B[r-1]->insert(*q);
      }
    }
    // Add the diagonal term
    Q.addel(-1.,r,r);
    // The row shouldn't be a `regular' row
    assert(fabs(Q.get(r,r))>tol);
  }

#ifdef DEBUG_ELIM
  // Print Q and B matrices
  Q.print("Q:");
  for (int k=0; k<nelim; k++) {
    printf("%d:",k+1);
    print(*B[k]);
    printf("\n");
  }
#endif

  // Compute inverse of Q
  iQ.inv(Q);

  // I don't know if this is needed really. 
  // We look for dependencies, i.e. the edofs in the
  // rows depend already on other edof's. 
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

  // This is not needed but may be is better to do. 
  Q.clear();
  iQ.clear();
  // Free memory in B
  for (int k=0; k<nelim; k++) delete B[k];
}
