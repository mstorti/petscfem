//__INSERT_LICENSE__
//$Id: dofmap2.cpp,v 1.14 2004/09/25 23:11:39 mstorti Exp $

#include <cassert>
#include <deque>

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
#include "fastmat2.h"

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Amplitude::eval(const TimeData *time_data) {
  PETSCFEM_ERROR0("Not implemented eval(const TimeData *time_data)");
  return 0.;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Eval the amplitude of the function at this time (needs node and field)
double Amplitude::eval(const TimeData *time_data,int node,int field) {
  PETSCFEM_ERROR0("Not implemented eval(const TimeData *time_data,"
		  "int node,int field)");
  return 0.;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Amplitude::needs_dof_field_q() { return 0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::nodf(int edof, int &node,int &field) const {
  field = (edof-1) % ndof +1;
  node = (edof-field)/ndof +1;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Dofmap::set_constraint"
int Dofmap::set_constraint(const Constraint &constraint) {
  //#define DEBUG_ELIM
  double tol = 1e-10;
  int length=constraint.size();
  row_t row,row0;
  
  // Create row from constraint and eliminate all already eliminated
  // variables. This is done by simply making a linear combination
  // of the corresponding rows in the current `idmap'
  for (int k=0; k<length; k++) {
    const constraint_entry *it = &(constraint[k]);
    int edoff = edof(it->node,it->field);
    // Look if row is regular (it doesn't point to other edof's)
    row.clear();
#ifdef DEBUG_ELIM
    printf("row item: node %d, dof %d, edof %d, coef %f\n",
	   it->node,it->field,edof(it->node,it->field),
	   it->coef);
#endif    
    get_row(it->node,it->field,row); 
    axpy(row0,it->coef,row);
  }

  // Look for the max abs.val. coefficient in the restriction
  int set_flag=0,all_non_null_flag=1;
  double cc, cmax; 
  row_t::iterator q, qmax, qe = row0.end();
  for (q=row0.begin(); q!=qe; q++) {
    if (fixed_dofs.find(q->first)!=fixed_dofs.end()) continue;
    all_non_null_flag=0;
    cc = fabs(q->second);
    if (!set_flag || cc>cmax) {
      qmax = q;
      cmax = cc;
      set_flag = 1;
    }
  }
  // If the row is null, then it means that the restriction is
  // linerly dependent on previous restrictions and do nothing
  if (all_non_null_flag || cmax<tol) return 1;
  if (!set_flag) return 2;
  // This is the `edof' to be eliminated
  double coef0 = qmax->second;
  int edof0 = qmax->first;
  row0.erase(qmax);

  // The column corresponding to this `edof'
  row_t col;
  id->get_col(edof0,col);
  
  // For each row that refers to `edof0' eliminate the `edof0' entry.
  qe = col.end();
  for (q=col.begin(); q!=qe; q++) {
    int node,field;
    nodf(q->first,node,field);
    row.clear();
    get_row(node,field,row);
    row_t::iterator r = row.find(edof0);
    double coef = r->second;
    row.erase(r);
    axpy(row,-coef/coef0,row0);
    row_set(node,field,row);
  }

  // Now, this is the dirty part... Fixations may have been entered as
  // an entry in the `fixed_dofs' array. That means that if we have
  // made a constraint that is a mere redirection of the unkown \phi_i -> \phi_j
  // then we have to meve the corresponding fixation. 
  map<int,int>::iterator w = fixed_dofs.find(edof0);
  if (w!=fixed_dofs.end()) {
    // `edof0' had a fixation, move...
    int node,field;
    nodf(edof0,node,field);
    get_row(node,field,row);
    // Check that it is really a simple mapping (only one entry in row with coef = 1.)
    PETSCFEM_ASSERT(row.size()==1,
		    "Constraint imposed on bad node/field combination\n"
		    "node %d, field %d\n",node,field);  
    row_t::iterator q = row.begin();
    assert(q->second==1.);
    // This is the new `edof'
    int edof = q->first;
    // printf("Moving edof %d -> %d!!\n",edof0,edof);
    // Move fixation
    int indx = w->second;
    fixed_dofs.erase(w);
    fixed_dofs[edof] = indx;
  }
  return 0;
}
