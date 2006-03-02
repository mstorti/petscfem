// $Id: DofMap.cpp,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

#include "DofMap.h"

#include <fem.h>
#include <dofmap.h>

PyPF::DofMap::DofMap() { }
PyPF::DofMap::~DofMap() { }

// PyPF::DofMap::DofMap(int nnod, int ndof)
//   : Ptr(new ::DofMap)
// { 
//   (*this)->nnod   = nnod;
//   (*this)->ndof   = ndof;
//   (*this)->id     = new idmap(nnod*ndof, NULL_MAP);
//   (*this)->tpwgts = new float[1];
  
//   for (int i=1; i<=nnod; i++)
//     for (int j=1; j<=ndof; j++) {
//       int idx = (*this)->edof(i,j);
//       (*this)->id->set_elem(idx,idx,1.);
//     }
  
// }


void
PyPF::DofMap::addFixations(int n, int node[], int field[], double value[])
{
  /* test */
  if ((*this) == NULL)  throw Error("null DofMap");
  if ((*this)->id == NULL) throw Error("null id map");

  row_t row;
  for (int i=0; i<n; i++) {
    (*this)->get_row(node[i]+1, field[i]+1, row);
    if (row.size()!=1) throw Error("bad fixation for node/field combination");
    int keq = row.begin()->first;
#if 0
    (*this)->fixed.push_back(fixation_entry(value[i]));
    (*this)->fixed_dofs[keq] = (*this)->fixed.size()-1;
#else
    typedef fixation_entry fixentry;
    std::vector<fixentry>& fixed      = (*this)->fixed;
    std::map<int,int>&     fixed_dofs = (*this)->fixed_dofs;
    if (fixed_dofs.find(keq) == fixed_dofs.end()) {
      fixed.push_back(fixentry(value[i]));
      fixed_dofs[keq] = fixed.size()-1;
    } else {
      fixed[fixed_dofs[keq]] = fixentry(value[i]);
    }
#endif
  }
}

void
PyPF::DofMap::addConstraints(int n, int node[], int field[], double coef[])
{
  /* test */
  /* test */
  if ((*this) == NULL)  throw Error("null DofMap");
  if ((*this)->id == NULL) throw Error("null id map");

  ::Constraint constraint;
  for (int i=0; i<n; i++) {
    constraint.add_entry(node[i]+1, field[i]+1, coef[i]);
  }
  (*this)->set_constraint(constraint);
}

void
PyPF::DofMap::getSize(int* nnod, int* ndof) 
{
  /* test */
  if ((*this) == NULL)  throw Error("null DofMap");
  if ((*this)->id == NULL) throw Error("null id map");
  
  *nnod = (*this)->nnod;
  *ndof = (*this)->ndof;
}

void
PyPF::DofMap::view()
{
  if ((*this) == NULL) return;
  if ((*this)->id == NULL) return;

  (*this)->id->print();
}
