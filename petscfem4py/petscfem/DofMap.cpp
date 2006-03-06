// $Id: DofMap.cpp,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

#include "DofMap.h"

#include <fem.h>
#include <dofmap.h>

PyPF::DofMap::~DofMap() 
{ }

PyPF::DofMap::DofMap()
{ }

PyPF::OptionTable*
PyPF::DofMap::get_opt_table(bool create)
{
  throw Error("no options for DofMap");
  return NULL;
}

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
  if ((*this)->id == NULL) throw Error("null id map");

  DofMap::Base* dofmap = *this;

  row_t row;
  for (int i=0; i<n; i++) {
    dofmap->get_row(node[i]+1, field[i]+1, row);
    if (row.size()!=1) throw Error("bad fixation for node/field combination");
    int keq = row.begin()->first;
#if 0
    dofmap->fixed.push_back(fixation_entry(value[i]));
    dofmap->fixed_dofs[keq] = dofmap->fixed.size()-1;
#else
    typedef fixation_entry fixentry;
    std::vector<fixentry>& fixed      = dofmap->fixed;
    std::map<int,int>&     fixed_dofs = dofmap->fixed_dofs;
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
  if ((*this)->id == NULL) throw Error("null id map");

  DofMap::Base* dofmap = *this;

  ::Constraint constraint;
  for (int i=0; i<n; i++) {
    constraint.add_entry(node[i]+1, field[i]+1, coef[i]);
  }
  dofmap->set_constraint(constraint);
}

void
PyPF::DofMap::getSize(int* nnod, int* ndof) 
{
  /* test */
  if ((*this)->id == NULL) throw Error("null id map");

  DofMap::Base* dofmap = *this;
  
  *nnod = dofmap->nnod;
  *ndof = dofmap->ndof;
}

void
PyPF::DofMap::view()
{
  /* test */
  DofMap::Base* dofmap = *this;
  if (dofmap == NULL || dofmap->id == NULL) return;

  dofmap->id->print();
}
