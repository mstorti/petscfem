// -*- c++ -*-

#ifndef PYPF_DOFMAP_H
#define PYPF_DOFMAP_H

#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(DofMap)
{

  PYPF_CONSTRUCTOR(DofMap)

 public:
  
  DofMap(int nnod, int ndof);
  ~DofMap();

  void addFixations(int n, int node[], int field[], double value[]);
  
  void addConstraints(int n, int node[], int field[], double coef[]);

  void print();
};


PYPF_NAMESPACE_END

#endif // PYPF_DOFMAP_H
