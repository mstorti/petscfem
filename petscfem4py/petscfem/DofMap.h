// -*- c++ -*-
// $Id: DofMap.h,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

#ifndef PYPF_DOFMAP_H
#define PYPF_DOFMAP_H

#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN


PYPF_CLASS(DofMap)
{
  PYPF_CTOR_FROM_PTR(DofMap)
  PYPF_OBJ_GETOPTTBL_DECL

 public:
  
  DofMap();
  ~DofMap();

  void addFixations(int n, int node[], int field[], double value[]);
  void addConstraints(int n, int node[], int field[], double coef[]);
  void getSize(int* nnod, int* ndof);


  void view();
};


PYPF_NAMESPACE_END

#endif // PYPF_DOFMAP_H
