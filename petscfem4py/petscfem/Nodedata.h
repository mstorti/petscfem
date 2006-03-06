// -*- c++ -*-

#ifndef PYPF_NODEDATA_H
#define PYPF_NODEDATA_H

#include <string>
#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

PYPF_CLASS(Nodedata)
{
  PYPF_CTOR_FROM_PTR(Nodedata)
  PYPF_OBJ_GETOPTTBL_DECL

 public:
  
  ~Nodedata();
  Nodedata();

  void getData(int* nnod, int* ndim, double* xyz[]);
  void setData(int  nnod, int  ndim, double  xyz[]);
  void getSize(int* nnod, int* ndim);

  void view();
  
  friend class Mesh;
};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
