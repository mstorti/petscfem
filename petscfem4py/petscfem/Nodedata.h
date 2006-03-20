// -*- c++ -*-

#ifndef PYPF_NODEDATA_H
#define PYPF_NODEDATA_H

#include <string>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Nodedata : SMARTPTR(Nodedata)
  public Object
{
  friend class Mesh;

protected: 
  OptionTable* get_opt_table() const;

#if !defined(SWIG)
 public:
  Nodedata(Nodedata::Base*);
#endif

 public:
  ~Nodedata();
  Nodedata();
  Nodedata(const Nodedata&);
  Nodedata(int nnod, int ndim);
  Nodedata(int nnod, int ndim, double xyz[]);

  void getSize(int* nnod, int* ndim) const;
  void getData(int* nnod, int* ndim, double* xyz[]) const;
  void setData(int  nnod, int  ndim, double  xyz[]);
  void delData();

  void setUp();
  void clear();
  void view() const;
  
};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
