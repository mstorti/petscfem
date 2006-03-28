// -*- c++ -*-

#ifndef PYPF_NODEDATA_H
#define PYPF_NODEDATA_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Nodedata : SMARTPTR(Nodedata)
  public Object
{
  friend class Mesh;
  friend class DofMap;

#if !defined(SWIG)
 public:
  Nodedata(Nodedata::Base*);
#endif

 protected:
  int nnod, ndim;
  std::vector<double> nodedata;

 public:
  ~Nodedata();
  Nodedata();
  Nodedata(const Nodedata&);
  Nodedata(int nnod, int ndim, const double xyz[]);

  void getSize(int* nnod, int* ndim) const;
  void getData(int* nnod, int* ndim, double* xyz[]) const;
  void setData(int  nnod, int  ndim, const double  xyz[]);

  typedef std::vector<double> Node;
  Node getNode(int i) const;
  void setNode(int i, const Node& node);

  void sync(int root = 0);

  void setUp();
  void clear();
  void view() const;

};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H
