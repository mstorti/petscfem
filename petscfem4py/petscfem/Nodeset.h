// -*- c++ -*-
// $Id: Nodeset.h,v 1.1.2.1 2006/03/30 15:40:05 rodrigop Exp $

#ifndef PYPF_NODESET_H
#define PYPF_NODESET_H

#include <vector>
#include "petscfem4py.h"
#include "Object.h"

PYPF_NAMESPACE_BEGIN

class Nodeset : SMARTPTR(Nodeset)
  public Object
{
  friend class Mesh;
  friend class DofMap;
  friend class Problem;

#if !defined(SWIG)
 public:
  Nodeset(Nodeset::Base*);
#endif

 protected:
  int                 nnod;
  int                 ndim;
  std::vector<double> nodedata;

 public:
  ~Nodeset();
  Nodeset();
  Nodeset(const Nodeset&);
  Nodeset(int nnod, int ndim, const double xyz[]);

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

#endif // PYPF_NODESET_H
