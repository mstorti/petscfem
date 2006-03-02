// -*- c++ -*-

#ifndef PYPF_NAVIERSTOKES_H
#define PYPF_NAVIERSTOKES_H


#include <string>
#include "petscfem4py.h"
#include "Problem.h"


PYPF_NAMESPACE_BEGIN

class NavierStokes: public Problem
{

public:
  NavierStokes(int nnod, int ndim, int ndof);
  ~NavierStokes();
  
  void computeResidual();
  void computeJacobian();
};


PYPF_NAMESPACE_END

#endif // PYPF_NAVIERSTOKES_H
