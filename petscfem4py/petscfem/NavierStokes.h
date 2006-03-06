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
  NavierStokes();
  NavierStokes(int nnod, int ndim);
  ~NavierStokes();
  
  void computeResidual();
  void computeJacobian();
};


PYPF_NAMESPACE_END

#endif // PYPF_NAVIERSTOKES_H
