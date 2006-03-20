// -*- c++ -*-

#ifndef PYPF_NAVIERSTOKES_H
#define PYPF_NAVIERSTOKES_H


#include <string>
#include "petscfem4py.h"
#include "Problem.h"


PYPF_NAMESPACE_BEGIN

class NavierStokes: public Problem
{

private:
  //NavierStokes();
  NavierStokes(const NavierStokes &);

protected:
  void* nsargs;

public:
  ~NavierStokes();
  NavierStokes();
  NavierStokes(int nnod, int ndim);
  
  // steady problems
  void assemble(Vec x, double t,
		Vec r, Mat J);

  // transient problems
  void assemble(Vec x_0, double t_0,
		Vec x_1, double t_1,
		Vec r, Mat J,
		double alpha=1.0);

};


PYPF_NAMESPACE_END

#endif // PYPF_NAVIERSTOKES_H
