// -*- c++ -*-
// $Id: NavierStokes.h,v 1.1.2.6 2006/03/30 15:40:05 rodrigop Exp $

#ifndef PYPF_NAVIERSTOKES_H
#define PYPF_NAVIERSTOKES_H


#include <string>
#include "petscfem4py.h"
#include "Problem.h"


PYPF_NAMESPACE_BEGIN

class NavierStokes:  public Problem
{

private:
  NavierStokes();
  NavierStokes(const NavierStokes&);

protected:
  void* nsargs;

public:
  ~NavierStokes();
  NavierStokes(Mesh*, DofMap*);


  //static NavierStokes* fromFile(const std::string& filename);

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
