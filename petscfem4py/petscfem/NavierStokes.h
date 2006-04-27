// $Id: NavierStokes.h,v 1.1.2.7 2006/04/27 19:09:17 rodrigop Exp $

#ifndef PYPF_NAVIERSTOKES_H
#define PYPF_NAVIERSTOKES_H


#include <petscvec.h>
#include <petscmat.h>
#include "petscfem4py.h"
#include "Problem.h"


PYPF_NAMESPACE_BEGIN

class NavierStokes
  : public Problem
{

private:
  NavierStokes();
  NavierStokes(const NavierStokes&);

protected:
  void* nsargs;

public:
  ~NavierStokes();
  NavierStokes(Mesh& mesh, DofMap& dofmap);
  NavierStokes(Nodeset& nodeset,
	       const std::vector<Elemset*>& elemsets,
	       Dofset& dofset);


  // steady problems
  void assemble(Vec x, double t,
		Vec r, Mat J);

  // transient problems
  void assemble(Vec x_0, double t_0,
		Vec x_1, double t_1,
		Vec r, Mat J,
		double alpha=1.0);

  // transient problems
  void assembleResidual(Vec x_0, double t_0,
			Vec x_1, double t_1,
			Vec r,   
			double alpha=1.0);


};


PYPF_NAMESPACE_END

#endif // PYPF_NAVIERSTOKES_H

// Local Variables:
// mode: C++
// End:
