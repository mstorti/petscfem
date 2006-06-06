// $Id: NavierStokes.h,v 1.1.2.8 2006/06/06 15:46:52 dalcinl Exp $

#ifndef PYPF_NAVIERSTOKES_H
#define PYPF_NAVIERSTOKES_H


#include "petscfem4py.h"
#include "Problem.h"


PYPF_NAMESPACE_BEGIN

class ArgsNS;

class NavierStokes
  : public Problem
{

private:
  NavierStokes();
  NavierStokes(const NavierStokes&);

protected:
  ArgsNS* args;

public:
  ~NavierStokes();
  NavierStokes(Mesh& mesh, Dofset& dofset);
  NavierStokes(Nodeset& nodeset,
	       const std::vector<Elemset*>& elemsets,
	       Dofset& dofset);


  // steady problems
  void assemble(Vec x, double t,
		Vec r, Mat J);

  // transient problems
  void assemble(Vec x0, double t0,
		Vec x1, double t1,
		Vec r, Mat J,
		double alpha=1.0, bool steady=false);

//   // transient problems
//   void assembleResidual(Vec x0, double t0,
// 			Vec x1, double t1,
// 			Vec r,   
// 			double alpha=1.0,bool steady=false);


};

PYPF_NAMESPACE_END

#endif // PYPF_NAVIERSTOKES_H

// Local Variables:
// mode: C++
// End:
