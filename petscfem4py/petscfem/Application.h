// $Id: Application.h,v 1.1.2.8 2006/11/29 22:35:09 dalcinl Exp $

#ifndef PYPF_APPLICATION_H
#define PYPF_APPLICATION_H

#include "petscfem4py.h"
#include "Object.h"
#include "Domain.h"
#include "ArgList.h"

PYPF_NAMESPACE_BEGIN

class Application
  : public Object
{
protected:
  Application();
  
protected:
  Domain*    domain;
  VecScatter scatter;
  Vec        state;

protected:
  static void assemble(const Application&, const ArgList&);

public:
  ~Application();
  Application(const Application& application);
  Application(Domain& domain);

public:

  Domain& getDomain() const;

  virtual void allocateState(Vec& x, const std::string&vec_type="") const;
  virtual void allocateResidual(Vec& r, const std::string&vec_type="") const;
  virtual void allocateJacobian(Mat& J, const std::string& mat_type="") const;
  virtual void allocateSolution(Vec& u) const;

  virtual void buildState(Vec solution, Vec state);
  virtual void buildSolution(double time, Vec state, Vec solution);
  
};

PYPF_NAMESPACE_END

#endif // PYPF_APPLICATION_H

// Local Variables:
// mode: C++
// End:
