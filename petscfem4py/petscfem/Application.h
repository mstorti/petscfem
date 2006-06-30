// $Id: Application.h,v 1.1.2.7 2006/06/30 18:33:33 dalcinl Exp $

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

  void buildState(Vec solution, Vec state);
  void buildSolution(double time, Vec state, Vec solution);
  
};

PYPF_NAMESPACE_END

#endif // PYPF_APPLICATION_H

// Local Variables:
// mode: C++
// End:
