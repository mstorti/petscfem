// $Id: Application.h,v 1.1.2.4 2006/05/30 20:13:55 dalcinl Exp $

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
  Domain* domain;

public:
  ~Application();
  Application(const Application& application);
  Application(Domain& domain);

public:
  Domain& getDomain() const;

  void getNodalValues(const double state[], double time,
		      double values[]) const;

  void getNodalValues(const double state[], double time,
		      int nn, const int nodes[],
		      int nf, const int fields[],
		      double values[]) const;

protected:
  static void assemble(const Application&, const ArgList&);
  
};

PYPF_NAMESPACE_END

#endif // PYPF_APPLICATION_H

// Local Variables:
// mode: C++
// End:
