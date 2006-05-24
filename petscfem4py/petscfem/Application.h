// $Id: Application.h,v 1.1.2.1 2006/05/24 21:02:33 dalcinl Exp $

#ifndef PYPF_APPLICATION_H
#define PYPF_APPLICATION_H

#include "petscfem4py.h"
#include "Object.h"
#include "Domain.h"

PYPF_NAMESPACE_BEGIN

class Application
  : public Object
{
private:
  Application();

protected:
  Domain* domain;

public:
  ~Application();
  Application(const Application&);
  Application(Domain& domain);

  Domain& getDomain() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_APPLICATION_H

// Local Variables:
// mode: C++
// End:
