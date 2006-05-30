// $Id: ArgList.h,v 1.1.2.1 2006/05/30 20:13:05 dalcinl Exp $

#ifndef PYPF_ARGLIST_H
#define PYPF_ARGLIST_H

#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

struct ArgList : public SmartPtr< ::ArgList >
{
  ~ArgList();
  ArgList();
  virtual const char* job()  const = 0;
  virtual const Time* time() const = 0;
};

PYPF_NAMESPACE_END

#endif // PYPF_ARGLIST_H

// Local Variables:
// mode: C++
// End:
