// $Id: ArgList.h,v 1.1.2.2 2006/06/08 16:14:21 dalcinl Exp $

#ifndef PYPF_ARGLIST_H
#define PYPF_ARGLIST_H

#include "namespace.h"
#include "forward.h"
#include "SmartPtr.h"

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
