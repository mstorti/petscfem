// $Id: Error.h,v 1.1.2.3 2006/04/27 19:09:17 rodrigop Exp $

#ifndef PYPF_ERROR_H
#define PYPF_ERROR_H

#include <string>
#include "namespace.h"

PYPF_NAMESPACE_BEGIN

class Error
{
 protected:
  std::string message;
 public:
  Error(const std::string& msg) : message(msg) { }
  operator std::string() const { return message; }
  operator const char*() const { return message.c_str(); }
};


PYPF_NAMESPACE_END

#endif // PYPF_NODEDATA_H

// Local Variables:
// mode: C++
// End:
