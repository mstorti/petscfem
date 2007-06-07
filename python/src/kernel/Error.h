// $Id$

#ifndef PF4PY_ERROR_H
#define PF4PY_ERROR_H

#include <string>
#include <sstream>
#include <petscerror.h>
#include "namespace.h"

PF4PY_NAMESPACE_BEGIN

class Error
{
 protected:
  std::string message;
 public:
  Error(const std::string& msg) : message(msg) { }
  Error(PetscErrorCode ierr, const std::string& msg)
    : message(msg)
  { 
    if (ierr) {
      const char* text = NULL;
      char*       spec = NULL;
      PetscErrorMessage(ierr, &text, &spec);
      std::stringstream sstream;
      if (text and text[0] != 0) 
	{ sstream << std::endl << text; }
      if (spec and spec[0] != 0) 
	{ sstream << std::endl << spec; }
      this->message += sstream.str();
    }
  }
  operator std::string() const { return message; }
  operator const char*() const { return message.c_str(); }
};


PF4PY_NAMESPACE_END

#endif // PF4PY_ERROR_H

// Local Variables:
// mode: C++
// End:
