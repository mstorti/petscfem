//__INSERT_LICENSE__
// $Id: generror.cpp,v 1.2 2003/12/06 15:11:09 mstorti Exp $

#include <string>

using namespace std;

#include <src/generror.h>
#include <src/autostr.h>

GenericError::
GenericError(char *s,va_list ap) {
  AutoString as;
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}

GenericError::
GenericError(char *s,...) {
  AutoString as;
  va_list ap;

  va_start(ap,s);
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}
