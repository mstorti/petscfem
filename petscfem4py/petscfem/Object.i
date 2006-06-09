// -*- c++ -*-
// $Id: Object.i,v 1.1.2.6 2006/06/09 14:23:32 dalcinl Exp $

%include Options.i
%include Comm.i

namespace std { }
%template() std::pair <std::string, std::string>;
%template() std::map  <std::string, std::string>;

PYPF_NAMESPACE_BEGIN
%feature("ref")   Object "$this->incref();"
%feature("unref") Object "$this->decref();"
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%ignore Object::~Object;
%ignore Object::getref;
%ignore Object::incref;
%ignore Object::decref;
PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Object {
  Comm getComm() { return self->getComm(); }
  %ignore getComm;
}
PYPF_NAMESPACE_END

PYPF_NAMESPACE_BEGIN
%extend Object {
  %pythoncode {
  comm    = property(getComm,    setComm,    doc="MPI communicator")
  options = property(getOptions, setOptions, doc="object options")
  }
}
PYPF_NAMESPACE_END


%include "Object.h"
