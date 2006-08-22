// -*- c++ -*-
// $Id: Object.i,v 1.1.2.7 2006/08/22 22:10:43 dalcinl Exp $

%include Options.i
%include Comm.i

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
  int __refcount__() { return self->getref(); }
}
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
  comm = property(getComm, setComm, doc="MPI communicator")
  }
}
PYPF_NAMESPACE_END


%include "Object.h"
