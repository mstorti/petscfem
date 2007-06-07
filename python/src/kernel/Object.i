// -*- c++ -*-
// $Id$

PF4PY_NAMESPACE_BEGIN
%feature("ref")   Object "$this->incref();"
%feature("unref") Object "$this->decref();"
%extend Object { 
  static int __allrefs__()  { return PETScFEM::Object::allrefs() ; }
}
%extend Object { 
  int __refcount__() { return $self->getref(); }
}
PF4PY_NAMESPACE_END

PF4PY_NAMESPACE_BEGIN
%ignore Object::~Object;
%ignore Object::getref;
%ignore Object::incref;
%ignore Object::decref;
PF4PY_NAMESPACE_END

%include "Object.h"
