// -*- c++ -*-
// $Id: Object.i,v 1.1.2.4 2006/04/27 19:09:17 rodrigop Exp $

namespace std { }
%template() std::pair<std::string, std::string>;
%template() std::map<std::string, std::string>;

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
  
// PYPF_NAMESPACE_BEGIN
// %feature("pythonprepend") Object::setOption
// %{args = (args[0], args[1], str(args[2]))%}

// %feature("pythonprepend") Object::setOptions
// %{args = (args[0], dict((k, str(v)) for k, v in args[1].iteritems()))%}
// PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%extend Object {
  %pythoncode {
  comm    = property(getComm,    setComm,    doc="MPI communicator")
  options = property(getOptions, setOptions, doc="object options")
  }
}
PYPF_NAMESPACE_END


%include "Object.h"
