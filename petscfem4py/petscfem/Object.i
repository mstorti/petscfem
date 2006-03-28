// -*- c++ -*-


PYPF_NAMESPACE_BEGIN

%ignore Object::~Object;
%ignore Object::get_opt_table;
%ignore Object::getref;
%ignore Object::incref;
%ignore Object::decref;
  
%feature("ref")   Object "$this->incref();"
%feature("unref") Object "$this->decref();"

// %feature("pythonprepend") Object::setOption
// %{args = (args[0], args[1], str(args[2]))%}

// %feature("pythonprepend") Object::setOptions
// %{args = (args[0], dict((k, str(v)) for k, v in args[1].iteritems()))%}

PYPF_NAMESPACE_END

%template() std::pair<std::string, std::string>;
%template() std::map<std::string, std::string>;

%init {
  TextHashTable* global_options = PyPF_NAMESPACE::Options::GLOBAL;
  global_options->register_name("global_options");
  global_options->set_as_global();
  GLOBAL_OPTIONS = global_options;
}

%include "Object.h"
