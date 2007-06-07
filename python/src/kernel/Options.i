// -*- c++ -*-
// $Id$

namespace std { }
%template() std::pair <std::string, std::string>;
%template() std::map  <std::string, std::string>;

PF4PY_NAMESPACE_BEGIN

%ignore Options::initGlobals;
%ignore Options::empty;

%rename(__len__)      Options::size;
%rename(__contains__) Options::has;
%rename(__getitem__)  Options::get;
%rename(__setitem__)  Options::set;
%rename(__delitem__)  Options::del;

%extend Options {
  const std::map<std::string,std::string>& todict() { return *self; }
}

PF4PY_NAMESPACE_END

%init %{ 
  PF4PY_NAMESPACE::Options::initGlobals();
%}

PF4PY_NAMESPACE_BEGIN
%implicitconv Options;
%typemap(freearg,noblock=1,implicitconv=1) 
  const Options& options
{ if (SWIG_IsNewObj(res$argnum)) %decref($1); }
PF4PY_NAMESPACE_END


%include "Options.h"
