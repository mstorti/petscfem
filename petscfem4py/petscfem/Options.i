// -*- c++ -*-
// $Id: Options.i,v 1.1.2.3 2006/08/22 22:10:43 dalcinl Exp $

namespace std { }
%template() std::pair <std::string, std::string>;
%template() std::map  <std::string, std::string>;

PYPF_NAMESPACE_BEGIN

%naturalvar   Options;
%implicitconv Options;

%ignore Options::operator=;
%ignore Options::empty;
%ignore Options::init;

%rename(__len__)      Options::size;
%rename(__contains__) Options::has;
%rename(__getitem__)  Options::get;
%rename(__setitem__)  Options::set;
%rename(__delitem__)  Options::del;

%extend Options {
  const std::map<std::string,std::string>& todict() { return *self; }
}

PYPF_NAMESPACE_END


%include "Options.h"


%init %{ 
  PYPF_NAMESPACE::Options::init();
%}
