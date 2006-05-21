// -*- c++ -*-
// $Id: Options.i,v 1.1.2.2 2006/05/21 03:48:11 dalcinl Exp $


PYPF_NAMESPACE_BEGIN
%ignore Options;
%ignore Options::Options;
%ignore Options::operator=;
%ignore Options::del;
PYPF_NAMESPACE_END


// PYPF_NAMESPACE_BEGIN
// %extend OPTIONS {
//   %pythoncode {
//     def __contains__(self, key):
//         return self.hasOption(key)
//     def __getitem__(self, key):
//         return self.getOption(key)
//     def __setitem__(self, key, value):
//         return self.setOption(key, value)
//     def __delitem__(self, key):
//         return self.delOption(key)
//   }
// }
// PYPF_NAMESPACE_END


PYPF_NAMESPACE_BEGIN
%ignore OPTIONS::GLOBAL;
%ignore OPTIONS::init;
PYPF_NAMESPACE_END


%include "Options.h"


%init %{
  PYPF_NAMESPACE::OPTIONS::init();
%}
