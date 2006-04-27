// -*- c++ -*-
// $Id: Options.i,v 1.1.2.1 2006/04/27 19:09:17 rodrigop Exp $


PYPF_NAMESPACE_BEGIN
%ignore Options;
%ignore Options::Options;
%ignore Options::operator=;
%ignore Options::del;
%ignore OPTIONS::GLOBAL;
%ignore OPTIONS::init;
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


%include "Options.h"


%init %{
  PYPF_NAMESPACE::OPTIONS::init();
%}
