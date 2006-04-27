// -*- c++ -*-

PYPF_NAMESPACE_BEGIN

%director Amplitude;
%director Amplitude::operator();

%feature("ref")   Amplitude "$this->incref();"
%feature("unref") Amplitude "$this->decref();"

PYPF_NAMESPACE_END


%include "Amplitude.h"
