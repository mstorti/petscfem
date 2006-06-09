// -*- c++ -*-
// $Id: Amplitude.i,v 1.1.2.4 2006/06/09 14:18:25 dalcinl Exp $

PYPF_NAMESPACE_BEGIN

%feature("ref")   Amplitude "$this->incref();"
%feature("unref") Amplitude "$this->decref();"

%director Amplitude;
%director Amplitude::operator();

PYPF_NAMESPACE_END


%include "Amplitude.h"
