// -*- c++ -*-
// $Id: Amplitude.i,v 1.1.2.2 2006/05/24 20:59:59 dalcinl Exp $

PYPF_NAMESPACE_BEGIN

%director Amplitude;
%director Amplitude::operator();

%feature("ref")   Amplitude "$this->incref();"
%feature("unref") Amplitude "$this->decref();"

PYPF_NAMESPACE_END


%include "Amplitude.h"
