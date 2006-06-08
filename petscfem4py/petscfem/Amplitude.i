// -*- c++ -*-
// $Id: Amplitude.i,v 1.1.2.3 2006/06/08 15:44:52 dalcinl Exp $

PYPF_NAMESPACE_BEGIN

%feature("ref")   Amplitude "$this->incref();"
%feature("unref") Amplitude "$this->decref();"

%director Amplitude;
%director Amplitude::operator();

%director AmpConstant;
%director AmpConstant::operator();

%director AmpTemporal;
%director AmpTemporal::operator();

%director AmpNodal;
%director AmpNodal::operator();


%director AmpGeneral;
%director AmpGeneral::operator();


PYPF_NAMESPACE_END


%include "Amplitude.h"
