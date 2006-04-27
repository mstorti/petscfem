// -*- c++ -*-

%include Object.i
%include Amplitude.i


PYPF_NAMESPACE_BEGIN

%pythonprepend Dofset::addFixations
%{args = len(args)==4 and (args[0],args[1:4], None) or (args[0],args[1:4], args[4])%}

%pythonprepend Dofset::addConstraints
%{args = (args[0],args[1:4])%}

PYPF_NAMESPACE_END

ARRAY_TRIAD(int n, const int node[], const int field[], const double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_TRIAD(int n, const int node[], const int field[], const double coeff[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)


%include "Dofset.h"

%clear (int n, const int node[], const int field[], const double value[]);
%clear (int n, const int node[], const int field[], const double coeff[]);
