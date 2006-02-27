// -*- c++ -*-

%{
#include "DofMap.h"
%}


ARRAY_TRIAD(int n, int node[], int field[], double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_TRIAD(int n, int node[], int field[], double coef[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

%include "DofMap.h"
