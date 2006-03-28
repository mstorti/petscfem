// -*- c++ -*-
// $Id: DofMap.i,v 1.1.2.4 2006/03/28 22:13:25 rodrigop Exp $

ARRAY_TRIAD(int n, int node[], int field[], double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_TRIAD(int n, int node[], int field[], double coeff[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_1D_NEW(int* rsize, int* ranges[], PyPF_INT)

PYPF_NAMESPACE_BEGIN

%pythonprepend DofMap::addFixations
%{args = (args[0],args[1:])%}

%pythonprepend DofMap::addConstraints 
%{args = (args[0],args[1:])%}

PYPF_NAMESPACE_END

%apply int* OUTPUT {int* local,  int* global};
%apply int* OUTPUT {int* start,  int* end};

%apply int* OUTPUT {int* nnod, int* ndof};
%apply int* OUTPUT {int* neq_fix, int* neq_tot};

%include "DofMap.h"

%clear (int n, int node[], int field[], double value[]);
%clear (int n, int node[], int field[], double coeff[]);
%clear (int* rsize, int* ranges[]);

%clear int* local, int* global;
%clear int* start, int* end;

%clear int* nnod, int* ndof;
%clear int* neq_fix, int* neq_tot;
