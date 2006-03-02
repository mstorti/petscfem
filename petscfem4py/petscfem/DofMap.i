// -*- c++ -*-
// $Id: DofMap.i,v 1.1.2.2 2006/03/02 21:37:12 rodrigop Exp $

ARRAY_TRIAD(int n, int node[], int field[], double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

ARRAY_TRIAD(int n, int node[], int field[], double coef[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

%apply int* OUTPUT {int* nnod, int* ndof};

%include "DofMap.h"
