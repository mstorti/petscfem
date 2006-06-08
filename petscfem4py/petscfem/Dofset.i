// -*- c++ -*-

%include Object.i
%include Amplitude.i


%typemap(doc,name="node, field, value",type="int[], int[], double[]")
  (int n, const int node[],const int field[],const double value[])
  "node, field, value: int[], int[], double[] array values";
ARRAY_TRIAD(int n, const int node[], const int field[], const double value[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

%typemap(doc,name="node, field, coeff",type="int[], int[], double[]")
  (int n, const int node[],const int field[],const double coeff[])
  "node, field, coeff: int[], int[], double[] array values";
ARRAY_TRIAD(int n, const int node[], const int field[], const double coeff[],
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_INT,
	    ARRAY_INPUT, PyPF_FLOAT)

%apply (int*, int*) { (int* nnod, int* ndof) };

PYPF_NAMESPACE_BEGIN
%pythonprepend Dofset::addFixations
%{args = (args[0],args[1:4]) + args[4:]%}
%pythonprepend Dofset::addConstraints
%{args = (args[0],args[1:4])%}
PYPF_NAMESPACE_END

%include "Dofset.h"

%clear (int* nnod, int* ndof);
%clear (int n, const int node[], const int field[], const double value[]);
%clear (int n, const int node[], const int field[], const double coeff[]);
